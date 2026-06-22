"""Energy extraction from solved worksheets (ph2 energy).

Single-file commands print one CSV line `E,E_ref,dE`: the estimate,
the trial-state energy E_ref = <psi|H|psi> computed directly from
the input, and dE = E - E_ref.  E_ref is nan when /state_prep is
absent.  The physics is ported verbatim from the retired
examples/scripts/energy_{fft,mc,ref}.py (validated against exact
diagonalisation of the bundled fixture) and extended to
coeff_matrix worksheets through the Slater-Condon oracle in
stprep.py.

fft -- the Trotter drivers record values[s] = <psi|U^(s+1)|psi> at
uniform spacing dt = delta, so each eigenstate contributes a tone
exp(i*lambda*delta*s); a peak at FFT frequency f maps to

    E = tau * f / (delta * normalization) + offset

np.fft.fftfreq returns signed frequencies, so a negative-eigenvalue
ground state lands at negative f -- no sign flip.  The DC bin
carries the off-resonant baseline and is excluded; the dominant
remaining peak is the largest-overlap eigenstate.  Resolution is
one bin = tau/(steps*delta*norm).

mc -- the randomised drivers reset the state every sample, so
values are independent draws of the overlap at one effective time
T (qdrift: depth*asin(step_size); cmpsit: steps*angle_det); the
channel reproduces exp(i*T*H_norm) in expectation, so

    E = arg(mean(values)) / (T * normalization) + offset

The estimate is reliable while the mean stays coherent; a warning
is emitted when |mean| < 0.1.

rpe -- robust phase estimation on the dyadic subsequence
values[2^i - 1], i = 0..J = floor(log2(len)): each level refines
the per-step phase theta = lambda*delta with the 2^-i ambiguity
resolved against the previous estimate, and

    E0 = sqrt(1-x^2)/x * tan(x*thetaJ) / norm / delta,  x = 2^-J

Resurrected from the retired trott_rpe.py / qdrift_rpe.py with
three fixes: the candidate-distance generator is an explicit
argmin, result files are opened read-only, and the sweep reads
norm/offset from PREFIX-0.  rpe-qdrift consumes the multi-run
depth sweep PREFIX-0..J (depth doubling at fixed step size
x = 2^-J); tan undoes the per-gate asin exactly, so there is no
/delta in the sweep formula.
"""

import cmath
import math
import sys
from math import tau

import numpy as np

from . import Ph2Error, open_ro
from .stprep import expand_coeff_matrix

_PHASE = (1.0 + 0j, 1j, -1.0 + 0j, -1j)


def _masks(row):
    """(flip, phase) bit masks of a per-qubit Pauli byte row."""
    flip = phase = 0
    for j, b in enumerate(row):
        if b in (1, 2):
            flip |= 1 << j
        if b in (2, 3):
            phase |= 1 << j
    return flip, phase


def _amplitudes(f):
    """Trial-state {index: complex amplitude}, or None."""
    if "state_prep" not in f:
        return None
    sp = f["state_prep"]
    if "multidet" in sp:
        md = sp["multidet"]
        amp = {}
        for (re_, im), occ in zip(md["coeffs"][...],
                                  md["dets"][...]):
            idx = int(sum(int(b) << j for j, b in enumerate(occ)))
            amp[idx] = complex(re_, im)
        return amp
    if "coeff_matrix" not in sp:
        return None
    cm = sp["coeff_matrix"]
    ns, na, nb = (int(cm.attrs[k])
                  for k in ("n_sites", "n_alpha", "n_beta"))
    closed = bool(cm.attrs["closed_shell"])
    tapered = bool(cm.attrs["tapered"])
    blocks = []
    if "csf" in cm:
        for k in range(int(cm["csf"].attrs["n_components"])):
            blk = cm[f"csf/{k}"]
            cb = None if closed else blk["C_beta"][...]
            blocks.append((float(blk.attrs["coefficient"]),
                           blk["C_alpha"][...], cb))
    else:
        cb = None if closed else cm["C_beta"][...]
        blocks.append((1.0, cm["C_alpha"][...], cb))
    amp = {}
    for weight, ca, cb in blocks:
        part = expand_coeff_matrix(ns, na, nb, ca, cb, tapered,
                                   weight=weight)
        for i, v in part.items():
            amp[i] = amp.get(i, 0.0) + v
    return {i: complex(v) for i, v in amp.items()}


def reference_energy(f):
    """<psi|H|psi> + offset, or nan without /state_prep.

    Each Pauli term acts on the stored amplitudes via flip/phase
    bit masks -- an O(nterms * ndets) walk, no diagonalisation.
    """
    amp = _amplitudes(f)
    if amp is None:
        return math.nan
    g = f["pauli_hamil"]
    coeffs = g["coeffs"][...]
    paulis = g["paulis"][...]
    offset = float(g.attrs["offset"])
    e = 0.0 + 0j
    for ck, row in zip(coeffs, paulis):
        flip, phase = _masks(row)
        for idx_b, c_b in amp.items():
            c_a = amp.get(idx_b ^ flip)
            if c_a is None:
                continue
            r4 = (bin(flip & phase).count("1")
                  + 2 * bin(idx_b & phase).count("1")) & 3
            e += ck * np.conjugate(c_a) * c_b * _PHASE[r4]
    return e.real + offset


def _norm_offset(f):
    g = f["pauli_hamil"]
    return (float(g.attrs["normalization"]),
            float(g.attrs["offset"]))


def _values(f, group):
    """Complex values of a results group, NaN suffix trimmed."""
    if group not in f:
        raise Ph2Error(f"no /{group} group (run the simulation"
                       " first?)")
    raw = f[group]["values"][...]
    v = raw[:, 0] + 1j * raw[:, 1]
    nan = np.isnan(raw).any(axis=1)
    if nan.any():
        v = v[:int(np.argmax(nan))]
    return v


def _print_result(e, e_ref):
    print(f"{e},{e_ref},{e - e_ref}")


def cmd_fft(args):
    with open_ro(args.filename) as f:
        v = _values(f, args.group)
        delta = float(f[args.group].attrs["delta"])
        norm, offset = _norm_offset(f)
        e_ref = reference_energy(f)
    n = len(v)
    amp = np.abs(np.fft.fft(v))
    energy = np.fft.fftfreq(n) / (norm * delta / tau) + offset
    amp[0] = 0.0  # off-resonant baseline
    dom = int(np.argmax(amp))
    _print_result(energy[dom], e_ref)
    if args.peaks:
        try:
            from scipy.signal import find_peaks
        except ImportError:
            raise Ph2Error(
                "energy fft --peaks requires scipy"
                " (pip install -e \".[examples]\")") from None
        for p in sorted(find_peaks(amp)[0],
                        key=lambda b: energy[b]):
            print(f"{energy[p]},{amp[p]}")
    return 0


def _effective_time(group, attrs):
    if group == "circ_qdrift":
        return float(attrs["depth"]) \
            * math.asin(float(attrs["step_size"]))
    return float(attrs["steps"]) * float(attrs["angle_det"])


def cmd_mc(args):
    with open_ro(args.filename) as f:
        v = _values(f, args.group)
        T = _effective_time(args.group, f[args.group].attrs)
        norm, offset = _norm_offset(f)
        e_ref = reference_energy(f)
    z = v.mean()
    if abs(z) < 0.1:
        print(f"ph2: |mean overlap| = {abs(z):.3f} < 0.1;"
              " decohered average, estimate unreliable",
              file=sys.stderr)
    e = cmath.phase(z) / (T * norm) + offset
    _print_result(e, e_ref)
    return 0


def cmd_ref(args):
    with open_ro(args.filename) as f:
        print(reference_energy(f))
    return 0


def _rpe_ladder(samples):
    """Dyadic phase refinement -> theta (per-step phase)."""
    thetas = [0.0]
    for i, z in enumerate(samples):
        phi = cmath.phase(z)
        if phi < 0:
            phi += tau
        cand = [2 ** -i * (tau * k + phi) for k in range(2 ** i)]
        best = min(range(2 ** i),
                   key=lambda k: abs(thetas[i] - cand[k]))
        thetas.append(cand[best])
    th = thetas[len(samples) - 1]
    if th > math.pi:
        th = -(tau - th)
    return th


def cmd_rpe(args):
    with open_ro(args.filename) as f:
        if args.group not in f:
            raise Ph2Error(f"no /{args.group} group")
        raw = f[args.group]["values"][...]
        delta = float(f[args.group].attrs["delta"])
        norm, offset = _norm_offset(f)
        e_ref = reference_energy(f)
    j = int(math.log2(len(raw)))
    samples = []
    for i in range(j + 1):
        row = raw[2 ** i - 1]
        if np.isnan(row).any():
            raise Ph2Error(f"values[{2 ** i - 1}] is NaN"
                           " (run incomplete)")
        samples.append(complex(row[0], row[1]))
    th = _rpe_ladder(samples)
    x = 2.0 ** -j
    e0 = math.sqrt(1 - x * x) / x * math.tan(x * th) / norm / delta
    _print_result(e0 + offset, e_ref)
    return 0


def cmd_rpe_qdrift(args):
    j = math.ceil(math.log2(args.delta / args.epsilon))
    x = 2.0 ** -j
    with open_ro(f"{args.prefix}-0") as f:
        norm, offset = _norm_offset(f)
    samples = []
    for i in range(j + 1):
        with open_ro(f"{args.prefix}-{i}") as f:
            v = _values(f, "circ_qdrift")
            depth = int(f["circ_qdrift"].attrs["depth"])
        z = v.mean()
        print(f"ph2: level {i}: depth={depth} t={depth * x:g}"
              f" |z|={abs(z):.3f}", file=sys.stderr)
        samples.append(z)
    th = _rpe_ladder(samples)
    e0 = math.sqrt(1 - x * x) / x * math.tan(x * th) / norm
    print(f"{args.delta},{args.epsilon},{e0},{e0 + offset}")
    return 0
