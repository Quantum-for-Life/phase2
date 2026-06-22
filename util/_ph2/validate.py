"""Worksheet validation against doc/simul-h5-specs.md (ph2 validate).

check() returns a list of violation strings "H5PATH: message"; the
CLI prefixes each with the file name.  The checks mirror the spec
verbatim: dataset dtypes and shapes, attribute presence and kinds,
the exactly-one state-prep dispatch rule, the CSF component contract,
and the NaN-suffix (crash-consistency) guarantee on result values.
/circ_qdrift "seed" is accepted as optional: current ph2run writes
it, older files predate it.
"""

import sys

import numpy as np

from . import RESULTS_GROUPS, open_ro

_PAULI_BYTES = frozenset((0, 1, 2, 3))


def _is_scalar(attrs, name, kinds):
    return (name in attrs
            and np.ndim(attrs[name]) == 0
            and np.asarray(attrs[name]).dtype.kind in kinds)


def _check_hamil(f, out):
    if "pauli_hamil" not in f:
        out.append("/pauli_hamil: missing")
        return None
    g = f["pauli_hamil"]
    nqb = None
    if ("coeffs" not in g or g["coeffs"].ndim != 1
            or g["coeffs"].dtype != np.float64):
        out.append("/pauli_hamil/coeffs: not float64 (N,)")
    if ("paulis" not in g or g["paulis"].ndim != 2
            or g["paulis"].dtype != np.uint8):
        out.append("/pauli_hamil/paulis: not uint8 (N, nqb)")
    else:
        paulis = g["paulis"][...]
        nqb = paulis.shape[1]
        if "coeffs" in g and g["coeffs"].shape[0] != paulis.shape[0]:
            out.append("/pauli_hamil/coeffs:"
                       " length differs from paulis rows")
        if not set(np.unique(paulis)) <= _PAULI_BYTES:
            out.append("/pauli_hamil/paulis: byte outside {0,1,2,3}")
        if paulis.size and not paulis.any(axis=1).all():
            out.append("/pauli_hamil/paulis: all-identity row present")
    if not (_is_scalar(g.attrs, "normalization", "f")
            and g.attrs["normalization"] > 0):
        out.append("/pauli_hamil@normalization:"
                   " missing or not a positive double")
    if not _is_scalar(g.attrs, "offset", "f"):
        out.append("/pauli_hamil@offset: missing or not a double")
    return nqb


def _check_multidet(md, nqb, out):
    path = "/state_prep/multidet"
    m = None
    if ("coeffs" not in md or md["coeffs"].ndim != 2
            or md["coeffs"].shape[1] != 2
            or md["coeffs"].dtype != np.float64):
        out.append(f"{path}/coeffs: not float64 (M, 2)")
    else:
        m = md["coeffs"].shape[0]
        if m == 0:
            out.append(f"{path}: empty (M == 0)")
    if ("dets" not in md or md["dets"].ndim != 2
            or md["dets"].dtype != np.uint8):
        out.append(f"{path}/dets: not uint8 (M, nqb)")
    else:
        dets = md["dets"][...]
        if m is not None and dets.shape[0] != m:
            out.append(f"{path}/dets: row count differs from coeffs")
        if not set(np.unique(dets)) <= {0, 1}:
            out.append(f"{path}/dets: byte outside {{0,1}}")
        if nqb is not None and dets.shape[1] != nqb:
            out.append(f"{path}: nqb differs from /pauli_hamil")


def _check_cmat(g, path, n_sites, n_occ, out):
    if (g.ndim != 2 or g.dtype != np.float64
            or g.shape != (n_sites, n_occ)):
        out.append(f"{path}: not float64 ({n_sites}, {n_occ})")


def _check_coeff_block(g, path, ns, na, nb, closed, out):
    if "C_alpha" not in g:
        out.append(f"{path}/C_alpha: missing")
    else:
        _check_cmat(g["C_alpha"], f"{path}/C_alpha", ns, na, out)
    if closed:
        if "C_beta" in g:
            out.append(f"{path}/C_beta: present on closed shell")
    else:
        if "C_beta" not in g:
            out.append(f"{path}/C_beta: missing on open shell")
        else:
            _check_cmat(g["C_beta"], f"{path}/C_beta", ns, nb, out)


def _check_coeff_matrix(cm, out):
    path = "/state_prep/coeff_matrix"
    for name in ("n_qubits", "n_sites", "n_alpha", "n_beta"):
        if not _is_scalar(cm.attrs, name, "iu"):
            out.append(f"{path}@{name}: missing or not an integer")
            return
    for name in ("closed_shell", "tapered"):
        if (not _is_scalar(cm.attrs, name, "iu")
                or int(cm.attrs[name]) not in (0, 1)):
            out.append(f"{path}@{name}: missing or not in {{0,1}}")
            return
    ns, na, nb = (int(cm.attrs[k])
                  for k in ("n_sites", "n_alpha", "n_beta"))
    closed = bool(cm.attrs["closed_shell"])
    tapered = int(cm.attrs["tapered"])
    if int(cm.attrs["n_qubits"]) != 2 * ns - 2 * tapered:
        out.append(f"{path}@n_qubits: != 2*n_sites - 2*tapered")
    _check_coeff_block(cm, path, ns, na, nb, closed, out)
    if "csf" not in cm:
        return
    csf = cm["csf"]
    if not _is_scalar(csf.attrs, "n_components", "iu"):
        out.append(f"{path}/csf@n_components: missing")
        return
    ncomp = int(csf.attrs["n_components"])
    if ncomp < 1:
        out.append(f"{path}/csf@n_components: == 0"
                   " (empty CSF superposition)")
        return
    if set(csf.keys()) != {str(k) for k in range(ncomp)}:
        out.append(f"{path}/csf: subgroups not named"
                   f" 0..{ncomp - 1}")
        return
    for k in range(ncomp):
        bpath = f"{path}/csf/{k}"
        blk = csf[str(k)]
        if not _is_scalar(blk.attrs, "coefficient", "f"):
            out.append(f"{bpath}@coefficient:"
                       " missing or not a double")
        _check_coeff_block(blk, bpath, ns, na, nb, closed, out)


def _check_state_prep(f, nqb, hamil_only, out):
    if "state_prep" not in f:
        if not hamil_only:
            out.append("/state_prep: missing")
        return
    sp = f["state_prep"]
    has_md = "multidet" in sp
    has_cm = "coeff_matrix" in sp
    if has_md and has_cm:
        out.append("/state_prep:"
                   " both multidet and coeff_matrix present")
        return
    if not has_md and not has_cm:
        out.append("/state_prep:"
                   " no subtype (multidet or coeff_matrix)")
        return
    if has_md:
        _check_multidet(sp["multidet"], nqb, out)
    else:
        _check_coeff_matrix(sp["coeff_matrix"], out)


def _check_results(f, out):
    for name, attrs in RESULTS_GROUPS.items():
        if name not in f:
            continue
        try:
            g = f[name]
        except KeyError:
            continue  # dangling link, already reported
        path = f"/{name}"
        for a in attrs:
            if a not in g.attrs:
                out.append(f"{path}@{a}: missing")
        if ("values" not in g or g["values"].ndim != 2
                or g["values"].shape[1] != 2
                or g["values"].dtype != np.float64):
            out.append(f"{path}/values: not float64 (N, 2)")
            continue
        nan = np.isnan(g["values"][...]).any(axis=1)
        if nan.size and (nan[:-1] & ~nan[1:]).any():
            out.append(f"{path}/values:"
                       " NaN row precedes a real row")


def _dangling(f):
    """Top-level names whose links cannot be resolved (legacy
    soft-link artifacts; ph2run unlinks them on group create)."""
    bad = []
    for name in f:
        try:
            f[name]
        except KeyError:
            bad.append(name)
    return bad


def check(f, hamil_only=False):
    """Violations for an open worksheet, as 'H5PATH: message'."""
    out = []
    for name in _dangling(f):
        out.append(f"/{name}: dangling link")
        if name in ("pauli_hamil", "state_prep"):
            return out  # nothing further checkable
    nqb = _check_hamil(f, out)
    _check_state_prep(f, nqb, hamil_only, out)
    _check_results(f, out)
    return out


def cmd_validate(args):
    rc = 0
    for path in args.files:
        try:
            with open_ro(path) as f:
                violations = check(f, hamil_only=args.hamil_only)
        except Exception as e:
            print(f"ph2: {e}", file=sys.stderr)
            rc = 2
            continue
        for v in violations:
            print(f"{path}: {v}")
        if violations and rc == 0:
            rc = 1
    return rc
