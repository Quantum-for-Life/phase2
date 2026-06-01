"""Reference-state energy <psi|H|psi> from a phase2 HDF5 file.

The analysis scripts report their estimate against the energy of the
trial state itself, so a reader sees both the number and how far the
evolution moved it.  H is the physical (identity-removed) Hamiltonian
stored in /pauli_hamil as raw `coeffs` and per-qubit `paulis`
(0=I,1=X,2=Y,3=Z); |psi> is the /state_prep/multidet expansion.  For a
determinant |b>, P_k|b> = z |b ^ flip_k> with z a fourth root of unity,
so

    <psi|H|psi> = sum_k coeff_k sum_{a,b} conj(c_a) c_b <a|P_k|b> + offset

is an O(nterms * ndets) walk -- no diagonalisation.  The result is real
for a Hermitian H with real coefficients; the imaginary part is dropped.
"""

import h5py
import numpy as np


def _masks(row):
    """(flip, phase) bit masks for a per-qubit Pauli byte row.

    flip bit set for X(1)/Y(2); phase bit set for Z(3)/Y(2).
    """
    flip = phase = 0
    for j, b in enumerate(row):
        if b in (1, 2):
            flip |= 1 << j
        if b in (2, 3):
            phase |= 1 << j
    return flip, phase


_PHASE = (1.0 + 0j, 1j, -1.0 + 0j, -1j)


def reference_energy(filename):
    """Return <psi|H|psi> + offset (Hartree) for the file's trial state."""
    with h5py.File(filename, "r") as f:
        coeffs = f["pauli_hamil/coeffs"][...]
        paulis = f["pauli_hamil/paulis"][...]
        offset = float(f["pauli_hamil"].attrs["offset"])
        det_cf = f["state_prep/multidet/coeffs"][...]
        dets = f["state_prep/multidet/dets"][...]

    amp = {}
    for occ, (re, im) in zip(dets, det_cf):
        idx = int(sum(int(b) << j for j, b in enumerate(occ)))
        amp[idx] = complex(re, im)

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
