#!/usr/bin/env python3
# build_coeff_fixtures.py
#
# Deterministic builder for the coeff_matrix HDF5 fixtures used
# by the phase2 unit and integration tests.
#
# Run once after the phase2 code changes land, then commit the
# resulting .h5 files.  Each fixture is generated from a fixed
# seed so the script is byte-stable across runs.
#
# Fixtures emitted (alongside this script):
#   N4_closed.h5            closed-shell, untapered, M=36
#   N4_open.h5              open-shell, n_alpha != n_beta
#   N4_csf.h5               two-component CSF, closed-shell
#   N8_untapered.h5         closed-shell, untapered, M=4900
#   N8_tapered.h5           closed-shell, tapered, n_qubits=14
#   N4_multidet.h5          same physical state as N4_closed.h5
#                           but encoded as /state_prep/multidet/
#   N4_both.h5              both /state_prep/multidet/ and
#                           /state_prep/coeff_matrix/ present
#                           (intentionally invalid; used for the
#                           ambiguity test)
#   N4_csf_empty.h5         coeff_matrix with an empty csf/
#                           subgroup (n_components=0).  Used to
#                           exercise the "explicit csf/ with no
#                           components" rejection path.
#
# Each closed/untapered file also carries a minimal /pauli_hamil/
# group so that the integration tests (t-circ_trott_coeff,
# t-circ_trott2_coeff, t-circ_prepst_coeff) can run circ_init()
# end to end without needing a separately shipped hamilop.

import math
import sys
from itertools import combinations
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent


def huckel_C(n_sites: int, n_occ: int, seed: int) -> np.ndarray:
    """Return an orthonormal (n_sites, n_occ) block of MO
    coefficients with deterministic structure."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n_sites, n_sites))
    A = A + A.T
    _, V = np.linalg.eigh(A)
    return V[:, :n_occ].copy()


def slater_condon_expand(
    n_sites: int,
    n_alpha: int,
    n_beta: int,
    C_alpha: np.ndarray,
    C_beta: np.ndarray | None,
    tapered: bool,
):
    """Reference Slater-Condon expansion.  Returns dict
    {index: amplitude}."""
    Cb = C_beta if C_beta is not None else C_alpha
    out: dict[int, float] = {}

    alpha_dets = []
    for occ in combinations(range(n_sites), n_alpha):
        sub = C_alpha[list(occ), :] if n_alpha else np.zeros((0, 0))
        d = 1.0 if n_alpha == 0 else float(np.linalg.det(sub))
        alpha_dets.append((occ, d))

    beta_dets = []
    for occ in combinations(range(n_sites), n_beta):
        sub = Cb[list(occ), :] if n_beta else np.zeros((0, 0))
        d = 1.0 if n_beta == 0 else float(np.linalg.det(sub))
        beta_dets.append((occ, d))

    for occ_a, da in alpha_dets:
        for occ_b, db in beta_dets:
            cf = da * db
            if abs(cf) < 1e-12:
                continue
            idx = 0
            for q in occ_a:
                idx |= 1 << q
            for q in occ_b:
                idx |= 1 << (n_sites + q)
            if tapered:
                low = (idx >> 1) & ((1 << (n_sites - 1)) - 1)
                high = idx >> (n_sites + 1)
                idx = low | (high << (n_sites - 1))
            out[idx] = out.get(idx, 0.0) + cf
    return out


def write_pauli_hamil(g, n_qubits, terms, offset=0.0, norm=1.0):
    """terms: list of (coeff, [pauli per qubit])."""
    n = len(terms)
    coeffs = np.array([t[0] for t in terms], dtype=np.float64)
    paulis = np.array([t[1] for t in terms], dtype=np.uint8)
    g.attrs["offset"] = np.float64(offset)
    g.attrs["normalization"] = np.float64(norm)
    g.create_dataset("coeffs", data=coeffs)
    g.create_dataset("paulis", data=paulis)


def attach_default_hamil(f, n_qubits, seed):
    """Attach a tiny deterministic Hamiltonian (3 Z terms)
    sufficient for the trott circuit but not physical."""
    rng = np.random.default_rng(seed)
    g = f.create_group("pauli_hamil")
    terms = []
    for k in range(3):
        coeff = float(rng.uniform(0.1, 0.5))
        paus = [0] * n_qubits
        paus[k % n_qubits] = 3
        terms.append((coeff, paus))
    write_pauli_hamil(g, n_qubits, terms)


def write_coeff_matrix(
    f,
    n_qubits,
    n_sites,
    n_alpha,
    n_beta,
    C_alpha,
    C_beta,
    closed_shell,
    tapered,
    csf=None,
):
    sp = f.create_group("state_prep")
    cm = sp.create_group("coeff_matrix")
    cm.attrs["n_qubits"] = np.uint32(n_qubits)
    cm.attrs["n_sites"] = np.uint32(n_sites)
    cm.attrs["n_alpha"] = np.uint32(n_alpha)
    cm.attrs["n_beta"] = np.uint32(n_beta)
    cm.attrs["closed_shell"] = np.uint8(1 if closed_shell else 0)
    cm.attrs["tapered"] = np.uint8(1 if tapered else 0)
    cm.create_dataset("C_alpha", data=np.asarray(C_alpha, dtype=np.float64))
    if not closed_shell:
        cm.create_dataset(
            "C_beta", data=np.asarray(C_beta, dtype=np.float64)
        )
    if csf is not None:
        csf_g = cm.create_group("csf")
        csf_g.attrs["n_components"] = np.uint32(len(csf))
        for k, (coef, Ca_k, Cb_k) in enumerate(csf):
            g = csf_g.create_group(str(k))
            g.attrs["coefficient"] = np.float64(coef)
            g.create_dataset(
                "C_alpha", data=np.asarray(Ca_k, dtype=np.float64)
            )
            if not closed_shell:
                g.create_dataset(
                    "C_beta", data=np.asarray(Cb_k, dtype=np.float64)
                )


def write_multidet(f, n_qubits, amps: dict[int, float]):
    sp = f.require_group("state_prep")
    md = sp.create_group("multidet")
    items = sorted(amps.items())
    n = len(items)
    coeffs = np.zeros((n, 2), dtype=np.float64)
    dets = np.zeros((n, n_qubits), dtype=np.uint8)
    for r, (idx, cf) in enumerate(items):
        coeffs[r, 0] = cf
        coeffs[r, 1] = 0.0
        for j in range(n_qubits):
            dets[r, j] = (idx >> j) & 1
    md.create_dataset("coeffs", data=coeffs)
    md.create_dataset("dets", data=dets)


def build_N4_closed():
    n_sites, n_alpha, n_beta = 4, 2, 2
    C = huckel_C(n_sites, n_alpha, seed=4001)
    amps = slater_condon_expand(n_sites, n_alpha, n_beta, C, None, False)
    path = HERE / "N4_closed.h5"
    with h5py.File(path, "w") as f:
        write_coeff_matrix(
            f,
            n_qubits=2 * n_sites,
            n_sites=n_sites,
            n_alpha=n_alpha,
            n_beta=n_beta,
            C_alpha=C,
            C_beta=None,
            closed_shell=True,
            tapered=False,
        )
        attach_default_hamil(f, 2 * n_sites, seed=4002)
    return path, amps


def build_N4_multidet(amps_from_closed):
    n_sites = 4
    path = HERE / "N4_multidet.h5"
    with h5py.File(path, "w") as f:
        write_multidet(f, 2 * n_sites, amps_from_closed)
        attach_default_hamil(f, 2 * n_sites, seed=4002)
    return path


def build_N4_open():
    n_sites, n_alpha, n_beta = 4, 2, 1
    Ca = huckel_C(n_sites, n_alpha, seed=4003)
    Cb = huckel_C(n_sites, n_beta, seed=4004)
    path = HERE / "N4_open.h5"
    with h5py.File(path, "w") as f:
        write_coeff_matrix(
            f,
            n_qubits=2 * n_sites,
            n_sites=n_sites,
            n_alpha=n_alpha,
            n_beta=n_beta,
            C_alpha=Ca,
            C_beta=Cb,
            closed_shell=False,
            tapered=False,
        )
    return path


def build_N4_csf():
    n_sites, n_alpha, n_beta = 4, 2, 2
    C0 = huckel_C(n_sites, n_alpha, seed=4101)
    C1 = huckel_C(n_sites, n_alpha, seed=4102)
    path = HERE / "N4_csf.h5"
    with h5py.File(path, "w") as f:
        write_coeff_matrix(
            f,
            n_qubits=2 * n_sites,
            n_sites=n_sites,
            n_alpha=n_alpha,
            n_beta=n_beta,
            C_alpha=C0,
            C_beta=None,
            closed_shell=True,
            tapered=False,
            csf=[(0.6, C0, None), (0.8, C1, None)],
        )
    return path


def build_N8_untapered():
    n_sites, n_alpha, n_beta = 8, 4, 4
    C = huckel_C(n_sites, n_alpha, seed=8001)
    path = HERE / "N8_untapered.h5"
    with h5py.File(path, "w") as f:
        write_coeff_matrix(
            f,
            n_qubits=2 * n_sites,
            n_sites=n_sites,
            n_alpha=n_alpha,
            n_beta=n_beta,
            C_alpha=C,
            C_beta=None,
            closed_shell=True,
            tapered=False,
        )
        attach_default_hamil(f, 2 * n_sites, seed=8002)
    return path


def build_N8_tapered():
    n_sites, n_alpha, n_beta = 8, 4, 4
    C = huckel_C(n_sites, n_alpha, seed=8101)
    path = HERE / "N8_tapered.h5"
    with h5py.File(path, "w") as f:
        write_coeff_matrix(
            f,
            n_qubits=2 * n_sites - 2,
            n_sites=n_sites,
            n_alpha=n_alpha,
            n_beta=n_beta,
            C_alpha=C,
            C_beta=None,
            closed_shell=True,
            tapered=True,
        )
        attach_default_hamil(f, 2 * n_sites - 2, seed=8102)
    return path


def build_N4_csf_empty():
    """Intentionally invalid fixture: csf/ subgroup present but
    n_components=0.  The reader must reject this."""
    n_sites, n_alpha, n_beta = 4, 2, 2
    C = huckel_C(n_sites, n_alpha, seed=4301)
    path = HERE / "N4_csf_empty.h5"
    with h5py.File(path, "w") as f:
        sp = f.create_group("state_prep")
        cm = sp.create_group("coeff_matrix")
        cm.attrs["n_qubits"] = np.uint32(2 * n_sites)
        cm.attrs["n_sites"] = np.uint32(n_sites)
        cm.attrs["n_alpha"] = np.uint32(n_alpha)
        cm.attrs["n_beta"] = np.uint32(n_beta)
        cm.attrs["closed_shell"] = np.uint8(1)
        cm.attrs["tapered"] = np.uint8(0)
        cm.create_dataset("C_alpha", data=np.asarray(C, dtype=np.float64))
        csf_g = cm.create_group("csf")
        csf_g.attrs["n_components"] = np.uint32(0)
    return path


def build_N4_both():
    """Intentionally invalid fixture: both subgroups present."""
    n_sites, n_alpha, n_beta = 4, 2, 2
    C = huckel_C(n_sites, n_alpha, seed=4201)
    amps = slater_condon_expand(n_sites, n_alpha, n_beta, C, None, False)
    path = HERE / "N4_both.h5"
    with h5py.File(path, "w") as f:
        write_coeff_matrix(
            f,
            n_qubits=2 * n_sites,
            n_sites=n_sites,
            n_alpha=n_alpha,
            n_beta=n_beta,
            C_alpha=C,
            C_beta=None,
            closed_shell=True,
            tapered=False,
        )
        write_multidet(f, 2 * n_sites, amps)
    return path


def main():
    p1, amps = build_N4_closed()
    p2 = build_N4_multidet(amps)
    p3 = build_N4_open()
    p4 = build_N4_csf()
    p5 = build_N8_untapered()
    p6 = build_N8_tapered()
    p7 = build_N4_both()
    p8 = build_N4_csf_empty()
    for p in (p1, p2, p3, p4, p5, p6, p7, p8):
        print(p)


if __name__ == "__main__":
    main()
