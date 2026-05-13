# coeff_matrix_reference.py
#
# Independent Python implementation of the Slater-Condon
# expansion used by `state_prep_coeff_expand` in the C path.
# Used by `t-ref-coeff_matrix.py` as a verification oracle:
# the C and Python expansions must agree bit-for-bit (within
# 1e-12) on every shipped fixture.
#
# `expand_coeff_matrix(n_sites, n_alpha, n_beta, C_alpha,
# C_beta_or_None, tapered)` returns a {int_index:
# float_amplitude} dict using the same bitstring convention as
# the C `occ_pair_to_idx`:
#   * alpha qubits at bits [0, n_sites)
#   * beta  qubits at bits [n_sites, 2*n_sites)
#   * tapered drops bit 0 and bit n_sites, shifts upper bits down
# Amplitude convention: c(occ_a, occ_b) = det(C_a[occ_a, :]) *
#   det(C_b[occ_b, :]).
# Sparsity prune at |c| < 1e-12 matches the C path.

from itertools import combinations

import numpy as np


def expand_coeff_matrix(
    n_sites: int,
    n_alpha: int,
    n_beta: int,
    C_alpha: np.ndarray,
    C_beta: np.ndarray | None,
    tapered: bool,
    weight: float = 1.0,
) -> dict[int, float]:
    Cb = C_beta if C_beta is not None else C_alpha
    out: dict[int, float] = {}

    alpha = []
    for occ in combinations(range(n_sites), n_alpha):
        if n_alpha == 0:
            d = 1.0
        else:
            sub = C_alpha[list(occ), :]
            d = float(np.linalg.det(sub))
        alpha.append((occ, d))

    beta = []
    for occ in combinations(range(n_sites), n_beta):
        if n_beta == 0:
            d = 1.0
        else:
            sub = Cb[list(occ), :]
            d = float(np.linalg.det(sub))
        beta.append((occ, d))

    for occ_a, da in alpha:
        for occ_b, db in beta:
            cf = weight * da * db
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
