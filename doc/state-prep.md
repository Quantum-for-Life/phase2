# State preparation: algorithm reference

This document is the algorithmic companion to
[simul-h5-specs.md](simul-h5-specs.md).  Where the spec is a
schema (which fields go where, what dtype, what shape), this
document describes how the simulator turns those fields into
a dense MPI-distributed state vector and why each step is what
it is.

## Purpose and audience

The audience is implementers porting phase2 to a new backend
(GPU, gate-based, sparse), reviewers vetting the v0.5 PPP
integration, and campaign designers needing to verify
reproducibility against the upstream Python reference.

Two state-prep subtypes are supported:

- `/state_prep/multidet/`: flat list of (bitstring, complex
  amplitude) pairs.  Used by every v0.4 Bendazzoli MO
  campaign; remains the default for small `M`.
- `/state_prep/coeff_matrix/`: a real `n_sites x n_occ`
  coefficient matrix from which the simulator expands all
  `M = C(n_sites, n_alpha) * C(n_sites, n_beta)` amplitudes
  at runtime via Slater-Condon.  Used by v0.5 PPP cells
  where `M` would otherwise be unmanageable on disk and in
  the pak.

## Bitstring-to-index convention

Each computational basis state of the `n_qubits`-qubit
register maps to an integer `idx` in `[0, 2^n_qubits)` by
setting bit `q` of `idx` to the occupation of qubit `q`
(`LSB-first`, identical to the v0.4 multidet path):

    idx = sum_q (occupation_q << q)

For the coeff_matrix path:

- qubits `[0, n_sites)` are the alpha-spin block,
- qubits `[n_sites, 2 * n_sites)` are the beta-spin block,
- `occ_pair_to_idx(occ_a, occ_b, n_sites)` sets bit `q` for
  each `q` in `occ_a` and bit `n_sites + q` for each `q` in
  `occ_b`.

Example: `n_sites = 4`, `occ_a = (0, 1)`, `occ_b = (0, 1)`,
`idx = 0b00110011 = 51`.

## MPI ownership decomposition

`struct qreg` splits the qubit index into `qb_lo` low bits
(amplitude index on this rank) and `qb_hi` high bits (rank
selector).  The rank that owns `idx` is

    rank = qreg_getihi(reg, idx) = (idx >> qb_lo) & ((1 << qb_hi) - 1)

and the local offset within that rank's amp buffer is

    i_lo = qreg_getilo(reg, idx) = idx & ((1 << qb_lo) - 1)

Each rank runs the full Slater-Condon outer product but only
writes the slice it owns: an inline `qreg_getihi(reg, idx) ==
my_rank` filter inside `state_prep_coeff_expand` replaces the
per-amplitude `MPI_Barrier` of the v0.4 `qreg_setamp` path
with one barrier at the end of the expansion.

## Slater-Condon expansion algorithm

The Slater-Condon formula for a half-filled real determinant
is

    c(occ_a, occ_b) = det(C_a[occ_a, :]) * det(C_b[occ_b, :])

where `C_a[occ_a, :]` is the `n_alpha x n_alpha` submatrix of
`C_alpha` selected by the alpha occupation pattern, and
similarly for beta.

The expansion (`state_prep_coeff_expand`) proceeds in three
phases:

1. **Precompute** `det_alpha[i]` and `det_beta[j]` arrays
   walking the canonical lex-ascending enumeration of
   k-subsets (`combinations_init`/`combinations_next`,
   Knuth TAOCP 4A §7.2.1.3 algorithm L in lex form).  The
   determinant routine `det_small` is a self-contained LU
   with partial pivoting; no LAPACK link, bounded for
   `n <= DET_SMALL_MAX_N`.
2. **Outer product**: iterate `(i, j)` over `M_alpha x
   M_beta`, form `cf = weight * det_alpha[i] *
   det_beta[j]`, apply the sparsity prune, compute the
   bitstring index, optionally drop the tapered bits, and
   filter on the rank owner.
3. **Barrier**: one `MPI_Barrier(MPI_COMM_WORLD)` after the
   outer product completes.

Spin-orbital ordering: the alpha block occupies qubits
`[0, n_sites)`, the beta block occupies qubits `[n_sites,
2 * n_sites)`.  The Slater-Condon product is taken in this
fixed order; no cross-block Jordan-Wigner sign appears
because alpha and beta blocks are encoded in disjoint qubit
ranges and never reorder.

## Sparsity-prune threshold

Amplitudes with `|cf| < 1e-12` are skipped.  This matches the
upstream Python reference at
`phase2/test/ref/psi_json_reference.py` and avoids polluting
the dense state with denormals from numerically-zero
determinants.  The threshold is exercised in
`t-state_prep_coeff_expand.c` (test `t_fixture` cross-checks
non-pruned amplitudes; the synthetic test inputs include
matrices that legitimately produce sub-threshold products).

If a future cell legitimately produces amplitudes smaller
than the threshold, the constant `SPARSITY_PRUNE` in
`phase2/state_prep_coeff.c` is the single tuning knob.

## Tapering

Tapered states drop two redundant qubits per the Z₂ + S_z
parity reduction used in the PPP encoding:

- bit `0` (the alpha lowest-site qubit), and
- bit `n_sites` (the beta lowest-site qubit).

`drop_two_bits(idx, n_sites)` performs the drop:

    low  = (idx >> 1) & ((1 << (n_sites - 1)) - 1)
    high = idx >> (n_sites + 1)
    out  = low | (high << (n_sites - 1))

The drop happens at expansion time, per generated bitstring;
the H5 file still carries the full `n_sites` MO matrix.  The
tapered convention is identical to the v0.4 `.inputst`
tapered file encoding, just applied at read time instead of
write time.

## CSF multi-block superposition

For some open-shell singlets, a single Slater determinant
cannot enforce S² = 0 simultaneously with the site-symmetry
constraints; a finite linear combination of CSF blocks is
required.  The schema reflects this via the optional
`csf/` subgroup; the dispatcher
`state_prep_coeff_expand_all` calls `state_prep_coeff_expand`
once per component with `accumulate = 1`, accumulating into
the same register.

Each component carries its own `coefficient` (real),
`C_alpha`, and (for open-shell) `C_beta`.  The superposition
is assumed to be normalised; individual components are not.

Component ordering is the natural integer order `0/`, `1/`,
..., `(n_components - 1)/`.

## Complexity

- **Precompute**: `O(M_alpha * n_alpha^3 + M_beta *
  n_beta^3)`, dominated by the per-subset LU.
- **Scatter**: `O(M_alpha * M_beta)` arithmetic per rank
  (sparsity prune trims the constant), plus a single MPI
  barrier.
- **Memory**: `O(M_alpha + M_beta)` for the determinant
  tables and `O(M_alpha * k_a + M_beta * k_b)` for the
  enumerated tuples.

`state_prep_coeff_inner` reuses the same precompute and the
same outer-product walk, so its cost is identical to the
expansion modulo the per-amplitude `qreg_setamp` write
being replaced by a `reg->amp[i_lo]` read.  An O(N)
measurement path that caches the expanded trial state in a
second buffer is feasible but not implemented; until then,
every inner-product evaluation re-runs the full outer
product.

At `n_sites = 18` tapered (`M_total ~= 2.36e9` after the
half-filling tableau), the per-rank walk costs roughly 70 s
once with the current arithmetic; this is a one-time cost
before the long Trotter loop.

## Numerical conditioning

`det_small` uses Gaussian elimination with partial pivoting,
sufficient when `C` is well-conditioned (PPP MO orthogonality
gives column-orthonormal blocks).  Condition number bounds at
`n_alpha <= 17` keep relative error well below `1e-12`;
`t-det_small.c` cross-checks the routine against an
independent naive implementation on 32 random matrices per
size.

If a future cell exhibits a near-singular `C`, the partial-
pivoting routine will lose digits but not crash; the
expansion will still complete with reduced precision.  The
sparsity-prune threshold acts as a coarse safety net for
underflow.

## Cross-validation

The vendored Python reference at
`phase2/test/ref/psi_json_reference.py` is a verbatim port of
the upstream tmm-chomp issue #14 implementation.
`t-ref-psi_json.py` is a strict term-by-term oracle: for
every shipped fixture (`N4_closed`, `N4_open`, `N4_csf`,
`N8_untapered`, `N8_tapered`) it

1. invokes `t-state_prep_coeff_expand --dump <out> <fixture>`
   to capture the C-side Slater-Condon expansion as
   `(idx, re, im)` lines,
2. runs `expand_coeff_matrix` from `psi_json_reference.py`
   on the same C matrices read straight from the fixture H5,
3. compares the two index sets exactly and asserts the
   per-index abs-diff is `<= 1e-12`.

On mismatch the script names the fixture, the worst-case
index, both amplitudes, and the diff, then exits non-zero.
The exit code is the gate; a failed comparison fails the
build.  Observed worst-case diff on the current fixtures is
`~ 6e-17`, well under the bound.

For end-to-end validation, the `N4_closed.h5` (coeff_matrix)
and `N4_multidet.h5` (multidet) fixtures encode the same
physical state.  `t-circ_trott_coeff.c` and
`t-circ_trott2_coeff.c` run one Trotter step on each and
verify `/circ_trott/values` agrees to `1e-13`.

The MPI scaling harness `phase2/test/run-mpi-scaling.sh`
verifies that the coeff_matrix path produces identical
results across `R in {1, 2, 4, 8}` ranks.

## Public API summary

Public symbols introduced by this change set:

- `enum stprep_kind` with `STPREP_MULTIDET`,
  `STPREP_COEFF_MATRIX`
- `data_state_prep_kind`
- `data_coeff_matrix_getnums`
- `data_coeff_matrix_read`
- `data_coeff_matrix_csf_count`
- `data_coeff_matrix_csf_read`
- `det_small`
- `combinations_init`, `combinations_next`
- `state_prep_coeff_expand`
- `state_prep_coeff_expand_all`
- `state_prep_coeff_inner`
- `circ_coeff_init`, `circ_coeff_free`

The coverage gate `phase2/test/check-docs.sh` greps for each
of these symbols in the documentation set and fails if any
are missing.
