#ifndef PHASE2_STATE_PREP_COEFF_H
#define PHASE2_STATE_PREP_COEFF_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/data.h"
#include "phase2/qreg.h"

/*
 * state_prep_coeff - Slater-Condon expansion of a coefficient-
 * matrix trial state into a dense MPI-distributed register.
 *
 * Public API:
 *   int state_prep_coeff_expand(struct qreg *reg,
 *           uint32_t n_sites, uint32_t n_alpha, uint32_t n_beta,
 *           const double *C_alpha, const double *C_beta,
 *           double weight, int tapered, int accumulate);
 *   int state_prep_coeff_expand_all(struct qreg *reg,
 *           const struct circ_coeff *cm);
 *   int circ_coeff_init(struct circ_coeff *cm, data_id fid);
 *   void circ_coeff_free(struct circ_coeff *cm);
 *
 * `state_prep_coeff_expand` walks the cartesian product of
 * alpha and beta k-subsets, evaluates the Slater-Condon product
 *
 *     c(occ_a, occ_b) = w * det(C_a[occ_a, :]) * det(C_b[occ_b, :])
 *
 * and writes (or accumulates) the amplitudes into the slice of
 * `reg->amp[]` owned by this MPI rank.  A single MPI_Barrier
 * is issued at the end.
 *
 * `accumulate == 0` clears each owned slot with the computed
 * coefficient.  `accumulate == 1` adds to the existing slot, as
 * required by the CSF superposition path.  Both modes apply the
 * sparsity prune |c| < 1e-12.
 *
 * `state_prep_coeff_expand_all` is the dispatcher used by
 * `circ_prepst()`: single-block files go through one call with
 * accumulate=0; CSF files zero the register and accumulate.
 *
 * Algorithm reference: phase2/doc/state-prep.md.
 */

struct circ_coeff_block {
	double cf;
	double *C_alpha;
	double *C_beta;
};

struct circ_coeff {
	uint32_t n_qubits;
	uint32_t n_sites;
	uint32_t n_alpha;
	uint32_t n_beta;
	int closed_shell;
	int tapered;
	size_t n_components;
	double *C_alpha;
	double *C_beta;
	struct circ_coeff_block *blocks;
};

int state_prep_coeff_expand(struct qreg *reg, uint32_t n_sites,
	uint32_t n_alpha, uint32_t n_beta, const double *C_alpha,
	const double *C_beta, double weight, int tapered, int accumulate);

int state_prep_coeff_expand_all(struct qreg *reg, const struct circ_coeff *cm);

/*
 * Inner product <trial | reg> over the coefficient-matrix trial
 * state described by (n_sites, n_alpha, n_beta, C_alpha, C_beta,
 * weight, tapered).
 *
 * Walks the same outer product used at expand time, sums
 * weight * conj(C-amplitude(idx)) * reg->amp[idx] over the
 * amplitudes owned by this rank, and MPI_Allreduces the
 * partial sum.  Returns the inner product on all ranks.
 *
 * The trial state is taken to be real (C is real), so the
 * conjugation is purely on the register amplitudes.  The
 * weight argument supports the CSF superposition: a CSF inner
 * product is the weighted sum over blocks.
 */
_Complex double state_prep_coeff_inner(struct qreg *reg, uint32_t n_sites,
	uint32_t n_alpha, uint32_t n_beta, const double *C_alpha,
	const double *C_beta, double weight, int tapered);

int circ_coeff_init(struct circ_coeff *cm, data_id fid);

void circ_coeff_free(struct circ_coeff *cm);

#endif /* PHASE2_STATE_PREP_COEFF_H */
