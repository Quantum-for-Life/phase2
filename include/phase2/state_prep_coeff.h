#ifndef PHASE2_STATE_PREP_COEFF_H
#define PHASE2_STATE_PREP_COEFF_H

#include <stddef.h>
#include <stdint.h>

#include "phase2/qreg.h"

struct data_coeff_matrix;

/*
 * Slater-Condon expansion of a coefficient-matrix trial
 * state into a dense MPI-distributed register, plus the
 * matching inner-product readout.
 *
 * Caller-supplied scratch carries the alpha / beta
 * k-subset tuples (pure combinatorics; filled once at
 * init) and the det workspaces (refilled per call from
 * the current block's C_alpha / C_beta).  One scratch
 * services any number of _expand / _inner calls sharing
 * the same (n_sites, n_alpha, n_beta).
 *
 * `_expand` writes (accumulate=0) or adds (accumulate=1)
 * weighted Slater-Condon products
 *     c = w * det(C_a[occ_a]) * det(C_b[occ_b])
 * into the owning rank's slice of reg->amp[].  Sparsity
 * prune: |c| < 1e-12 is skipped.
 *
 * `_expand_all` is the dispatcher used by circ_prepst():
 * always zeros the register first, then either runs one
 * single-block expand (accumulate=0) or accumulates over
 * cm->blocks[].
 *
 * `_inner` writes <trial | reg> through `*out`.  The
 * trial state is real (C is real), so conjugation is
 * purely on the register amplitudes.  Weight supports
 * CSF superposition: the full inner product is the
 * weighted sum over blocks.
 */

struct state_prep_coeff_scratch {
	uint32_t n_sites, n_alpha, n_beta;
	uint32_t *tup_a, *tup_b;
	double   *det_a, *det_b;
	size_t    Ma, Mb;
};

/* Allocate scratch buffers and fill the (n_sites,
 * n_alpha) / (n_sites, n_beta) k-subset tuples.  These
 * are pure combinatorics and are reused across every
 * _expand / _inner call against this scratch.  Returns
 * 0 on success, -1 on size-bound violation or alloc
 * failure. */
int state_prep_coeff_scratch_init(struct state_prep_coeff_scratch *sc,
	uint32_t n_sites, uint32_t n_alpha, uint32_t n_beta);

void state_prep_coeff_scratch_free(struct state_prep_coeff_scratch *sc);

/* Returns 0 on success, -1 on error. */
int state_prep_coeff_expand(struct qreg *reg,
	struct state_prep_coeff_scratch *sc,
	const double *C_alpha, const double *C_beta,
	double weight, int tapered, int accumulate);

int state_prep_coeff_expand_all(struct qreg *reg,
	struct state_prep_coeff_scratch *sc,
	const struct data_coeff_matrix *cm);

int state_prep_coeff_inner(struct qreg *reg,
	struct state_prep_coeff_scratch *sc,
	const double *C_alpha, const double *C_beta,
	double weight, int tapered, _Complex double *out);

#endif /* PHASE2_STATE_PREP_COEFF_H */
