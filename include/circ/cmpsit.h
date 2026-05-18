#ifndef CMPSIT_H
#define CMPSIT_H

#include <stddef.h>

#include "phase2.h"
#include "phase2/step_writer.h"
#include "xoshiro256ss.h"

/*
 * Composite (partially randomised) 2nd-order Trotter.
 *
 * The Hamiltonian is split into two disjoint parts:
 *
 *   - deterministic: the top `length` terms by |c_k|,
 *     sorted lexicographically and applied with step size
 *     `angle_det`;
 *   - randomised: the remaining terms, sampled qDRIFT-style
 *     with step size `angle_rand` to `depth` rotations.
 *
 * One step applies a forward half-sweep of one composite
 * sample at omega = 0.5, then a reverse half-sweep of an
 * independent sample — yielding the symmetric S_2 integrator.
 *
 * Output is the per-sample overlap series.  `seed` must be
 * non-zero.
 */

struct cmpsit_data {
	size_t samples;		/* number of independent samples */
	size_t steps;		/* number of Trotter steps */
	size_t length;		/* deterministic term count */
	size_t depth;		/* randomised term count per step */
	double angle_det;	/* deterministic step size */
	double angle_rand;	/* randomised step size */
	uint64_t seed;		/* PRNG seed; must be non-zero */
};

/* Working state for one sampled composite circuit. */
struct cmpsit_ranct {
	struct circ_hamil hm_det, hm_ran, hm_smpl;
	struct prob_cdf cdf;
	double lambda_r;
	size_t depth;
	double angle_rand;
};

struct cmpsit {
	struct circ ct;
	struct cmpsit_data dt;
	struct cmpsit_ranct ranct;
	struct xoshiro256ss rng;
	struct phase2_step_writer *sw;
};

/* Adopt the Hamiltonian + state-prep into a fresh
 * cmpsit context, split the Hamiltonian into
 * deterministic (top `dt->length` by |c_k|) and
 * randomised parts, build the importance-sampling CDF
 * on the randomised pool, and seed the per-context
 * PRNG from `dt->seed`.  Ownership of `hm` and
 * `*sp_data` transfers (see circ_init).  `sw` is held
 * by reference; its lifetime must outlive
 * cmpsit_simul.  Returns 0 on success, -1 on error. */
int cmpsit_init(struct cmpsit *cp, const struct cmpsit_data *dt,
	struct circ_hamil hm, enum stprep_kind sp_kind, const void *sp_data,
	struct phase2_step_writer *sw);

/* Release all resources held by `cp`. */
void cmpsit_free(struct cmpsit *cp);

/* Run `dt.samples` independent composite samples; each
 * sample applies `dt.steps` symmetric Trotter steps
 * over the deterministic + randomised pools.  Per-
 * sample overlap is stored in `ct.vals` and forwarded
 * through `cp->sw` when non-NULL.  Returns 0 on
 * success, -1 on error. */
int cmpsit_simul(struct cmpsit *cp);

#endif // CMPSIT_H
