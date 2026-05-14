#ifndef CMPSIT_H
#define CMPSIT_H

#include <stddef.h>

#include "phase2.h"
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
 * Output is the per-sample overlap series written to
 * /circ_cmpsit/values.  `seed` must be non-zero.
 */

struct cmpsit_data {
	uint64_t seed;		/* PRNG seed; must be non-zero */
	size_t length;		/* deterministic term count */
	size_t depth;		/* randomised term count per step */
	size_t steps;		/* number of Trotter steps */
	double angle_det;	/* deterministic step size */
	double angle_rand;	/* randomised step size */
	size_t samples;		/* number of independent samples */
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
	data_id fid;	/* output file; 0 means "no per-step writes" */
};

/* Load the Hamiltonian and initial state, split into
 * deterministic / randomised parts, sort the deterministic
 * part lexicographically, build the randomised CDF, seed the
 * PRNG, create the /circ_cmpsit output group with NaN-padded
 * values dataset and the scalar attributes.  Returns 0 on
 * success, -1 on error. */
int cmpsit_init(struct cmpsit *cp, const struct cmpsit_data *dt, data_id fid);

/* Release all resources held by `cp`. */
void cmpsit_free(struct cmpsit *cp);

/* Run `dt.samples` independent samples; per-sample overlap
 * is stored in `ct.vals` and written to
 * /circ_cmpsit/values[i] one row at a time (rank-0-only,
 * fflush per sample).  Returns 0 on success, -1 on error. */
int cmpsit_simul(struct cmpsit *cp);

#endif // CMPSIT_H
