#ifndef QDRIFT_H
#define QDRIFT_H

#include <stddef.h>

#include "phase2/circ.h"
#include "phase2/data.h"
#include "phase2/prob.h"
#include "xoshiro256ss.h"

/*
 * qDRIFT randomised product formula (Campbell, 2019).
 *
 * For each of `samples` independent runs, draw `depth`
 * Hamiltonian terms i.i.d. from the probability distribution
 * p_k = |c_k| / sum_j |c_j| and apply
 *
 *   exp(i * asin(step_size) * sign(c_k) * P_k)
 *
 * in draw order.  Output is the per-sample overlap
 * <psi|U|psi>, written to /circ_qdrift/values.
 *
 * The PRNG is xoshiro256** seeded by `seed` and split
 * deterministically across MPI ranks; `seed` must be
 * non-zero.
 */

struct qdrift_data {
	size_t depth;		/* terms drawn per sample */
	size_t samples;		/* number of independent samples */
	double step_size;	/* qDRIFT step size */
	uint64_t seed;		/* PRNG seed; must be non-zero */
};

struct qdrift_ranct {
	struct circ_hamil hm_ran;
	struct prob_cdf cdf;
};

struct qdrift {
	struct circ ct;
	struct qdrift_data dt;
	struct qdrift_ranct ranct;
	struct xoshiro256ss rng;
};

/* Load Hamiltonian and initial state, build the |c_k| CDF,
 * seed the PRNG.  Returns 0 on success, -1 on error. */
int qdrift_init(struct qdrift *qd, const struct qdrift_data *dt, data_id fid);

/* Release all resources held by `qd`. */
void qdrift_free(struct qdrift *qd);

/* Run `dt.samples` independent samples of depth `dt.depth`
 * and store per-sample overlaps in `ct.vals`.  Returns 0 on
 * success, -1 on error. */
int qdrift_simul(struct qdrift *qd);

/* Write step_size, depth, num_samples, seed attributes and
 * the overlap series to /circ_qdrift.  Returns 0 on success,
 * -1 on error. */
int qdrift_write_res(struct qdrift *qd, data_id fid);

#endif // QDRIFT_H
