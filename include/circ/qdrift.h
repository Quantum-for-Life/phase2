#ifndef QDRIFT_H
#define QDRIFT_H

#include <stddef.h>

#include "phase2/circ.h"
#include "phase2/step_writer.h"
#include "prob.h"
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
 * <psi|U|psi>.
 *
 * The PRNG is xoshiro256** seeded by `seed` and split
 * deterministically across MPI ranks; `seed` must be
 * non-zero.
 */

struct qdrift_data {
	size_t samples;		/* number of independent samples */
	size_t depth;		/* terms drawn per sample */
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
	struct phase2_step_writer *sw;
};

/* Adopt the Hamiltonian + state-prep into a fresh
 * qdrift context, build the importance-sampling CDF
 * from |c_k|, and seed the per-context PRNG from
 * `dt->seed`.  Ownership of `hm` and `*sp_data`
 * transfers (see circ_init).  `sw` is held by
 * reference; its lifetime must outlive qdrift_simul.
 * Returns 0 on success, -1 on error. */
int qdrift_init(struct qdrift *qd, const struct qdrift_data *dt,
	struct circ_hamil hm, enum stprep_kind sp_kind, const void *sp_data,
	struct phase2_step_writer *sw);

/* Release all resources held by `qd`. */
void qdrift_free(struct qdrift *qd);

/* Run `dt.samples` independent qDRIFT samples.  Each
 * sample's overlap is stored in `ct.vals` and, when
 * `qd->sw` is non-NULL, forwarded through the step
 * writer.  Returns 0 on success, -1 on error. */
int qdrift_simul(struct qdrift *qd);

#endif // QDRIFT_H
