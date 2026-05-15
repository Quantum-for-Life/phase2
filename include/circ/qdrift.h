#ifndef QDRIFT_H
#define QDRIFT_H

#include <stddef.h>

#include "phase2/circ.h"
#include "phase2/prob.h"
#include "phase2/step_writer.h"
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
	struct phase2_step_writer *sw;
};

int qdrift_init(struct qdrift *qd, const struct qdrift_data *dt,
	struct circ_hamil hm, enum stprep_kind sp_kind, const void *sp_data,
	struct phase2_step_writer *sw);

void qdrift_free(struct qdrift *qd);

int qdrift_simul(struct qdrift *qd);

#endif // QDRIFT_H
