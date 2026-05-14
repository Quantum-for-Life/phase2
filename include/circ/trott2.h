#ifndef TROTT2_H
#define TROTT2_H

#include "phase2/circ.h"
#include "phase2/data.h"

/*
 * Symmetric (Strang) 2nd-order Trotter product formula:
 *
 *   S_2(delta) = prod_{k=1..K} exp(-i h_k delta/2)
 *              * prod_{k=K..1} exp(-i h_k delta/2)
 *
 * One full step is a forward sweep at delta/2 followed by a
 * reverse sweep at delta/2 — both implemented via the
 * existing circ_step / circ_step_reverse primitives.  Output
 * is the overlap <psi| S_2(delta)^s |psi> for s = 1..steps,
 * written to /circ_trott2/values.
 */

struct trott2_data {
	double delta;	/* step size */
	size_t steps;	/* number of Trotter steps */
};

struct trott2 {
	struct circ ct;
	struct trott2_data dt;
	struct data_circ_writer wr;
};

/* Load Hamiltonian and initial state, allocate the register,
 * sort the Hamiltonian lexicographically, create the
 * /circ_trott2 output group with NaN-padded values dataset
 * and the delta attribute.  Returns 0 on success, -1 on
 * error. */
int trott2_init(struct trott2 *t2, const struct trott2_data *dt, data_id fid);

/* Release all resources held by `t2`. */
void trott2_free(struct trott2 *t2);

/* Run `dt.steps` symmetric Trotter steps; per-step overlap
 * goes to `ct.vals` in memory and to /circ_trott2/values[i]
 * on disk (rank 0 only, fflush per step).  Returns 0 on
 * success, -1 on error. */
int trott2_simul(struct trott2 *t2);

#endif // TROTT2_H
