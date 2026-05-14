#ifndef TROTT_H
#define TROTT_H

#include "phase2/circ.h"
#include "phase2/data.h"

/*
 * 1st-order Lie-Trotter product formula:
 *
 *   T_1(delta) = prod_{k=1..K} exp(-i h_k delta)
 *
 * One step applies all Hamiltonian terms once in
 * lexicographic order.  Output is the overlap
 * <psi| T_1(delta)^s |psi> for s = 1..steps,
 * written to /circ_trott/values.
 */

struct trott_data {
	double delta;	/* step size */
	size_t steps;	/* number of Trotter steps */
};

struct trott {
	struct circ ct;
	struct trott_data dt;
	struct data_circ_writer wr;
};

/* Load Hamiltonian and initial state, allocate the register,
 * sort the Hamiltonian lexicographically, create the
 * /circ_trott output group with NaN-padded values dataset
 * and the delta attribute.  Caches the open dataset in
 * tt->wr so trott_simul() can write one row per step.
 * Returns 0 on success, -1 on error. */
int trott_init(struct trott *tt, const struct trott_data *dt, data_id fid);

/* Release all resources held by `tt`. */
void trott_free(struct trott *tt);

/* Run `dt.steps` Trotter steps.  Each step's overlap is
 * stored in `ct.vals` for in-memory use and atomically
 * written (rank-0-only) to /circ_trott/values[i].  Returns
 * 0 on success, -1 on error. */
int trott_simul(struct trott *tt);

#endif // TROTT_H
