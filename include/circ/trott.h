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
};

/* Load Hamiltonian and initial state from `fid`, allocate
 * the register, sort the Hamiltonian lexicographically.
 * Returns 0 on success, -1 on error. */
int trott_init(struct trott *tt, const struct trott_data *dt, data_id fid);

/* Release all resources held by `tt`. */
void trott_free(struct trott *tt);

/* Run `dt.steps` Trotter steps and store overlaps in
 * `ct.vals`.  Returns 0 on success, -1 on error. */
int trott_simul(struct trott *tt);

/* Write `delta` attribute and overlap series to the
 * /circ_trott group.  Returns 0 on success, -1 on error. */
int trott_write_res(struct trott *tt, data_id fid);

#endif // TROTT_H
