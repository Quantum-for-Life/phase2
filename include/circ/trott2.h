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
};

/* Load Hamiltonian and initial state from `fid`, allocate
 * the register, sort the Hamiltonian lexicographically.
 * Returns 0 on success, -1 on error. */
int trott2_init(struct trott2 *t2, const struct trott2_data *dt, data_id fid);

/* Release all resources held by `t2`. */
void trott2_free(struct trott2 *t2);

/* Run `dt.steps` symmetric Trotter steps and store overlaps
 * in `ct.vals`.  Returns 0 on success, -1 on error. */
int trott2_simul(struct trott2 *t2);

/* Write `delta` attribute and overlap series to the
 * /circ_trott2 group.  Returns 0 on success, -1 on error. */
int trott2_write_res(struct trott2 *t2, data_id fid);

#endif // TROTT2_H
