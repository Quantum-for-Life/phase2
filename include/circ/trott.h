#ifndef TROTT_H
#define TROTT_H

#include "phase2/circ.h"
#include "phase2/step_writer.h"

/*
 * 1st-order Lie-Trotter product formula:
 *
 *   T_1(delta) = prod_{k=1..K} exp(-i h_k delta)
 *
 * One step applies all Hamiltonian terms once in
 * lexicographic order.  Output is the overlap
 * <psi| T_1(delta)^s |psi> for s = 1..steps.
 */

struct trott_data {
	size_t steps;	/* number of Trotter steps */
	double delta;	/* step size */
};

struct trott {
	struct circ ct;
	struct trott_data dt;
	struct phase2_step_writer *sw;	/* per-step output sink;
					 * may be NULL */
};

/* Adopt hm and *sp_data (ownership transfers, see
 * circ_init); sw is held by reference, must outlive
 * trott_simul.  Returns 0 or -1. */
int trott_init(struct trott *tt, const struct trott_data *dt,
	struct circ_hamil hm, enum stprep_kind sp_kind,
	const void *sp_data, struct phase2_step_writer *sw);

void trott_free(struct trott *tt);

/* Run dt.steps steps; overlap per step into ct.vals
 * and, if non-NULL, tt->sw.  Returns 0 or -1. */
int trott_simul(struct trott *tt);

#endif // TROTT_H
