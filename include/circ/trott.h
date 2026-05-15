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
	double delta;	/* step size */
	size_t steps;	/* number of Trotter steps */
};

struct trott {
	struct circ ct;
	struct trott_data dt;
	struct phase2_step_writer *sw;	/* per-step output sink;
					 * may be NULL */
};

/* Adopt the pre-loaded Hamiltonian + state-prep into a
 * fresh trott context.  Ownership of `hm` and `*sp_data`
 * transfers (see circ_init).  `sw` is held by reference;
 * its lifetime must outlive trott_simul.  Returns 0 on
 * success, -1 on error. */
int trott_init(struct trott *tt, const struct trott_data *dt,
	struct circ_hamil hm, enum stprep_kind sp_kind,
	const void *sp_data, struct phase2_step_writer *sw);

/* Release all resources held by `tt`. */
void trott_free(struct trott *tt);

/* Run `dt.steps` Trotter steps.  Each step's overlap is
 * stored in `ct.vals` for in-memory use and, when `tt->sw`
 * is non-NULL, also forwarded through the step writer.
 * Returns 0 on success, -1 on error. */
int trott_simul(struct trott *tt);

#endif // TROTT_H
