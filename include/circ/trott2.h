#ifndef TROTT2_H
#define TROTT2_H

#include "phase2/circ.h"
#include "phase2/step_writer.h"

/*
 * Symmetric (Strang) 2nd-order Trotter product formula:
 *
 *   S_2(delta) = prod_{k=1..K} exp(-i h_k delta/2)
 *              * prod_{k=K..1} exp(-i h_k delta/2)
 *
 * One full step is a forward sweep at delta/2 followed by a
 * reverse sweep at delta/2 — both implemented via the
 * existing circ_step / circ_step_reverse primitives.  Output
 * is the overlap <psi| S_2(delta)^s |psi> for s = 1..steps.
 */

struct trott2_data {
	size_t steps;	/* number of Trotter steps */
	double delta;	/* step size */
};

struct trott2 {
	struct circ ct;
	struct trott2_data dt;
	struct phase2_step_writer *sw;
};

/* Adopt hm and *sp_data (ownership transfers, see
 * circ_init); sw is held by reference, must outlive
 * trott2_simul.  Returns 0 or -1. */
int trott2_init(struct trott2 *t2, const struct trott2_data *dt,
	struct circ_hamil hm, enum stprep_kind sp_kind,
	const void *sp_data, struct phase2_step_writer *sw);

void trott2_free(struct trott2 *t2);

/* Run dt.steps Strang steps; overlap per step into
 * ct.vals and, if non-NULL, t2->sw.  Returns 0 or -1. */
int trott2_simul(struct trott2 *t2);

#endif // TROTT2_H
