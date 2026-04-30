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
 * One full step is forward sweep at delta/2 followed by
 * reverse sweep at delta/2 — both implemented via existing
 * circ_step / circ_step_reverse primitives.
 */

struct trott2_data {
	double delta;
	size_t steps;
};

struct trott2 {
	struct circ ct;
	struct trott2_data dt;
};

int trott2_init(struct trott2 *t2, const struct trott2_data *dt, data_id fid);

void trott2_free(struct trott2 *t2);

int trott2_simul(struct trott2 *t2);

int trott2_write_res(struct trott2 *t2, data_id fid);

#endif // TROTT2_H
