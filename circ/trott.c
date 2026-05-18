/*
 * circ/trott.c -- 1st-order Lie-Trotter product
 * formula.  See include/circ/trott.h for the public
 * API (struct trott, trott_data, init / simul / free).
 *
 * The simulation loop runs `dt.steps` Trotter steps;
 * each step calls circ_step once (forward sweep over
 * the lex-sorted Hamiltonian at angle `dt.delta`)
 * followed by circ_measure.  The per-step overlap is
 * stored in `ct.vals` and, when `tt->sw` is non-NULL,
 * forwarded through the step writer.
 *
 * Error order O(delta^2) per step.  See doc/phase2.md
 * §5.1 for the algorithmic statement.
 */

#define LOG_SUBSYS "trott"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "log.h"
#include "phase2.h"

#include "circ/trott.h"

int trott_init(struct trott *tt, const struct trott_data *dt,
	struct circ_hamil hm, const enum stprep_kind sp_kind,
	const void *sp_data, struct phase2_step_writer *sw)
{
	if (circ_init(&tt->ct, hm, sp_kind, sp_data, dt->steps) < 0) {
		log_error("trott_init: circ_init failed");
		return -1;
	}

	tt->dt = *dt;
	tt->sw = sw;

	circ_hamil_sort_lex(&tt->ct.hm);

	log_debug("trott_init: delta=%g steps=%zu", dt->delta, dt->steps);
	return 0;
}

void trott_free(struct trott *tt)
{
	circ_free(&tt->ct);
}

int trott_simul(struct trott *tt)
{
	struct circ *ct = &tt->ct;
	struct circ_values *vals = &ct->vals;

	circ_prepst(ct);
	for (size_t i = 0; i < vals->len; i++) {
		log_debug("step %zu/%zu", i + 1, vals->len);
		if (circ_step(ct, &ct->hm, tt->dt.delta) < 0) {
			log_error("trott_simul: step %zu failed", i);
			return -1;
		}
		vals->z[i] = circ_measure(ct);

		if (tt->sw && tt->sw->write(tt->sw->ctx, i, vals->z[i]) < 0) {
			log_error("trott_simul: write_step %zu failed", i);
			return -1;
		}
	}

	return 0;
}
