/*
 * circ/trott2.c -- 2nd-order symmetric (Strang)
 * Trotter product formula.  See
 * include/circ/trott2.h for the public API.
 *
 * Each Trotter step is a forward sweep at `delta/2`
 * (circ_step) followed by a reverse sweep at the same
 * angle (circ_step_reverse).  The reverse sweep
 * traverses the lex-sorted Hamiltonian in reverse
 * order; the two sweeps together yield the symmetric
 * S_2 integrator.  Per-step overlap is stored in
 * `ct.vals` and forwarded through `t2->sw` when set.
 *
 * Error order O(delta^3) per step.  See doc/phase2.md
 * §5.2 for the algorithmic statement.
 */

#define LOG_SUBSYS "trott2"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "log.h"
#include "phase2.h"

#include "circ/trott2.h"

int trott2_init(struct trott2 *t2, const struct trott2_data *dt,
	struct circ_hamil hm, const enum stprep_kind sp_kind,
	const void *sp_data, struct phase2_step_writer *sw)
{
	if (circ_init(&t2->ct, hm, sp_kind, sp_data, dt->steps) < 0) {
		log_error("trott2_init: circ_init failed");
		return -1;
	}

	t2->dt = *dt;
	t2->sw = sw;

	circ_hamil_sort_lex(&t2->ct.hm);

	log_debug("trott2_init: delta=%g steps=%zu", dt->delta, dt->steps);
	return 0;
}

void trott2_free(struct trott2 *t2)
{
	circ_free(&t2->ct);
}

int trott2_simul(struct trott2 *t2)
{
	struct circ *ct = &t2->ct;
	struct circ_values *vals = &ct->vals;

	const double half = t2->dt.delta / 2.0;

	circ_prepst(ct);
	for (size_t i = 0; i < vals->len; i++) {
		log_debug("step %zu/%zu fwd", i + 1, vals->len);
		if (circ_step(ct, &ct->hm, half) < 0) {
			log_error("trott2_simul: fwd sweep failed at step %zu",
				i);
			return -1;
		}
		log_debug("step %zu/%zu rev", i + 1, vals->len);
		if (circ_step_reverse(ct, &ct->hm, half) < 0) {
			log_error("trott2_simul: rev sweep failed at step %zu",
				i);
			return -1;
		}
		vals->z[i] = circ_measure(ct);

		if (t2->sw && t2->sw->write(t2->sw->ctx, i, vals->z[i]) < 0) {
			log_error("trott2_simul: write_step %zu failed", i);
			return -1;
		}
	}

	return 0;
}
