#define LOG_SUBSYS "trott2"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "log.h"
#include "phase2.h"

#include "circ/trott2.h"

int trott2_init(struct trott2 *t2, const struct trott2_data *dt, const data_id fid)
{
	if (circ_init(&t2->ct, fid, dt->steps) < 0) {
		log_error("trott2_init: circ_init failed");
		return -1;
	}

	t2->dt = *dt;

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

	struct circ_prog prog;
	circ_prog_init(&prog, vals->len, "step");

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

		circ_prog_tick(&prog);
		circ_prog_emit(&prog, LOG_SUBSYS);
	}

	return 0;
}

int trott2_write_res(struct trott2 *t2, data_id fid)
{
	int rt = -1;

	if (data_grp_create(fid, DATA_CIRCTROTT2) < 0)
		goto data_grp_create;
	if (data_attr_write_dbl(fid, DATA_CIRCTROTT2, DATA_CIRCTROTT2_DELTA,
		    t2->dt.delta) < 0)
		goto data_attr_write;
	if (data_res_write(fid, DATA_CIRCTROTT2, DATA_CIRCTROTT2_VALUES,
		    t2->ct.vals.z, t2->ct.vals.len) < 0)
		goto data_res_write;

	rt = 0;
data_res_write:
data_attr_write:
data_grp_create:
	return rt;
}
