#define LOG_SUBSYS "trott"

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "log.h"
#include "phase2.h"

#include "circ/trott.h"

int trott_init(struct trott *tt, const struct trott_data *dt, const data_id fid)
{
	if (circ_init(&tt->ct, fid, dt->steps) < 0) {
		log_error("trott_init: circ_init failed");
		return -1;
	}

	tt->dt = *dt;

	circ_hamil_sort_lex(&tt->ct.hm);

	if (data_circ_writer_init(fid, DATA_CIRCTROTT, dt->steps, &tt->wr)
		< 0) {
		log_error("trott_init: data_circ_writer_init(%s) failed",
			DATA_CIRCTROTT);
		circ_free(&tt->ct);
		return -1;
	}
	if (data_attr_write_dbl(fid, DATA_CIRCTROTT, DATA_CIRCTROTT_DELTA,
		    dt->delta) < 0) {
		log_error("trott_init: write delta attribute failed");
		data_circ_writer_close(&tt->wr);
		circ_free(&tt->ct);
		return -1;
	}

	log_debug("trott_init: delta=%g steps=%zu", dt->delta, dt->steps);
	return 0;
}

void trott_free(struct trott *tt)
{
	data_circ_writer_close(&tt->wr);
	circ_free(&tt->ct);
}

int trott_simul(struct trott *tt)
{
	struct circ *ct = &tt->ct;
	struct circ_values *vals = &ct->vals;

	struct circ_prog prog;
	circ_prog_init(&prog, vals->len, "step");

	circ_prepst(ct);
	for (size_t i = 0; i < vals->len; i++) {
		log_debug("step %zu/%zu", i + 1, vals->len);
		if (circ_step(ct, &ct->hm, tt->dt.delta) < 0) {
			log_error("trott_simul: step %zu failed", i);
			return -1;
		}
		vals->z[i] = circ_measure(ct);

		if (data_circ_write_step(&tt->wr, i, vals->z[i]) < 0) {
			log_error("trott_simul: write_step %zu failed", i);
			return -1;
		}

		circ_prog_tick(&prog);
		circ_prog_emit(&prog, LOG_SUBSYS);
	}

	return 0;
}
