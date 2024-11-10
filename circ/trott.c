#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2.h"

#include "circ/trott.h"

static int trott_steps_init(struct trott_steps *stp, size_t steps)
{
	_Complex double *z = malloc(sizeof(_Complex double) * steps);
	if (!z)
		goto malloc_z;

	stp->z = z;
	stp->len = steps;

	return 0;

	// free(steps);
malloc_z:
	return -1;
}

static void trott_steps_free(struct trott_steps *stp)
{
	free(stp->z);
}

int trott_init(struct trott *tt, const struct trott_data *dt, const data_id fid)
{
	if (circ_init(&tt->ct, fid) < 0)
		goto err_circ_init;

	tt->dt = *dt;

	circ_hamil_sort_lex(&tt->ct.hm);

	if (trott_steps_init(&tt->stp, dt->steps) < 0)
		goto err_trott_res_init;

	return 0;

	// trott_steps_destroy(tt);
err_trott_res_init:
	circ_free(&tt->ct);
err_circ_init:
	return -1;
}

void trott_free(struct trott *tt)
{
	circ_free(&tt->ct);
	trott_steps_free(&tt->stp);
}

int trott_simul(struct trott *tt)
{
	struct circ *ct = &tt->ct;
	struct circ_prog prog;

	circ_prog_init(&prog, tt->stp.len);

	circ_prepst(ct);
	for (size_t i = 0; i < tt->stp.len; i++) {
		if (circ_step(ct, &ct->hm, tt->dt.delta) < 0)
			return -1;
		tt->stp.z[i] = circ_measure(ct);

		circ_prog_tick(&prog);
	}

	return 0;
}

int trott_write_res(struct trott *tt, data_id fid)
{
	int rt = -1;

	if (data_grp_create(fid, DATA_CIRCTROTT) < 0)
		goto data_grp_create;
	if (data_attr_write_dbl(fid, DATA_CIRCTROTT, DATA_CIRCTROTT_DELTA,
		    tt->dt.delta) < 0)
		goto data_attr_write;
	if (data_res_write(fid, DATA_CIRCTROTT, DATA_CIRCTROTT_VALUES,
		    tt->stp.z, tt->stp.len) < 0)
		goto data_res_write;

	rt = 0;
data_res_write:
data_attr_write:
data_grp_create:
	return rt;
}
