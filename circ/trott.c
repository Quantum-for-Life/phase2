#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2.h"

#include "circ/trott.h"

int trott_init(struct trott *tt, const struct trott_data *dt, const data_id fid)
{
	if (circ_init(&tt->ct, fid, dt->steps) < 0)
		return -1;

	tt->dt = *dt;

	circ_hamil_sort_lex(&tt->ct.hm);

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

	struct circ_prog prog;
	circ_prog_init(&prog, vals->len);

	circ_prepst(ct);
	for (size_t i = 0; i < vals->len; i++) {
		if (circ_step(ct, &ct->hm, tt->dt.delta) < 0)
			return -1;
		vals->z[i] = circ_measure(ct);

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
		    tt->ct.vals.z, tt->ct.vals.len) < 0)
		goto data_res_write;

	rt = 0;
data_res_write:
data_attr_write:
data_grp_create:
	return rt;
}
