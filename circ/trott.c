#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "container_of.h"
#include "phase2.h"

#include "circ/trott.h"

int trott_write_res(struct circ *ct, data_id fid);
int trott_simulate(struct circ *ct);

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

static void trott_steps_destroy(struct trott_steps *stp)
{
	free(stp->z);
}

int trott_init(struct trott *tt, const struct trott_data *dt, const data_id fid)
{
	struct circ *c = &tt->ct;
	if (circ_init(c, fid) < 0)
		goto err_circ_init;
	c->simulate = trott_simulate;
	c->write_res = trott_write_res;

	tt->dt = *dt;

	circ_hamil_sort_lex(&tt->ct.hamil);

	if (trott_steps_init(&tt->stp, dt->steps) < 0)
		goto err_trott_res_init;

	return 0;

	// trott_steps_destroy(tt);
err_trott_res_init:
	circ_destroy(&tt->ct);
err_circ_init:
	return -1;
}

void trott_destroy(struct trott *tt)
{
	circ_destroy(&tt->ct);
	trott_steps_destroy(&tt->stp);
}

static int trott_prepst(struct trott *tt)
{
	const struct circ_muldet *md = &tt->ct.muldet;

	qreg_zero(&tt->ct.reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&tt->ct.reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void trott_flush(struct paulis code_hi, struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct trott *tt = data;
	qreg_paulirot(&tt->ct.reg, code_hi, codes_lo, phis, ncodes);
}

static int trott_step(struct trott *tt, const double omega)
{
	const struct circ_hamil *hamil = &tt->ct.hamil;
	struct circ_cache *cache = &tt->ct.cache;

	for (size_t i = 0; i < hamil->nterms; i++) {
		const double phi = omega * hamil->terms[i].cf;
		const struct paulis code = hamil->terms[i].op;

		if (circ_cache_insert(cache, code, phi) == 0)
			continue;

		log_trace("paulirot, term: %zu, num_codes: %zu", i, cache->n);
		circ_cache_flush(cache, trott_flush, tt);
		if (circ_cache_insert(cache, code, phi) < 0)
			return -1;
	}
	log_trace("paulirot, last term group, num_codes: %zu", cache->n);
	circ_cache_flush(cache, trott_flush, tt);

	return 0;
}

static int trott_effect(struct trott *tt)
{
	const double delta = tt->dt.delta;
	if (isnan(delta))
		return -1;
	if (fabs(delta) < DBL_EPSILON)
		return 0;

	return trott_step(tt, delta);
}

static _Complex double trott_measure(struct trott *tt)
{
	const struct circ_muldet *md = &tt->ct.muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&tt->ct.reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}

	return pr;
}

int trott_simulate(struct circ *ct)
{
	int rt = -1;

	size_t prog_pc = 0;
	struct trott *tt = container_of(ct, struct trott, ct);

	trott_prepst(tt);
	for (size_t i = 0; i < tt->stp.len; i++) {
		size_t pc = i * 100 / tt->stp.len;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (trott_step: %zu)", pc, i);
		}

		if (trott_effect(tt) < 0)
			goto ex_trott_effect;
		tt->stp.z[i] = trott_measure(tt);
	}

	rt = 0; /* Success. */
ex_trott_effect:
	return rt;
}

int trott_write_res(struct circ *ct, data_id fid)
{
	int rt = -1;
	struct trott *tt = container_of(ct, struct trott, ct);

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
