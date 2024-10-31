#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "circ/trott.h"
#include "container_of.h"
#include "phase2.h"

int trott_write_res(struct circ *c, data_id fid);
int trott_simulate(struct circ *c);

static int trott_res_init(struct circ_trott *tt, size_t nsteps)
{
	_Complex double *steps = malloc(sizeof(_Complex double) * nsteps);
	if (!steps)
		goto malloc_steps;

	tt->res.steps = steps;
	tt->res.nsteps = nsteps;

	return 0;

	// free(steps);
malloc_steps:
	return -1;
}

static void trott_res_destroy(struct circ_trott *tt)
{
	free(tt->res.steps);
}

int circ_trott_init(
	struct circ_trott *tt, struct circ_trott_data *data, data_id fid)
{
	struct circ *c = &tt->circ;
	if (circ_init(c, fid) < 0)
		goto err_circ_init;
	c->simulate = trott_simulate;
	c->write_res = trott_write_res;

	tt->delta = data->delta;
	if (trott_res_init(tt, data->nsteps) < 0)
		goto err_trott_res_init;

	return 0;

	// trott_res_destroy(tt);
err_trott_res_init:
	circ_destroy(&tt->circ);
err_circ_init:
	return -1;
}

void circ_trott_destroy(struct circ_trott *tt)
{
	circ_destroy(&tt->circ);
	trott_res_destroy(tt);
}

static int trott_prepst(struct circ_trott *tt)
{
	const struct circ_muldet *md = &tt->circ.muldet;

	qreg_zero(&tt->circ.reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&tt->circ.reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void trott_flush(struct paulis code_hi, struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct circ_trott *tt = data;
	qreg_paulirot(&tt->circ.reg, code_hi, codes_lo, phis, ncodes);
}

static int trott_step(struct circ_trott *tt, const double omega)
{
	const struct circ_hamil *hamil = &tt->circ.hamil;
	struct circ_cache *cache = &tt->circ.cache;

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

static int trott_effect(struct circ_trott *tt)
{
	const double delta = tt->delta;
	if (isnan(delta))
		return -1;
	if (fabs(delta) < DBL_EPSILON)
		return 0;

	return trott_step(tt, delta);
}

static _Complex double trott_measure(struct circ_trott *tt)
{
	const struct circ_muldet *md = &tt->circ.muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&tt->circ.reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}

	return pr;
}

int trott_write_res(struct circ *c, data_id fid)
{
	int rt = -1;
	struct circ_trott *tt = container_of(c, struct circ_trott, circ);

	if (data_grp_create(fid, DATA_CIRCTROTT) < 0)
		goto data_grp_create;
	if (data_attr_write_dbl(
		    fid, DATA_CIRCTROTT, DATA_CIRCTROTT_DELTA, tt->delta) < 0)
		goto data_attr_write;
	if (data_res_write(fid, DATA_CIRCTROTT, DATA_CIRCTROTT_VALUES,
		    tt->res.steps, tt->res.nsteps) < 0)
		goto data_res_write;

	rt = 0;
data_res_write:
data_attr_write:
data_grp_create:
	return rt;
}

int trott_simulate(struct circ *c)
{
	int rt = -1;

	size_t prog_pc = 0;
	struct circ_trott *tt = container_of(c, struct circ_trott, circ);

	trott_prepst(tt);
	for (size_t i = 0; i < tt->res.nsteps; i++) {
		size_t pc = i * 100 / tt->res.nsteps;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (trott_step: %zu)", pc, i);
		}

		if (trott_effect(tt) < 0)
			goto ex_trott_effect;
		tt->res.steps[i] = trott_measure(tt);
	}

	rt = 0; /* Success. */
ex_trott_effect:
	return rt;
}
