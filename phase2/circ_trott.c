#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/circ_trott.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#include "circ.h"

#define MAX_CACHE_CODES UINT64_C(0x10000)

struct trott {
	struct qreg reg;

	_Complex double prod;

	struct code_cache {
		struct paulis code_hi;
		struct paulis *codes_lo;
		double *angles;
		size_t ncodes;
	} cache;
};

static int trott_init(struct trott *tt, struct circ *c)
{
	if (qreg_init(&tt->reg, c->hamil.nqb) < 0)
		goto qreg_init;

	struct paulis *codes_lo =
		malloc(sizeof(struct paulis) * MAX_CACHE_CODES);
	if (!codes_lo)
		goto malloc_codes_lo;
	double *angles = malloc(sizeof(double) * MAX_CACHE_CODES);
	if (!angles)
		goto malloc_angles;

	tt->cache.codes_lo = codes_lo;
	tt->cache.angles = angles;

	return 0;

	// free(angles);
malloc_angles:
	free(codes_lo);
malloc_codes_lo:
	qreg_destroy(&tt->reg);
qreg_init:
	return -1;
}

static void trott_destroy(struct trott *tt)
{
	qreg_destroy(&tt->reg);
	free(tt->cache.codes_lo);
	free(tt->cache.angles);
}

static int trott_prepst(struct trott *tt, struct circ *c)
{
	const struct circ_muldet *md = &c->muldet;

	qreg_zero(&tt->reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&tt->reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void trott_step(struct trott *tt, struct circ *c, const double omega)
{
	const struct circ_hamil *hamil = &c->hamil;

	struct code_cache *cache = &tt->cache;
	cache->ncodes = 0;
	for (size_t i = 0; i < hamil->nterms; i++) {
		const double angle = omega * hamil->terms[i].cf;
		const struct paulis code = hamil->terms[i].op;

		struct paulis code_hi, code_lo;
		paulis_split(
			code, tt->reg.qb_lo, tt->reg.qb_hi, &code_lo, &code_hi);

		if (cache->ncodes == 0) {
			cache->code_hi = code_hi;
			cache->codes_lo[0] = code_lo;
			cache->angles[0] = angle;
			cache->ncodes++;
			continue;
		}

		if (paulis_eq(cache->code_hi, code_hi) &&
			cache->ncodes < MAX_CACHE_CODES) {
			const size_t k = cache->ncodes++;
			cache->codes_lo[k] = code_lo;
			cache->angles[k] = angle;
			continue;
		}

		/* Flush the cache. */
		log_trace("paulirot, term: %zu, num_codes: %zu", i,
			cache->ncodes);
		qreg_paulirot(&tt->reg, cache->code_hi, cache->codes_lo,
			cache->angles, cache->ncodes);

		cache->ncodes = 1;
		cache->code_hi = code_hi;
		cache->codes_lo[0] = code_lo;
		cache->angles[0] = angle;
	}

	log_trace("paulirot, last term group, num_codes: %zu", cache->ncodes);

	if (cache->ncodes > 0)
		qreg_paulirot(&tt->reg, cache->code_hi, cache->codes_lo,
			cache->angles, cache->ncodes);
}

static int trott_effect(struct trott *tt, struct circ *c)
{
	const struct circ_trott_data *data = c->data;
	const double delta = data->delta;
	if (isnan(delta))
		return -1;
	if (fabs(delta) < DBL_EPSILON)
		return 0;
	trott_step(tt, c, delta);

	return 0;
}

static int trott_measure(struct trott *tt, struct circ *c)
{
	const struct circ_muldet *md = &c->muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&tt->reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}
	tt->prod = pr;

	return 0;
}

int circ_res_init(struct circ *c)
{
	struct circ_trott_data *data = c->data;
	size_t nsteps = data->nsteps;

	struct circ_trott_res *res = malloc(sizeof(struct circ_trott_res));
	if (!res)
		goto malloc_res;
	_Complex double *steps = malloc(sizeof(_Complex double) * nsteps);
	if (!steps)
		goto malloc_steps;

	res->steps = steps;
	res->nsteps = nsteps;
	c->res = res;

	return 0;

	// free(steps);
malloc_steps:
	free(res);
malloc_res:
	return -1;
}

void circ_res_destroy(struct circ *c)
{
	struct circ_trott_res *res = c->res;

	free(res->steps);
	free(c->res);
}

int circ_res_write(struct circ *c, data_id fid)
{
	struct circ_trott_res *res = c->res;

	return data_write_vals(fid, DATA_CIRCTROTT, DATA_CIRCTROTT_VALUES,
		res->steps, res->nsteps);
}

int circ_simulate(struct circ *c)
{
	int rt = -1;

	size_t prog_pc = 0;

	struct trott tt;
	if (trott_init(&tt, c) < 0)
		goto trott_init;

	trott_prepst(&tt, c);

	struct circ_trott_res *res = c->res;
	for (size_t i = 0; i < res->nsteps; i++) {
		size_t pc = i * 100 / res->nsteps;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (trott_step: %zu)", pc, i);
		}

		if (trott_effect(&tt, c) < 0)
			goto trott_effect;
		trott_measure(&tt, c);
		res->steps[i] = tt.prod;
	}

	rt = 0; /* Success. */

trott_effect:
	trott_destroy(&tt);
trott_init:

	return rt;
}
