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

#define MAX_CACHE_CODES UINT64_C(0x0400)

struct circ_cache {
	uint32_t qb_lo, qb_hi;
	struct paulis *codes_lo, code_hi;
	double *phis;
	size_t n;
};

int circ_cache_init(struct circ_cache *ch, size_t qb_lo, size_t qb_hi)
{
	struct paulis *lo =
		malloc(sizeof(struct paulis) * MAX_CACHE_CODES);
	if (!lo)
		goto err_lo;
	double *angles = malloc(sizeof(double) * MAX_CACHE_CODES);
	if (!angles)
		goto err_angles;

	ch->codes_lo = lo;
	ch->phis = angles;
	ch->qb_lo = qb_lo;
	ch->qb_hi = qb_hi;
	ch->n = 0;

	return 0;

	// free(phs);
err_angles:
	free(lo);
err_lo:
	return -1;
}

void circ_cache_destroy(struct circ_cache *ch)
{
	free(ch->codes_lo);
	free(ch->phis);
}

int circ_cache_insert(struct circ_cache *ch, struct paulis code, double angle)
{
	struct paulis lo, hi;
	paulis_split(code, ch->qb_lo, ch->qb_hi, &lo, &hi);
	if (ch->n == 0) {
		ch->code_hi = hi;
		ch->codes_lo[0] = lo;
		ch->phis[0] = angle;
		ch->n = 1;
		return 0;
	}

	if (ch->n < MAX_CACHE_CODES && paulis_eq(ch->code_hi, hi)) {
		const size_t k = ch->n++;
		ch->codes_lo[k] = lo;
		ch->phis[k] = angle;
		return 0;
	}

	return -1;
}

void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data)
{
	if (ch->n > 0 && op)
		op(ch->code_hi, ch->codes_lo, ch->phis, ch->n, data);
	ch->n = 0;
}

struct trott {
	struct qreg reg;
	_Complex double prod;
	struct circ_cache cache;
};

static int trott_init(struct trott *tt, struct circ *c)
{
	if (qreg_init(&tt->reg, c->hamil.nqb) < 0)
		goto err_qreg_init;
	if (circ_cache_init(&tt->cache, tt->reg.qb_lo, tt->reg.qb_hi) < 0)
		goto err_circ_cache_init;

	return 0;

	// circ_cache_destroy(&tt->cache);
err_circ_cache_init:
	qreg_destroy(&tt->reg);
err_qreg_init:
	return -1;
}

static void trott_destroy(struct trott *tt)
{
	qreg_destroy(&tt->reg);
	circ_cache_destroy(&tt->cache);
}

static int trott_prepst(struct trott *tt, struct circ *c)
{
	const struct circ_muldet *md = &c->muldet;

	qreg_zero(&tt->reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&tt->reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void trott_flush(struct paulis code_hi, struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct trott *tt = data;
	qreg_paulirot(&tt->reg, code_hi, codes_lo, phis, ncodes);
}

static int trott_step(struct trott *tt, struct circ *c, const double omega)
{
	const struct circ_hamil *hamil = &c->hamil;
	struct circ_cache *cache = &tt->cache;

	circ_cache_flush(cache, nullptr, nullptr);
	for (size_t i = 0; i < hamil->nterms; i++) {
		const double phi = omega * hamil->terms[i].cf;
		const struct paulis code = hamil->terms[i].op;

		if (circ_cache_insert(cache, code, phi) == 0)
			continue;

		/* Flush the cache. */
		log_trace("paulirot, term: %zu, num_codes: %zu", i,
			cache->n);
		circ_cache_flush(cache, trott_flush, tt);

		if (circ_cache_insert(cache, code, phi) < 0)
			return -1;
	}

	log_trace("paulirot, last term group, num_codes: %zu", cache->n);

	circ_cache_flush(cache, trott_flush, tt);
	return 0;
}

static int trott_effect(struct trott *tt, struct circ *c)
{
	const struct circ_trott_data *data = c->data;
	const double delta = data->delta;
	if (isnan(delta))
		return -1;
	if (fabs(delta) < DBL_EPSILON)
		return 0;

	return trott_step(tt, c, delta);
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
		goto ex_trott_init;

	trott_prepst(&tt, c);

	struct circ_trott_res *res = c->res;
	for (size_t i = 0; i < res->nsteps; i++) {
		size_t pc = i * 100 / res->nsteps;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (trott_step: %zu)", pc, i);
		}

		if (trott_effect(&tt, c) < 0)
			goto ex_trott_effect;
		trott_measure(&tt, c);
		res->steps[i] = tt.prod;
	}

	rt = 0; /* Success. */
ex_trott_effect:
	trott_destroy(&tt);
ex_trott_init:
	return rt;
}
