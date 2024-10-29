#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/circ_qdrift.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#define SEED UINT64_C(0xeccd9dcc749fcdca)

struct qdrift {
	struct qreg reg;
	_Complex double prod;
	struct circ_cache cache;

	struct xoshiro256ss rng;
	size_t *smpl;
};

static int qdrift_init(struct qdrift *qd, struct circ *c)
{
	if (qreg_init(&qd->reg, c->hamil.nqb) < 0)
		goto err_qreg_init;
	if (circ_cache_init(&qd->cache, qd->reg.qb_lo, qd->reg.qb_hi) < 0)
		goto err_circ_cache_init;

	xoshiro256ss_init(&qd->rng, SEED);

	struct circ_qdrift_data *data = c->data;
	size_t *smpl = malloc(sizeof(size_t) * data->depth);
	if (!smpl)
		goto err_smpl;
	qd->smpl = smpl;

	return 0;

	// free(smpl);
err_smpl:
	circ_cache_destroy(&qd->cache);
err_circ_cache_init:
	qreg_destroy(&qd->reg);
err_qreg_init:
	return -1;
}

static void qdrift_destroy(struct qdrift *qd)
{
	qreg_destroy(&qd->reg);
	circ_cache_destroy(&qd->cache);
}

static int qdrift_prepst(struct qdrift *qd, struct circ *c)
{
	const struct circ_muldet *md = &c->muldet;

	qreg_zero(&qd->reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&qd->reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void qdrift_flush(struct paulis code_hi, struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct qdrift *qd = data;
	qreg_paulirot(&qd->reg, code_hi, codes_lo, phis, ncodes);
}

static int qdrift_step(struct qdrift *qd, struct circ *c, const double omega)
{
	const struct circ_hamil *hamil = &c->hamil;
	const struct circ_qdrift_data *data = c->data;
	struct circ_cache *cache = &qd->cache;

	for (size_t i = 0; i < data->depth; i++) {
		const double phi = omega;
		const size_t i_smpl = qd->smpl[i];
		const struct paulis code = hamil->terms[i_smpl].op;

		if (circ_cache_insert(cache, code, phi) == 0)
			continue;

		log_trace("paulirot, term: %zu, num_codes: %zu", i, cache->n);
		circ_cache_flush(cache, qdrift_flush, qd);
		if (circ_cache_insert(cache, code, phi) < 0)
			return -1;
	}
	log_trace("paulirot, last term group, num_codes: %zu", cache->n);
	circ_cache_flush(cache, qdrift_flush, qd);

	return 0;
}

static int qdrift_effect(struct qdrift *qd, struct circ *c)
{
	const struct circ_qdrift_data *data = c->data;
	const double t = data->step_size;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;
	const double theta = asin(t);

	return qdrift_step(qd, c, theta);
}

static int qdrift_measure(struct qdrift *qd, struct circ *c)
{
	const struct circ_muldet *md = &c->muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&qd->reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}
	qd->prod = pr;

	return 0;
}

static size_t sample_invcdf(struct qdrift *qd, struct circ *c, double x)
{
	(void)qd;
	size_t i = 0;
	double cdf = 0;
	while (cdf <= x)
		cdf += fabs(c->hamil.terms[i++].cf);
	return i - 1; /* Never again make the same off-by-one error! */
}

/* TODO: Move to xoshiro256ss.h and document. */
#define rand_dbl01(rng) ((double)(xoshiro256ss_next(rng) >> 11) * 0x1.0p-53)

static void sample_terms(struct qdrift *qd, struct circ *c)
{
	const struct circ_qdrift_data *data = c->data;
	for (size_t i = 0; i < data->depth; i++) {
		double x = rand_dbl01(&qd->rng);
		qd->smpl[i] = sample_invcdf(qd, c, x);
	}
}

int circ_res_init(struct circ *c)
{
	struct circ_qdrift_data *data = c->data;
	size_t nsamples = data->nsamples;

	struct circ_qdrift_res *res = malloc(sizeof(struct circ_qdrift_res));
	if (!res)
		goto err_res;

	_Complex double *samples = malloc(sizeof(_Complex double) * nsamples);
	if (!samples)
		goto err_samples;

	res->samples = samples;
	res->nsamples = nsamples;
	c->res = res;

	return 0;

	// free(samples);
err_samples:
	free(res);
err_res:
	return -1;
}

void circ_res_destroy(struct circ *c)
{
	const struct circ_qdrift_res *res = c->res;

	free(res->samples);
	free(c->res);
}

int circ_res_write(struct circ *c, data_id fid)
{
	int rt = -1;

	const struct circ_qdrift_data *data = c->data;
	const struct circ_qdrift_res *res = c->res;

	if (data_grp_create(fid, DATA_CIRCQDRIFT) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_STEPSIZE,
		    data->step_size) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_DEPTH,
		    data->depth) < 0)
		goto data_res_write;
	if (data_res_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_VALUES,
		    res->samples, res->nsamples) < 0)
		goto data_res_write;

	rt = 0;

data_res_write:
	return rt;
}

int circ_simulate(struct circ *c)
{
	int rt = -1;

	size_t prog_pc = 0;

	struct qdrift qd;
	if (qdrift_init(&qd, c) < 0)
		goto ex_qdrift_init;

	struct circ_qdrift_res *res = c->res;
	for (size_t i = 0; i < res->nsamples; i++) {
		size_t pc = i * 100 / res->nsamples;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (samples: %zu)", pc, i);
		}

		sample_terms(&qd, c);
		qdrift_prepst(&qd, c);
		if (qdrift_effect(&qd, c) < 0)
			goto ex_qdrift_effect;
		qdrift_measure(&qd, c);
		res->samples[i] = qd.prod;
	}

	rt = 0; /* Success. */

ex_qdrift_effect:
	qdrift_destroy(&qd);
ex_qdrift_init:

	return rt;
}
