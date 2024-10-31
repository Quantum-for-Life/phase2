#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "circ/qdrift.h"
#include "container_of.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#define SEED UINT64_C(0xeccd9dcc749fcdca)

int qdrift_write_res(struct circ *c, data_id fid);
int qdrift_simulate(struct circ *c);

static int qdrift_res_init(struct circ_qdrift *qd, size_t nsamples)
{
	_Complex double *samples = malloc(sizeof(_Complex double) * nsamples);
	if (!samples)
		goto err_samples;
	qd->res.samples = samples;
	qd->res.nsamples = nsamples;

	return 0;

// free(samples);
err_samples:
	return -1;
}

static void qdrift_res_destroy(struct circ_qdrift *qd)
{
	free(qd->res.samples);
}

int circ_qdrift_init(
	struct circ_qdrift *qd, struct circ_qdrift_data *data, data_id fid)
{
	struct circ *c = &qd->circ;
	if (circ_init(c, fid) < 0)
		goto err_circ_init;
	c->simulate = qdrift_simulate;
	c->write_res = qdrift_write_res;

	xoshiro256ss_init(&qd->rng, SEED);
	size_t *smpl = malloc(sizeof(size_t) * data->depth);
	if (!smpl)
		goto err_smpl;
	qd->smpl = smpl;

	if (qdrift_res_init(qd, data->nsamples) < 0)
		goto err_res_init;

	qd->depth = data->depth;
	qd->step_size = data->step_size;

	return 0;

err_res_init:
	free(smpl);
err_smpl:
	circ_destroy(c);
err_circ_init:
	return -1;
}

void circ_qdrift_destroy(struct circ_qdrift *qd)
{
	circ_destroy(&qd->circ);
	qdrift_res_destroy(qd);
}

static int qdrift_prepst(struct circ_qdrift *qd)
{
	const struct circ_muldet *md = &qd->circ.muldet;

	qreg_zero(&qd->circ.reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&qd->circ.reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void qdrift_flush(struct paulis code_hi, struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct circ_qdrift *qd = data;
	qreg_paulirot(&qd->circ.reg, code_hi, codes_lo, phis, ncodes);
}

static int qdrift_step(struct circ_qdrift *qd, const double omega)
{
	const struct circ_hamil *hamil = &qd->circ.hamil;
	struct circ_cache *cache = &qd->circ.cache;

	for (size_t i = 0; i < qd->depth; i++) {
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

static int qdrift_effect(struct circ_qdrift *qd)
{
	const double t = qd->step_size;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;
	const double theta = asin(t);

	return qdrift_step(qd, theta);
}

static _Complex double qdrift_measure(struct circ_qdrift *qd)
{
	const struct circ_muldet *md = &qd->circ.muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&qd->circ.reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}

	return pr;
}

static size_t sample_invcdf(struct circ_qdrift *qd, double x)
{
	(void)qd;
	size_t i = 0;
	double cdf = 0;
	while (cdf <= x)
		cdf += fabs(qd->circ.hamil.terms[i++].cf);
	return i - 1; /* Never again make the same off-by-one error! */
}

/* TODO: Move to xoshiro256ss.h and document. */
#define rand_dbl01(rng) ((double)(xoshiro256ss_next(rng) >> 11) * 0x1.0p-53)

static void sample_terms(struct circ_qdrift *qd)
{
	for (size_t i = 0; i < qd->depth; i++) {
		double x = rand_dbl01(&qd->rng);
		qd->smpl[i] = sample_invcdf(qd, x);
	}
}

int qdrift_write_res(struct circ *c, data_id fid)
{
	int rt = -1;

	struct circ_qdrift *qd = container_of(c, struct circ_qdrift, circ);

	if (data_grp_create(fid, DATA_CIRCQDRIFT) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_STEPSIZE,
		    qd->step_size) < 0)
		goto data_res_write;
	if (data_attr_write(
		    fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_DEPTH, qd->depth) < 0)
		goto data_res_write;
	if (data_res_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_VALUES,
		    qd->res.samples, qd->res.nsamples) < 0)
		goto data_res_write;

	rt = 0;

data_res_write:
	return rt;
}

int qdrift_simulate(struct circ *c)
{
	int rt = -1;

	size_t prog_pc = 0;
	struct circ_qdrift *qd = container_of(c, struct circ_qdrift, circ);

	for (size_t i = 0; i < qd->res.nsamples; i++) {
		size_t pc = i * 100 / qd->res.nsamples;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (samples: %zu)", pc, i);
		}

		sample_terms(qd);
		qdrift_prepst(qd);
		if (qdrift_effect(qd) < 0)
			goto ex_qdrift_effect;
		qd->res.samples[i] = qdrift_measure(qd);
	}

	rt = 0; /* Success. */
ex_qdrift_effect:
	return rt;
}
