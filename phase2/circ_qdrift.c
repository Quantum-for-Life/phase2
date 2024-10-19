#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/circ_qdrift.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#define SEED UINT64_C(0xeccd9dcc749fcdca)

#define MAX_CACHE_CODES UINT64_C(0x10000)

struct qdrift {
	struct qreg reg;

	_Complex double prod;

	struct code_cache {
		struct paulis code_hi;
		struct paulis *codes_lo;
		double *angles;
		size_t ncodes;
	} cache;

	struct xoshiro256ss rng;
	size_t *smpl;
};

static int qdrift_init(struct qdrift qd, struct circ *c)
{
	if (qreg_init(&qd->reg, c->hamil.nqb) < 0)
		goto err_qreg_init;

	struct paulis *codes_lo =
		malloc(sizeof(struct paulis) * MAX_CACHE_CODES);
	if (!codes_lo)
		goto err_codes_lo;
	double *angles = malloc(sizeof(double) * MAX_CACHE_CODES);
	if (!angles)
		goto err_angles;
	qd->cache.codes_lo = codes_lo;
	qd->cache.angles = angles;

	xoshiro256ss_init(&qd->rng, SEED);

	struct circ_qdrift_data *data = c->data;
	size_t *smpl = malloc(sizeof(size_t) * data->depth);
	if (!smpl)
		goto err_smpl;
	qd->smpl = smpl;

	return 0;

	// free(smpl);
err_smpl:
	free(angles);
err_angles:
	free(codes_lo);
err_codes_lo:
err_qreg_init:
	return 1;
}

static void qdrift_destroy(struct qdrift *qd)
{
	qreg_destroy(&qd->reg);
	free(qd->smpl);
	free(qd->angles);
	free(qd->codes_lo);
}

static int qdrift_prepst(struct qdrift *qd, struct circ *c)
{
	const struct circ_muldet *md = &c->muldet;

	qreg_zero(&qd->reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&c->reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void qdrift_step(struct qdrift *qd, struct circ *c, const double omega)
{
	const struct circ_hamil *hamil = &c->hamil;
	const struct circ_qdrift_data *data = c->data;

	struct code_cache *cache = &qd->cache;
	cache->ncodes = 0;
	for (size_t i = 0; i < data->depth; i++) {
		const double angle = omega;
		const size_t i_smpl = qd->smpl[i];
		const struct paulis code = hamil->terms[i_smpl].op;

		struct paulis code_hi, code_lo;
		paulis_split(code, qd->reg.nqb_lo, qd->reg.nqb_hi, &code_lo,
			&code_hi);

		if (cache.num_codes == 0) {
			cache->code_hi = code_hi;
			cache->codes_lo[0] = code_lo;
			cache->angles[0] = angle;
			cache->num_codes++;
			continue;
		}

		if (paulis_eq(cache->code_hi, code_hi) &&
			cache->ncodes < MAX_CACHE_CODES) {
			const size_t k = cache->ncodes++;
			cache->codes_lo[k] = code_lo;
			cache->angles[k] = angle;
			continue;
		}

		log_trace("paulirot, term: %zu, num_codes: %zu", i,
			cache->ncodes);
		qreg_paulirot(&qd->reg, cache->code_hi, cache->codes_lo,
			cache->angles, cache->ncodes);

		cache->ncodes = 1;
		cache->code_hi = code_hi;
		cache->codes_lo[0] = code_lo;
		cache->angles[0] = angle;
	}

	log_trace("paulirot, last term group, num_codes: %zu", cache->ncodes);

	if (cache->ncodes > 0)
		qreg_paulirot(&qd->reg, cache->code_hi, cache->codes_lo,
			cache->angles, cache->ncodes);
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
	qdrift_step(c, theta);

	return 0;
}

static int qdrift_measure(struct qdrift *qd, struct circ *c)
{
	const struct circ_muldet *md = &c->muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&c->reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}
	qd->prod = pr;

	return 0;
}

static size_t sample_invcdf(struct qdrift *qd, struct circ *c, double x)
{
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
		double x = rand_dbl01(&qd->rng)
				   c->smpl[i] = sample_invcdf(qd, c, x);
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
	struct circ_qdrift_res *res = c->res;

	free(res->samples);
	free(c->res);
}

int circ_res_write(struct circ *c, data_id fid)
{
	struct circ_qdrift_res *res = c->res;

	return data_circ_qdrift_write_values(fid, res->samples, res->nsamples);
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

		sample_terms(&qd, &c);
		qdrift_prepst(&qd, &c);
		if (qdrift_effect(&qd, &c) < 0)
			goto ex_qdrift_effect;
		qdrift_measure(&qd, &c);
		res->samples[i] = qd.prod;
	}

	rt = 0; /* Success. */

ex_qdrift_effect:
	qdrift_destroy(&qd);
ex_qdrift_init:

	return rt;
}
