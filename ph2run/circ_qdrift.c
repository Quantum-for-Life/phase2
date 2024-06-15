#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "circ.h"
#include "log.h"
#include "xoshiro256ss.h"

#define MAX_CACHE_CODES (1024)

#define PRNG_SEED (0x235eac32)

struct circ_qdrift {
	size_t	    num_qb;
	struct qreg reg;

	const struct circ_qdrift_data *data;

	double prod[2];

	struct code_cache {
		struct paulis code_hi;
		struct paulis codes_lo[MAX_CACHE_CODES];
		double	      angles[MAX_CACHE_CODES];
		size_t	      num_codes;
	} cache;

	struct xoshiro256ss rng;
	size_t		   *sampled_idx;
};

static int circ_create(struct circ_qdrift *c,
	const struct circ_qdrift_data *data, const size_t num_qubits)
{
	struct qreg reg;
	if (qreg_init(&reg, num_qubits) < 0)
		return -1;

	c->num_qb = num_qubits;
	c->data	  = data;
	c->reg	  = reg;

	xoshiro256ss_init(&c->rng, PRNG_SEED);
	size_t *sampled_idx = malloc(sizeof(size_t) * data->depth);
	if (sampled_idx == NULL)
		return -1;
	c->sampled_idx = sampled_idx;

	return 0;
}

static void circ_destroy(struct circ_qdrift *c)
{
	qreg_destroy(&c->reg);
	if (c->sampled_idx)
		free(c->sampled_idx);
}

int circ_data_from_file(struct circ_qdrift_data *cd, const data_id fid)
{
	int rc = circ_hamil_from_file(&cd->hamil, fid);
	rc |= circ_multidet_from_file(&cd->multidet, fid);
	rc |= data_circ_qdrift_getattrs(
		fid, &cd->num_samples, &cd->step_size, &cd->depth);

	return rc;
}

int circ_qdrift_data_init(struct circ_qdrift_data *cd, data_id fid)
{
	circ_hamil_init(&cd->hamil);
	circ_multidet_init(&cd->multidet);

	if (circ_data_from_file(cd, fid) < 0)
		return -1;
	cd->samples[0] = malloc(sizeof(double) * 2 * cd->num_samples);
	if (cd->samples[0] == NULL)
		return -1;
	cd->samples[1] = cd->samples[0] + cd->num_samples;

	return 0;
}

void circ_qdrift_data_destroy(struct circ_qdrift_data *cd)
{
	circ_multidet_destroy(&cd->multidet);
	circ_hamil_destroy(&cd->hamil);

	free(cd->samples[0]);
}

static int circuit_prepst(struct circ_qdrift *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	qreg_zero(&c->reg);
	for (size_t i = 0; i < md->num_dets; i++) {
		const _Complex double coeff = md->dets[i].coeff[0] +
			 _Complex_I * md->dets[i].coeff[1];
		qreg_setamp(&c->reg, md->dets[i].idx, coeff);
	}

	return 0;
}

static void trott_step(struct circ_qdrift *c, const double omega)
{
	const struct circ_hamil *hamil = &c->data->hamil;

	struct code_cache cache = c->cache;
	cache.num_codes		= 0;

	for (size_t i = 0; i < c->data->depth; i++) {
		const double	    angle     = omega;
		const size_t	    i_sampled = c->sampled_idx[i];
		const struct paulis code      = hamil->paulis[i_sampled];

		struct paulis code_hi, code_lo;
		paulis_split(
			code, c->reg.qb_lo, c->reg.qb_hi, &code_lo, &code_hi);
		paulis_shr(&code_hi, c->reg.qb_lo);

		if (cache.num_codes == 0) {
			cache.code_hi	  = code_hi;
			cache.codes_lo[0] = code_lo;
			cache.angles[0]	  = angle;
			cache.num_codes++;
			continue;
		}

		if (paulis_eq(cache.code_hi, code_hi) &&
			cache.num_codes < MAX_CACHE_CODES) {
			const size_t k	  = cache.num_codes++;
			cache.codes_lo[k] = code_lo;
			cache.angles[k]	  = angle;
			continue;
		}

		log_trace("paulirot, term: %zu, num_codes: %zu", i,
			cache.num_codes);
		qreg_paulirot(&c->reg, cache.code_hi, cache.codes_lo,
			cache.angles, cache.num_codes);

		cache.num_codes	  = 1;
		cache.code_hi	  = code_hi;
		cache.codes_lo[0] = code_lo;
		cache.angles[0]	  = angle;
	}

	log_trace("paulirot, last term group, num_codes: %zu", cache.num_codes);

	if (cache.num_codes > 0)
		qreg_paulirot(&c->reg, cache.code_hi, cache.codes_lo,
			cache.angles, cache.num_codes);
}

static int circ_effect(struct circ_qdrift *c)
{
	const double t = c->data->step_size;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	const double theta = asin(t);
	trott_step(c, theta);

	return 0;
}

static int circ_measure(struct circ_qdrift *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->num_dets; i++) {
		_Complex double a;
		qreg_getamp(&c->reg, md->dets[i].idx, &a);

		const _Complex double damp = md->dets[i].coeff[0] +
			_Complex_I * md->dets[i].coeff[1];
		pr += a * conj(damp);
	}
	c->prod[0] = creal(pr);
	c->prod[1] = cimag(pr);

	return 0;
}

static size_t circ_sample_invcdf(struct circ_qdrift *c, double x)
{
	size_t i   = 0;
	double cdf = 0;
	while (cdf <= x)
		cdf += fabs(c->data->hamil.coeffs[i++]);
	return i - 1; /* Never again make the same off-by-one error! */
}

static void circ_sample_terms(struct circ_qdrift *c)
{
	for (size_t i = 0; i < c->data->depth; i++) {
		double x =
			(double)(xoshiro256ss_next(&c->rng) >> 11) * 0x1.0p-53;
		c->sampled_idx[i] = circ_sample_invcdf(c, x);
	}
}

int circ_qdrift_simulate(const struct circ_qdrift_data *cd)
{
	int ret = 0;

	const size_t num_qb = cd->hamil.num_qubits;

	struct circ_qdrift c;
	if (circ_create(&c, cd, num_qb) < 0)
		goto error;

	for (size_t i = 0; i < cd->num_samples; i++) {
		circ_sample_terms(&c);
		circuit_prepst(&c);
		if (circ_effect(&c) < 0)
			goto error;
		circ_measure(&c);
		cd->samples[0][i] = c.prod[0];
		cd->samples[1][i] = c.prod[1];
	}

	goto exit;
error:
	ret = -1;
exit:
	circ_destroy(&c);

	return ret;
}
