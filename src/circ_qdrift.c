#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "circ_qdrift.h"
#include "log.h"
#include "qreg.h"
#include "xoshiro256starstar.h"

#define MAX_CACHE_CODES (1024)

#define PRNG_SEED (0x235eac32)

struct circ {
	size_t	    num_qb;
	struct qreg reg;

	const struct circ_data *data;

	double prod[2];

	struct code_cache {
		struct paulis code_hi;
		struct paulis codes_lo[MAX_CACHE_CODES];
		fl	      angles[MAX_CACHE_CODES];
		size_t	      num_codes;
	} cache;

	struct xoshiro256starstar rng;
	size_t			 *sampled_idx;
};

static int
circ_create(
	struct circ *c, const struct circ_data *data, const size_t num_qubits)
{
	struct qreg reg;
	if (qreg_init(&reg, num_qubits) < 0)
		return -1;

	c->num_qb = num_qubits;
	c->data	  = data;
	c->reg	  = reg;

	xoshiro256starstar_init(&c->rng, PRNG_SEED);
	size_t *sampled_idx = malloc(sizeof(size_t) * data->depth);
	if (sampled_idx == NULL)
		return -1;
	c->sampled_idx = sampled_idx;

	return 0;
}

static void
circ_destroy(struct circ *c)
{
	qreg_destroy(&c->reg);
	if (c->sampled_idx)
		free(c->sampled_idx);
}

static void
circ_hamil_init(struct circ_hamil *h)
{
	h->num_qubits = 0;
	h->num_terms  = 0;
	h->coeffs     = NULL;
	h->paulis     = NULL;
}

static void
circ_hamil_destroy(struct circ_hamil *h)
{
	if (h->paulis) {
		free(h->paulis);
		h->paulis = NULL;
	}
	if (h->coeffs) {
		free(h->coeffs);
		h->coeffs = NULL;
	}
	h->num_terms  = 0;
	h->num_qubits = 0;
}

struct hamil_iter_data {
	size_t	       idx;
	size_t	       num_qubits;
	double	       norm;
	double	      *coeffs;
	struct paulis *paulis;
};

static int
hamil_iter(double coeff, unsigned char *paulis, void *iter_data)
{
	struct hamil_iter_data *idat	   = iter_data;
	const size_t		i	   = idat->idx++;
	const size_t		num_qubits = idat->num_qubits;

	idat->coeffs[i] = coeff * idat->norm;
	struct paulis p = paulis_new();
	for (size_t j = 0; j < num_qubits; j++) {
		paulis_set(&p, paulis[j], j);
	}
	idat->paulis[i] = p;

	return 0;
}

static int
circ_hamil_from_file(struct circ_hamil *h, const data_id fid)
{
	size_t num_qubits, num_terms;
	double norm;

	if (data_hamil_getnums(fid, &num_qubits, &num_terms) < 0)
		return -1;
	if (data_hamil_getnorm(fid, &norm) < 0)
		return -1;

	double	      *coeffs = malloc(sizeof *coeffs * num_terms);
	struct paulis *paulis = malloc(sizeof(struct paulis) * num_terms);
	if (!(coeffs && paulis))
		goto err;

	struct hamil_iter_data idat = { .idx = 0,
		.num_qubits		     = num_qubits,
		.norm			     = norm,
		.coeffs			     = coeffs,
		.paulis			     = paulis };
	if (data_hamil_foreach(fid, hamil_iter, &idat) != 0)
		goto err;

	h->num_qubits = num_qubits;
	h->num_terms  = num_terms;
	h->coeffs     = coeffs;
	h->paulis     = paulis;

	return 0;
err:
	free(coeffs);
	free(paulis);
	return -1;
}

static void
circ_multidet_init(struct circ_multidet *md)
{
	md->num_dets = 0;
	md->dets     = NULL;
}

static void
circ_multidet_destroy(struct circ_multidet *md)
{
	if (md->dets) {
		free(md->dets);
		md->dets = NULL;
	}
	md->num_dets = 0;
}

struct iter_multidet_data {
	size_t		      i;
	struct circ_multidet *md;
};

static int
iter_multidet(double coeff[2], const uint64_t idx, void *op_data)
{
	struct iter_multidet_data *imd = op_data;

	imd->md->dets[imd->i].coeff[0] = coeff[0];
	imd->md->dets[imd->i].coeff[1] = coeff[1];
	imd->md->dets[imd->i].idx      = idx;
	imd->i++;

	return 0;
}

static int
circuit_multidet_from_data(struct circ_multidet *md, const data_id fid)
{
	size_t num_qubits, num_dets;
	if (data_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		return -1;
	md->dets = malloc(sizeof *md->dets * num_dets);
	if (md->dets == NULL)
		return -1;

	struct iter_multidet_data imd;
	imd.i  = 0;
	imd.md = md;
	if (data_multidet_foreach(fid, iter_multidet, &imd) < 0)
		goto error;

	md->num_dets = num_dets;

	return 0;
error:
	free(md->dets);
	return -1;
}

int
circ_data_from_file(struct circ_data *cd, const data_id fid)
{
	int rc = circ_hamil_from_file(&cd->hamil, fid);
	rc |= circuit_multidet_from_data(&cd->multidet, fid);
	data_trotter_get_factor(fid, &cd->step_size);
	data_trotter_get_num_samples(fid, &cd->num_samples);
	data_trotter_get_depth(fid, &cd->depth);

	return rc;
}

int
circ_data_init(struct circ_data *cd, const data_id fid)
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

void
circ_data_destroy(struct circ_data *cd)
{
	circ_multidet_destroy(&cd->multidet);
	circ_hamil_destroy(&cd->hamil);

	free(cd->samples[0]);
}

static int
circuit_prepst(struct circ *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	qreg_zero(&c->reg);
	for (size_t i = 0; i < md->num_dets; i++) {
		const fl coeff[2] = {
			/* We cast each element to a (possibly) lower precision
			 * floating point number: fl */
			md->dets[i].coeff[0], md->dets[i].coeff[1]
		};
		qreg_setamp(&c->reg, md->dets[i].idx, coeff);
	}

	return 0;
}

static void
trotter_step(struct circ *c, const double omega)
{
	const struct circ_hamil *hamil = &c->data->hamil;

	struct code_cache cache = c->cache;
	cache.num_codes		= 0;

	for (size_t i = 0; i < c->data->depth; i++) {
		const fl	    angle     = omega;
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

static int
circ_effect(struct circ *c)
{
	const double t = c->data->step_size;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	const double theta = asin(t);
	trotter_step(c, theta);

	return 0;
}

static int
circ_measure(struct circ *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	double pr[2] = { 0.0, 0.0 };
	for (size_t i = 0; i < md->num_dets; i++) {
		fl amp[2];
		qreg_getamp(&c->reg, md->dets[i].idx, &amp);

		const double damp_re = md->dets[i].coeff[0];
		const double damp_im = md->dets[i].coeff[1];
		/* inner product with damp complex-conjugated */
		pr[0] += damp_re * amp[0] + damp_im * amp[1];
		pr[1] += damp_re * amp[1] - damp_im * amp[0];
	}
	c->prod[0] = pr[0];
	c->prod[1] = pr[1];

	return 0;
}

static size_t
circ_sample_invcdf(struct circ *c, double x)
{
	size_t i   = 0;
	double cdf = 0;
	while (cdf <= x)
		cdf += fabs(c->data->hamil.coeffs[i++]);
	return i - 1; /* Never again make the same off-by-one error! */
}

static void
circ_sample_terms(struct circ *c)
{
	for (size_t i = 0; i < c->data->depth; i++) {
		double x = (double)(xoshiro256starstar_next(&c->rng) >> 11) *
			   0x1.0p-53;
		c->sampled_idx[i] = circ_sample_invcdf(c, x);
	}
}

int
circ_simulate(const struct circ_data *cd)
{
	int ret = 0;

	const size_t num_qb = cd->hamil.num_qubits;

	struct circ c;
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
