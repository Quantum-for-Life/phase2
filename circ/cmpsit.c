#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "container_of.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#include "circ/cmpsit.h"

#define SEED UINT64_C(0xafb424901446f21f)

static int cmpsit_write_res(struct circ *c, data_id fid);
static int cmpsit_simulate(struct circ *c);

/* Sort in _descending_ order */
static int cmp_hamil_term_by_cf(const void *a, const void *b)
{
	const struct circ_hamil_term ta = *(const struct circ_hamil_term *)a;
	const struct circ_hamil_term tb = *(const struct circ_hamil_term *)b;

	const double x = fabs(ta.cf);
	const double y = fabs(tb.cf);

	if (x > y)
		return -1;
	if (x < y)
		return 1;
	return 0;
}

static int cmp_hamil_term_by_lex(const void *a, const void *b)
{
	const struct paulis x = ((const struct circ_hamil_term *)a)->op;
	const struct paulis y = ((const struct circ_hamil_term *)b)->op;

	return paulis_cmp(x, y);
}

/*
 * Sort the Hamiltonian by absolute val of coefficients, in descending order.
 *
 * The first data->length terms are "deterministic" and they can be sorted
 * back by lexicographical order of Pauli strings to improve performance of
 * the distributed computation.
 *
 * The remaining terms will be randomly sampled from.
 */
static void cmpsit_hamil_rearrange(
	const struct cmpsit *ct, const struct cmpsit_data *data)
{
	qsort(ct->circ.hamil.terms, ct->circ.hamil.nterms,
		sizeof(struct circ_hamil_term), cmp_hamil_term_by_cf);
	qsort(ct->circ.hamil.terms, data->length,
		sizeof(struct circ_hamil_term), cmp_hamil_term_by_lex);
}

static int cmpsit_pd_init(struct cmpsit_pd *pd,
	const struct circ_hamil_term *src, const size_t len)
{
	double *const x = malloc(sizeof(double) * len);
	if (!x)
		return -1;

	double lambda_r = 0.0;
	for (size_t i = 0; i < len; i++)
		lambda_r += fabs(src[i].cf);
	for (size_t i = 0; i < len; i++)
		x[i] = fabs(src[i].cf) / lambda_r;

	pd->x = x;
	pd->len = len;
	pd->lambda_r = lambda_r;

	return 0;
}

static void cmpsit_pd_destroy(struct cmpsit_pd *pd)
{
	free(pd->x);
}

static int cmpsit_rct_init(struct cmpsit_rct *rct, const size_t len)
{
	struct circ_hamil_term *trm =
		malloc(sizeof(struct circ_hamil_term) * len);
	if (!trm)
		return -1;

	rct->trm = trm;
	rct->len = len;

	return 0;
}

static void cmpsit_rct_destroy(struct cmpsit_rct *rct)
{
	free(rct->trm);
}

static int cmpsit_samples_init(struct cmpsit_samples *samples, size_t len)
{
	_Complex double *a = malloc(sizeof(_Complex double) * len);
	if (!a)
		return -1;
	samples->a = a;
	samples->len = len;

	return 0;
}

static void cmpsit_samples_destroy(struct cmpsit_samples *samples)
{
	free(samples->a);
}

int cmpsit_init(struct cmpsit *ct,
	const struct cmpsit_data *data, const data_id fid)
{
	struct circ *c = &ct->circ;
	if (circ_init(c, fid) < 0)
		goto err_circ_init;
	c->simulate = cmpsit_simulate;
	c->write_res = cmpsit_write_res;

	ct->data = *data;

	cmpsit_hamil_rearrange(ct, data);

	const struct circ_hamil_term *src = ct->circ.hamil.terms + data->length;
	if (cmpsit_pd_init(&ct->pd, src, data->length) < 0)
		goto err_pd_init;

	/*
	 * Calculate the sampled circuit size (no. of pauli rotations).
	 * This is 2nd order Trotter formula. Update it, if you want to
	 * introduce higher orders.
	 */
	const size_t rct_len =
		data->length + data->depth * (CMPSIT_TRUNC_DIST + 2);
	if (cmpsit_rct_init(&ct->rct, rct_len) < 0)
		goto err_rct_init;

	if (cmpsit_samples_init(&ct->samples, data->samples) < 0)
		goto err_samples_init;

	/* TODO: seed it with user-supplied seed */
	xoshiro256ss_init(&ct->rng, SEED);

	return 0;

	// cmpsit_samples_destroy(&ct->samples);
err_samples_init:
	cmpsit_rct_destroy(&ct->rct);
err_rct_init:
	cmpsit_pd_destroy(&ct->pd);
err_pd_init:
	circ_destroy(c);
err_circ_init:
	return -1;
}

void cmpsit_destroy(struct cmpsit *ct)
{
	cmpsit_samples_destroy(&ct->samples);
	cmpsit_rct_destroy(&ct->rct);
	cmpsit_pd_destroy(&ct->pd);
	circ_destroy(&ct->circ);
}

static int cmpsit_prepst(struct cmpsit *ct)
{
	const struct circ_muldet *md = &ct->circ.muldet;

	qreg_zero(&ct->circ.reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&ct->circ.reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void cmpsit_flush(struct paulis code_hi, struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct cmpsit *ct = data;
	qreg_paulirot(&ct->circ.reg, code_hi, codes_lo, phis, ncodes);
}

static int cmpsit_step(struct cmpsit *ct, const double omega)
{
	const struct circ_hamil *hamil = &ct->circ.hamil;
	struct circ_cache *cache = &ct->circ.cache;

	for (size_t i = 0; i < ct->data.depth; i++) {
		const double phi = omega;
		/* const size_t i_smpl = ct->smpl[i];
		const struct paulis code = hamil->terms[i_smpl].op;

		if (circ_cache_insert(cache, code, phi) == 0)
			continue;

		log_trace("paulirot, term: %zu, num_codes: %zu", i, cache->n);
		circ_cache_flush(cache, cmpsit_flush, ct);
		if (circ_cache_insert(cache, code, phi) < 0)
			return -1;
		*/
	}
	log_trace("paulirot, last term group, num_codes: %zu", cache->n);
	circ_cache_flush(cache, cmpsit_flush, ct);

	return 0;
}

static int cmpsit_effect(struct cmpsit *ct)
{
	const double t = ct->data.step_size;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;
	const double theta = asin(t);

	return cmpsit_step(ct, theta);
}

static _Complex double cmpsit_measure(struct cmpsit *ct)
{
	const struct circ_muldet *md = &ct->circ.muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&ct->circ.reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}

	return pr;
}

static size_t sample_invcdf(struct cmpsit *ct, double x)
{
	(void)ct;
	size_t i = 0;
	double cdf = 0;
	while (cdf <= x)
		cdf += fabs(ct->circ.hamil.terms[i++].cf);
	return i - 1; /* Never again make the same off-by-one error! */
}

/* TODO: Move to xoshiro256ss.h and document. */
#define rand_dbl01(rng) ((double)(xoshiro256ss_next(rng) >> 11) * 0x1.0p-53)

static void sample_terms(struct cmpsit *ct)
{
	for (size_t i = 0; i < ct->data.depth; i++) {
		double x = rand_dbl01(&ct->rng);
		/* ct->smpl[i] = sample_invcdf(ct, x); */
	}
}

static int cmpsit_simulate(struct circ *c)
{
	int rt = -1;

	size_t prog_pc = 0;
	struct cmpsit *ct = container_of(c, struct cmpsit, circ);

	for (size_t i = 0; i < ct->samples.len; i++) {
		size_t pc = i * 100 / ct->samples.len;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (samples: %zu)", pc, i);
		}

		sample_terms(ct);
		cmpsit_prepst(ct);
		if (cmpsit_effect(ct) < 0)
			goto ex_cmpsit_effect;
		ct->samples.a[i] = cmpsit_measure(ct);
	}

	rt = 0; /* Success. */
ex_cmpsit_effect:
	return rt;
}

static int cmpsit_write_res(struct circ *c, data_id fid)
{
	int rt = -1;

	struct cmpsit *ct = container_of(c, struct cmpsit, circ);

	if (data_grp_create(fid, DATA_CIRCCMPSIT) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_DEPTH,
		    ct->data.depth) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_LENGTH,
		    ct->data.length) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_STEPSIZE,
		    ct->data.step_size) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_STEPS,
		    ct->data.steps) < 0)
		goto data_res_write;
	if (data_res_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_VALUES,
		    ct->samples.a, ct->samples.len) < 0)
		goto data_res_write;

	rt = 0;

data_res_write:
	return rt;
}
