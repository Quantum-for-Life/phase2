#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "circ/cmpsit.h"
#include "container_of.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#define SEED UINT64_C(0xeccd9dcc749fcdca)

int cmpsit_write_res(struct circ *c, data_id fid);
int cmpsit_simulate(struct circ *c);

static int cmpsit_res_init(struct circ_cmpsit *ct, size_t nsamples)
{
	_Complex double *samples = malloc(sizeof(_Complex double) * nsamples);
	if (!samples)
		goto err_samples;
	ct->res.samples = samples;
	ct->res.nsamples = nsamples;

	return 0;

// free(samples);
err_samples:
	return -1;
}

static void cmpsit_res_destroy(struct circ_cmpsit *ct)
{
	free(ct->res.samples);
}

int circ_cmpsit_init(
	struct circ_cmpsit *ct, struct circ_cmpsit_data *data, data_id fid)
{
	struct circ *c = &ct->circ;
	if (circ_init(c, fid) < 0)
		goto err_circ_init;
	c->simulate = cmpsit_simulate;
	c->write_res = cmpsit_write_res;

	xoshiro256ss_init(&ct->rng, SEED);

	/*
	 * Calculate the sampled circuit size (no. of pauli rotations).
	 * This is 2nd order Trotter formula. Update it, if you want to
	 * introduce higher orders.
	 */
	size_t nsmpl_ct =
		data->length + data->depth * (CIRC_CMPSIT_TRUNC_DIST + 2);
	ct->smpl_ct = malloc(sizeof *ct->smpl_ct * nsmpl_ct);
	if (!ct->smpl_ct)
		goto err_smplct;

	if (cmpsit_res_init(ct, data->samples) < 0)
		goto err_res_init;

	ct->depth = data->depth;
	ct->length = data->length;
	ct->step_size = data->step_size;
	ct->steps = data->steps;

	return 0;

err_res_init:
	free(ct->smpl_ct);
err_smplct:
	circ_destroy(c);
err_circ_init:
	return -1;
}

void circ_cmpsit_destroy(struct circ_cmpsit *ct)
{
	circ_destroy(&ct->circ);
	free(ct->smpl_ct);
	cmpsit_res_destroy(ct);
}

static int cmpsit_prepst(struct circ_cmpsit *ct)
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
	struct circ_cmpsit *ct = data;
	qreg_paulirot(&ct->circ.reg, code_hi, codes_lo, phis, ncodes);
}

static int cmpsit_step(struct circ_cmpsit *ct, const double omega)
{
	const struct circ_hamil *hamil = &ct->circ.hamil;
	struct circ_cache *cache = &ct->circ.cache;

	for (size_t i = 0; i < ct->depth; i++) {
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

static int cmpsit_effect(struct circ_cmpsit *ct)
{
	const double t = ct->step_size;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;
	const double theta = asin(t);

	return cmpsit_step(ct, theta);
}

static _Complex double cmpsit_measure(struct circ_cmpsit *ct)
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

static size_t sample_invcdf(struct circ_cmpsit *ct, double x)
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

static void sample_terms(struct circ_cmpsit *ct)
{
	for (size_t i = 0; i < ct->depth; i++) {
		double x = rand_dbl01(&ct->rng);
		/* ct->smpl[i] = sample_invcdf(ct, x); */
	}
}

int cmpsit_write_res(struct circ *c, data_id fid)
{
	int rt = -1;

	struct circ_cmpsit *ct = container_of(c, struct circ_cmpsit, circ);

	if (data_grp_create(fid, DATA_CIRCCMPSIT) < 0)
		goto data_res_write;
	if (data_attr_write(
		    fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_DEPTH, ct->depth) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_LENGTH,
		    ct->length) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_STEPSIZE,
		    ct->step_size) < 0)
		goto data_res_write;
	if (data_attr_write(
		    fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_STEPS, ct->steps) < 0)
		goto data_res_write;
	if (data_res_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_VALUES,
		    ct->res.samples, ct->res.nsamples) < 0)
		goto data_res_write;

	rt = 0;

data_res_write:
	return rt;
}

int cmpsit_simulate(struct circ *c)
{
	int rt = -1;

	size_t prog_pc = 0;
	struct circ_cmpsit *ct = container_of(c, struct circ_cmpsit, circ);

	for (size_t i = 0; i < ct->res.nsamples; i++) {
		size_t pc = i * 100 / ct->res.nsamples;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (samples: %zu)", pc, i);
		}

		sample_terms(ct);
		cmpsit_prepst(ct);
		if (cmpsit_effect(ct) < 0)
			goto ex_cmpsit_effect;
		ct->res.samples[i] = cmpsit_measure(ct);
	}

	rt = 0; /* Success. */
ex_cmpsit_effect:
	return rt;
}
