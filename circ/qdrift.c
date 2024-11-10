#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "container_of.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#include "circ/qdrift.h"

#define SEED UINT64_C(0xeccd9dcc749fcdca)

static int qdrift_simul(struct circ *ct);

static int qdrift_rct_init(struct qdrift_randct *rct, const uint32_t qb,
	const size_t depth, const size_t cdf_len)
{
	if (prob_cdf_init(&rct->cdf, cdf_len) < 0)
		return -1;
	if (circ_hamil_init(&rct->hm, qb, depth) < 0)
		return -1;

	return 0;
}

static void qdrift_rct_free(struct qdrift_randct *rct)
{
	circ_hamil_free(&rct->hm);
	prob_cdf_free(&rct->cdf);
}

struct get_smpl_data {
	size_t i;
	struct circ_hamil_term *terms;
};

static double get_smpl(void *data)
{
	struct get_smpl_data *dt = data;

	return dt->terms[dt->i++].cf;
}

static void qdrift_rct_calc_pd(
	struct qdrift_randct *rct, struct circ_hamil_term *terms)
{
	struct get_smpl_data data = { .i = 0, .terms = terms };
	prob_cdf_from_samples(&rct->cdf, get_smpl, &data);
}

static int qdrift_samples_init(struct qdrift_samples *smp, size_t samples)
{
	_Complex double *z = malloc(sizeof(_Complex double) * samples);
	if (!z)
		return -1;
	smp->z = z;
	smp->len = samples;

	return 0;
}

static void qdrift_samples_free(struct qdrift_samples *smp)
{
	free(smp->z);
}

int qdrift_init(
	struct qdrift *qd, const struct qdrift_data *dt, const data_id fid)
{
	if (circ_init(&qd->ct, fid, qdrift_simul) < 0)
		goto err_circ_init;

	qd->dt = *dt;

	if (qdrift_rct_init(
		    &qd->randct, qd->ct.hm.qb, dt->depth, qd->ct.hm.len) < 0)
		goto err_rct_init;
	qdrift_rct_calc_pd(&qd->randct, qd->ct.hm.terms);

	if (qdrift_samples_init(&qd->smpl, dt->samples) < 0)
		goto err_samples_init;

	xoshiro256ss_init(&qd->rng, SEED);

	return 0;

	// qdrift_samples_destroy(&qd->smpl);
err_samples_init:
	qdrift_rct_free(&qd->randct);
err_rct_init:
	circ_free(&qd->ct);
err_circ_init:
	return -1;
}

void qdrift_free(struct qdrift *qd)
{
	qdrift_samples_free(&qd->smpl);
	qdrift_rct_free(&qd->randct);
	circ_free(&qd->ct);
}

static void qdrift_rct_sample_terms(struct qdrift *qd)
{
	for (size_t i = 0; i < qd->randct.hm.len; i++) {
		const double x = xoshiro256ss_dbl01(&qd->rng);
		const size_t idx = prob_cdf_inverse(&qd->randct.cdf, x);
		qd->randct.hm.terms[i].cf = 1.0; // qd->ct.hm.terms[idx].cf;
		qd->randct.hm.terms[i].op = qd->ct.hm.terms[idx].op;
	}
}

static int qdrift_simul(struct circ *ct)
{
	struct qdrift *qd = container_of(ct, struct qdrift, ct);
	struct circ_prog prog;

	circ_prog_init(&prog, qd->smpl.len);
	for (size_t i = 0; i < qd->smpl.len; i++) {
		qdrift_rct_sample_terms(qd);

		circ_prepst(ct);
		if (circ_step(ct, &qd->randct.hm, asin(qd->dt.step_size)) < 0)
			return -1;
		qd->smpl.z[i] = circ_measure(ct);

		circ_prog_tick(&prog);
	}

	return 0;
}

int qdrift_write_res(struct qdrift *qd, data_id fid)
{
	int rt = -1;

	if (data_grp_create(fid, DATA_CIRCQDRIFT) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_STEPSIZE,
		    qd->dt.step_size) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_DEPTH,
		    qd->dt.depth) < 0)
		goto data_res_write;
	if (data_res_write(fid, DATA_CIRCQDRIFT, DATA_CIRCQDRIFT_VALUES,
		    qd->smpl.z, qd->smpl.len) < 0)
		goto data_res_write;

	rt = 0;

data_res_write:
	return rt;
}