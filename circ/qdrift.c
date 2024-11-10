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

static int qdrift_ranct_init(struct qdrift_ranct *rct, const uint32_t qb,
	const size_t depth, const size_t cdf_len)
{
	if (prob_cdf_init(&rct->cdf, cdf_len) < 0)
		return -1;
	if (circ_hamil_init(&rct->hm_ran, qb, depth) < 0)
		return -1;

	return 0;
}

static void qdrift_ranct_free(struct qdrift_ranct *rct)
{
	circ_hamil_free(&rct->hm_ran);
	prob_cdf_free(&rct->cdf);
}

struct get_vals_data {
	size_t i;
	struct circ_hamil_term *terms;
};

static double get_vals(void *data)
{
	struct get_vals_data *dt = data;

	return dt->terms[dt->i++].cf;
}

static void qdrift_ranct_calc_cdf(
	struct qdrift_ranct *rct, struct circ_hamil_term *terms)
{
	struct get_vals_data data = { .i = 0, .terms = terms };
	prob_cdf_from(&rct->cdf, get_vals, &data);
}

static int qdrift_smpl_init(struct qdrift_smpl *smp, size_t samples)
{
	_Complex double *z = malloc(sizeof(_Complex double) * samples);
	if (!z)
		return -1;
	smp->z = z;
	smp->len = samples;

	return 0;
}

static void qdrift_smpl_free(struct qdrift_smpl *smp)
{
	free(smp->z);
}

int qdrift_init(
	struct qdrift *qd, const struct qdrift_data *dt, const data_id fid)
{
	if (circ_init(&qd->ct, fid) < 0)
		goto err_circ_init;

	qd->dt = *dt;

	if (qdrift_ranct_init(
		    &qd->ranct, qd->ct.hm.qb, dt->depth, qd->ct.hm.len) < 0)
		goto err_rct_init;
	qdrift_ranct_calc_cdf(&qd->ranct, qd->ct.hm.terms);

	if (qdrift_smpl_init(&qd->smpl, dt->samples) < 0)
		goto err_samples_init;

	xoshiro256ss_init(&qd->rng, SEED);

	return 0;

	// qdrift_samples_destroy(&qd->smpl);
err_samples_init:
	qdrift_ranct_free(&qd->ranct);
err_rct_init:
	circ_free(&qd->ct);
err_circ_init:
	return -1;
}

void qdrift_free(struct qdrift *qd)
{
	qdrift_smpl_free(&qd->smpl);
	qdrift_ranct_free(&qd->ranct);
	circ_free(&qd->ct);
}

static void qdrift_ranct_sample(struct qdrift *qd)
{
	for (size_t i = 0; i < qd->ranct.hm_ran.len; i++) {
		const double x = xoshiro256ss_dbl01(&qd->rng);
		const size_t idx = prob_cdf_inverse(&qd->ranct.cdf, x);
		qd->ranct.hm_ran.terms[i].cf = 1.0; // qd->ct.hm.terms[idx].cf;
		qd->ranct.hm_ran.terms[i].op = qd->ct.hm.terms[idx].op;
	}
}

int qdrift_simul(struct qdrift *qd)
{
	struct circ *ct = &qd->ct;
	struct circ_prog prog;

	circ_prog_init(&prog, qd->smpl.len);
	for (size_t i = 0; i < qd->smpl.len; i++) {
		qdrift_ranct_sample(qd);

		circ_prepst(ct);
		if (circ_step(ct, &qd->ranct.hm_ran, asin(qd->dt.step_size)) < 0)
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