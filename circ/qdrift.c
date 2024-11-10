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

int qdrift_simulate(struct circ *ct);

static int qdrift_rct_init(struct qdrift_rct *rct, const uint32_t qb,
	const size_t depth, const size_t pd_len)
{
	if (prob_pd_init(&rct->pd, pd_len) < 0)
		return -1;
	if (circ_hamil_init(&rct->hm, qb, depth) < 0)
		return -1;

	return 0;
}

static void qdrift_rct_free(struct qdrift_rct *rct)
{
	circ_hamil_free(&rct->hm);
	prob_pd_free(&rct->pd);
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
	struct qdrift_rct *rct, struct circ_hamil_term *terms)
{
	struct get_smpl_data data = { .i = 0, .terms = terms };
	prob_pdf_from_samples(&rct->pd, get_smpl, &data);
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
	if (circ_init(&qd->ct, fid, qdrift_simulate) < 0)
		goto err_circ_init;

	qd->dt = *dt;

	if (qdrift_rct_init(
		    &qd->rct, qd->ct.hamil.qb, dt->depth, qd->ct.hamil.len) < 0)
		goto err_rct_init;
	qdrift_rct_calc_pd(&qd->rct, qd->ct.hamil.terms);

	if (qdrift_samples_init(&qd->smpl, dt->samples) < 0)
		goto err_samples_init;

	xoshiro256ss_init(&qd->rng, SEED);

	return 0;

	// qdrift_samples_destroy(&qd->smpl);
err_samples_init:
	qdrift_rct_free(&qd->rct);
err_rct_init:
	circ_free(&qd->ct);
err_circ_init:
	return -1;
}

void qdrift_free(struct qdrift *qd)
{
	qdrift_samples_free(&qd->smpl);
	qdrift_rct_free(&qd->rct);
	circ_free(&qd->ct);
}

static void qdrift_rct_sample_terms(struct qdrift *qd)
{
	for (size_t i = 0; i < qd->rct.hm.len; i++) {
		const double x = xoshiro256ss_dbl01(&qd->rng);
		const size_t idx = prob_cdf_inverse(&qd->rct.pd, x);
		qd->rct.hm.terms[i] = qd->ct.hamil.terms[idx];
	}
}

int qdrift_simulate(struct circ *ct)
{
	int rt = -1;

	size_t prog_pc = 0;
	struct qdrift *qd = container_of(ct, struct qdrift, ct);

	for (size_t i = 0; i < qd->smpl.len; i++) {
		size_t pc = i * 100 / qd->smpl.len;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (samples: %zu)", pc, i);
		}

		circ_prepst(ct);
		qdrift_rct_sample_terms(qd);
		if (circ_step(ct, &qd->rct.hm, asin(qd->dt.step_size)) < 0)
			goto ex_qdrift_effect;
		qd->smpl.z[i] = circ_measure(ct);
	}

	rt = 0; /* Success. */
ex_qdrift_effect:
	return rt;
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