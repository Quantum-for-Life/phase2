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

static int qdrift_rct_init(
	struct qdrift_rct *rct, const uint32_t qb, const size_t depth)
{
	return circ_hamil_init(&rct->rhm, qb, depth);
}

static void qdrift_rct_destroy(struct qdrift_rct *rct)
{
	circ_hamil_free(&rct->rhm);
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

static void qdrift_samples_destroy(struct qdrift_samples *smp)
{
	free(smp->z);
}

int qdrift_init(struct qdrift *qd, const struct qdrift_data *dt, data_id fid)
{
	struct circ *c = &qd->ct;
	if (circ_init(c, fid, qdrift_simulate) < 0)
		goto err_circ_init;

	qd->dt = *dt;

	if (qdrift_rct_init(&qd->rct, qd->ct.hamil.qb, dt->depth) < 0)
		goto err_rct_init;
	if (qdrift_samples_init(&qd->smp, dt->samples) < 0)
		goto err_samples_init;

	xoshiro256ss_init(&qd->rng, SEED);

	return 0;

	// qdrift_samples_destroy(&qd->smp);
err_samples_init:
	qdrift_rct_destroy(&qd->rct);
err_rct_init:
	circ_free(c);
err_circ_init:
	return -1;
}

void qdrift_free(struct qdrift *qd)
{
	qdrift_samples_destroy(&qd->smp);
	qdrift_rct_destroy(&qd->rct);
	circ_free(&qd->ct);
}

static size_t sample_invcdf(struct qdrift *qd, double x)
{
	(void)qd;
	size_t i = 0;
	double cdf = 0;
	while (cdf <= x)
		cdf += fabs(qd->ct.hamil.terms[i++].cf);
	return i - 1; /* Never again make the same off-by-one error! */
}

static void sample_terms(struct qdrift *qd)
{
	for (size_t i = 0; i < qd->rct.rhm.len; i++) {
		const double x = xoshiro256ss_dbl01(&qd->rng);
		const size_t idx = sample_invcdf(qd, x);
		qd->rct.rhm.terms[i] = qd->ct.hamil.terms[idx];
	}
}

int qdrift_simulate(struct circ *ct)
{
	int rt = -1;

	size_t prog_pc = 0;
	struct qdrift *qd = container_of(ct, struct qdrift, ct);

	for (size_t i = 0; i < qd->smp.len; i++) {
		size_t pc = i * 100 / qd->smp.len;
		if (pc > prog_pc) {
			prog_pc = pc;
			log_info("Progress: %zu\% (samples: %zu)", pc, i);
		}

		circ_prepst(ct);
		sample_terms(qd);
		if (circ_step(ct, &qd->rct.rhm, asin(qd->dt.step_size)) < 0)
			goto ex_qdrift_effect;
		qd->smp.z[i] = circ_measure(ct);
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
		    qd->smp.z, qd->smp.len) < 0)
		goto data_res_write;

	rt = 0;

data_res_write:
	return rt;
}