#define LOG_SUBSYS "qdrift"

#include "c23_compat.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "container_of.h"
#include "log.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#include "circ/qdrift.h"

static uint64_t SEED = UINT64_C(0xeccd9dcc749fcdca);

static int ranct_init(struct qdrift_ranct *rct, const uint32_t qb,
	const size_t depth, const size_t cdf_len)
{
	if (prob_cdf_init(&rct->cdf, cdf_len) < 0)
		return -1;
	if (circ_hamil_init(&rct->hm_ran, qb, depth) < 0)
		return -1;

	return 0;
}

static void ranct_free(struct qdrift_ranct *rct)
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

static void ranct_calc_cdf(
	struct qdrift_ranct *rct, struct circ_hamil_term *terms)
{
	struct get_vals_data data = { .i = 0, .terms = terms };
	prob_cdf_from_iter(&rct->cdf, get_vals, &data);
}

int qdrift_init(struct qdrift *qd, const struct qdrift_data *dt,
	struct circ_hamil hm, const enum stprep_kind sp_kind,
	const void *sp_data, struct phase2_step_writer *sw)
{
	if (circ_init(&qd->ct, hm, sp_kind, sp_data, dt->samples) < 0)
		goto err_circ_init;

	qd->dt = *dt;
	qd->sw = sw;

	if (ranct_init(&qd->ranct, qd->ct.hm.qb, dt->depth, qd->ct.hm.len) < 0)
		goto err_rct_init;
	ranct_calc_cdf(&qd->ranct, qd->ct.hm.terms);

	if (qd->dt.seed != 0)
		SEED = qd->dt.seed;
	else
		qd->dt.seed = SEED;
	xoshiro256ss_init(&qd->rng, qd->dt.seed);

	return 0;

err_rct_init:
	circ_free(&qd->ct);
err_circ_init:
	return -1;
}

void qdrift_free(struct qdrift *qd)
{
	ranct_free(&qd->ranct);
	circ_free(&qd->ct);
}

static double signof(double a)
{
	const double f = fabs(a);
	if (f < DBL_EPSILON)
		return 0.0;
	return a < f ? -1.0 : 1.0;
}

static void ranct_sample(struct qdrift *qd)
{
	for (size_t i = 0; i < qd->ranct.hm_ran.len; i++) {
		const double x = xoshiro256ss_dbl01(&qd->rng);
		const size_t idx = prob_cdf_inverse(&qd->ranct.cdf, x);
		qd->ranct.hm_ran.terms[i].cf = signof(qd->ct.hm.terms[idx].cf);
		qd->ranct.hm_ran.terms[i].op = qd->ct.hm.terms[idx].op;
	}
}

int qdrift_simul(struct qdrift *qd)
{
	struct circ *ct = &qd->ct;
	struct circ_values *vals = &ct->vals;

	log_debug("simul: samples=%zu depth=%zu step_size=%g seed=%lu"
		  " cdf_len=%zu",
		vals->len, qd->dt.depth, qd->dt.step_size,
		(unsigned long)qd->dt.seed, qd->ranct.cdf.len);

	struct circ_prog prog;
	circ_prog_init(&prog, vals->len, "sample");
	for (size_t i = 0; i < vals->len; i++) {
		log_debug("sample %zu/%zu", i + 1, vals->len);
		circ_prepst(ct);

		ranct_sample(qd);
		if (circ_step(ct, &qd->ranct.hm_ran, asin(qd->dt.step_size)) <
			0) {
			log_error("simul: circ_step failed at sample %zu", i);
			return -1;
		}
		vals->z[i] = circ_measure(ct);

		if (qd->sw && qd->sw->write(qd->sw->ctx, i, vals->z[i]) < 0) {
			log_error("qdrift_simul: write_step %zu failed", i);
			return -1;
		}

		circ_prog_tick(&prog);
		circ_prog_emit(&prog, LOG_SUBSYS);
	}

	return 0;
}