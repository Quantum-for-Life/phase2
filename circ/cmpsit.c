/*
 * Composite (partially-randomised) symmetric Trotter.
 * Init splits H into hm_det (top `length` terms by
 * |c_k|, lex-sorted, step `angle_det`) and hm_ran
 * (rest, qDRIFT-sampled at step `angle_rand`).  Each
 * Trotter step: forward half-sweep of one composite
 * sample at omega=0.5, reverse half-sweep of an
 * independent draw at omega=0.5.  Overlap goes to
 * ct.vals and, if non-NULL, cp->sw.  Error hybrid:
 * O(delta^3) deterministic + O(1/sqrt(samples))
 * stochastic.
 */

#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define LOG_SUBSYS "cmpsit"
#include "log.h"
#include "phase2.h"
#include "xoshiro256ss.h"

#include "circ/cmpsit.h"

#include "internal.h"

static uint64_t SEED = UINT64_C(0xafb424901446f21f);

#ifndef FRAC_PI_2
#define FRAC_PI_2 1.57079632679489661923132169163975144
#endif

/* qsort comparator: |cf| descending. */
static int hamil_term_cmp_abscf_desc(const void *a, const void *b)
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

static void ranct_calc_cdf(
	struct cmpsit_ranct *rct, struct circ_hamil_term *terms)
{
	/* cf is the first field of struct circ_hamil_term, so &terms[0].cf
	 * equals (double *)terms; stride is the term-record size.
	 * out_lambda receives the L1 norm used downstream. */
	prob_cdf_from_array_strided(&rct->cdf, &terms[0].cf,
		sizeof terms[0], &rct->lambda_r);
}

/* Split H by |c_k| desc: top dt->length terms into
 * hm_det (then lex-sorted), rest into hm_ran with a
 * CDF over their |c_k|. */
static int ranct_init(struct cmpsit_ranct *rct, const struct circ_hamil *hm,
	const struct cmpsit_data *dt)
{
	const uint32_t qb = hm->qb;

	struct circ_hamil hm_tmp;
	if (circ_hamil_init(&hm_tmp, qb, hm->len) < 0)
		goto err_hm_tmp;
	for (size_t i = 0; i < hm->len; i++)
		hm_tmp.terms[i] = hm->terms[i];
	qsort(hm_tmp.terms, hm_tmp.len, sizeof(struct circ_hamil_term),
		hamil_term_cmp_abscf_desc);

	if (circ_hamil_init(&rct->hm_det, qb, dt->length) < 0)
		goto err_hm_det;
	for (size_t i = 0; i < dt->length; i++)
		rct->hm_det.terms[i] = hm_tmp.terms[i];
	circ_hamil_sort_lex(&rct->hm_det);

	if (circ_hamil_init(&rct->hm_ran, qb, hm->len - dt->length) < 0)
		goto err_hm_ran;
	for (size_t i = dt->length; i < hm->len; i++)
		rct->hm_ran.terms[i - dt->length] = hm_tmp.terms[i];

	if (prob_cdf_init(&rct->cdf, hm->len - dt->length) < 0)
		goto err_cdf;

	ranct_calc_cdf(rct, hm_tmp.terms + dt->length);
	log_debug("ranct.lambda_r: %.9f", rct->lambda_r);
	rct->depth = dt->depth;
	log_debug("ranct.depth: %zu", rct->depth);
	rct->angle_rand = dt->angle_rand;
	log_debug("ranct.tau: %.9f", rct->angle_rand);

	circ_hamil_free(&hm_tmp);

	return 0;

err_cdf:
	circ_hamil_free(&rct->hm_ran);
err_hm_ran:
	circ_hamil_free(&rct->hm_det);
err_hm_det:
	circ_hamil_free(&hm_tmp);
err_hm_tmp:
	return -1;
}

static void cmpsit_ranct_free(struct cmpsit_ranct *rct)
{
	circ_hamil_free(&rct->hm_ran);
	circ_hamil_free(&rct->hm_det);
	prob_cdf_free(&rct->cdf);
}

int cmpsit_init(struct cmpsit *cp, const struct cmpsit_data *dt,
	struct circ_hamil hm, const enum stprep_kind sp_kind,
	const void *sp_data, struct phase2_step_writer *sw)
{
	if (circ_init(&cp->ct, hm, sp_kind, sp_data, dt->samples) < 0)
		goto err_circ_init;
	cp->dt = *dt;
	cp->sw = sw;
	if (ranct_init(&cp->ranct, &cp->ct.hm, dt) < 0)
		goto err_ranct_init;

	if (cp->dt.seed != 0)
		SEED = cp->dt.seed;
	else
		cp->dt.seed = SEED;
	xoshiro256ss_init(&cp->rng, cp->dt.seed);

	return 0;

err_ranct_init:
	circ_free(&cp->ct);
err_circ_init:
	return -1;
}

void cmpsit_free(struct cmpsit *cp)
{
	cmpsit_ranct_free(&cp->ranct);
	circ_free(&cp->ct);
}

static int hm_sample(struct cmpsit *cp)
{
	int rt = -1;

	const size_t len_max = cp->dt.length + cp->ranct.depth;

	if (circ_hamil_init(&cp->ranct.hm_smpl, cp->ct.hm.qb, len_max) < 0)
		goto hm_smpl_init;
	for (size_t i = 0; i < cp->dt.length; i++) {
		cp->ranct.hm_smpl.terms[i].op = cp->ranct.hm_det.terms[i].op;
		cp->ranct.hm_smpl.terms[i].cf =
			cp->ranct.hm_det.terms[i].cf * cp->dt.angle_det;
	}

	const double tau = cp->ranct.angle_rand;
	for (size_t i = cp->dt.length; i < len_max; i++) {
		const double x = xoshiro256ss_dbl01(&cp->rng);
		const size_t idx = prob_cdf_inverse(&cp->ranct.cdf, x);
		cp->ranct.hm_smpl.terms[i].op = cp->ranct.hm_ran.terms[idx].op;
		cp->ranct.hm_smpl.terms[i].cf =
			signof(cp->ranct.hm_ran.terms[idx].cf) * tau;
	}

	rt = 0;
hm_smpl_init:
	return rt;
}

static void ranct_hmsmpl_free(struct cmpsit *cp)
{
	circ_hamil_free(&cp->ranct.hm_smpl);
}

/* Draw a fresh composite sample, apply forward or
 * reverse half-step at omega=0.5, free the sample. */
static int half_step(struct cmpsit *cp, size_t s, bool reverse)
{
	const char *tag = reverse ? "rev" : "fwd";
	log_trace("step %zu/%zu %s half-sweep", s + 1, cp->dt.steps, tag);
	(void)s; /* unused when log_trace expands to no-op */

	if (hm_sample(cp) < 0) {
		log_error("simul: hm_sample failed");
		return -1;
	}
	const int rc = reverse
		? circ_step_reverse(&cp->ct, &cp->ranct.hm_smpl, 0.5)
		: circ_step(&cp->ct, &cp->ranct.hm_smpl, 0.5);
	if (rc < 0) {
		log_error("simul: %s circ_step failed", tag);
		ranct_hmsmpl_free(cp);
		return -1;
	}
	ranct_hmsmpl_free(cp);
	return 0;
}

int cmpsit_simul(struct cmpsit *cp)
{
	struct circ *ct = &cp->ct;
	struct circ_values *vals = &ct->vals;

	log_debug("simul: samples=%zu steps=%zu length=%zu depth=%zu"
		  " angle_det=%g angle_rand=%g",
		vals->len, cp->dt.steps, cp->dt.length, cp->dt.depth,
		cp->dt.angle_det, cp->dt.angle_rand);

	for (size_t i = 0; i < vals->len; i++) {
		log_debug("sample %zu/%zu", i + 1, vals->len);
		circ_prepst(ct);
		for (size_t s = 0; s < cp->dt.steps; s++) {
			if (half_step(cp, s, false) < 0)
				return -1;
			if (half_step(cp, s, true) < 0)
				return -1;
		}
		vals->z[i] = circ_measure(ct);
		if (cp->sw && cp->sw->write(cp->sw->ctx, i, vals->z[i]) < 0) {
			log_error("cmpsit_simul: write_step %zu failed", i);
			return -1;
		}
	}

	return 0;
}
