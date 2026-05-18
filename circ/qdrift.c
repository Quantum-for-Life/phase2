/*
 * circ/qdrift.c -- qDRIFT randomised product formula
 * (Campbell, 2019).  See include/circ/qdrift.h for the
 * public API.
 *
 * The driver builds a normalised CDF from the
 * Hamiltonian's |c_k| weights and, for each of
 * `dt.samples` independent runs, draws `dt.depth`
 * terms i.i.d. from that distribution and applies one
 * rotation per draw at angle asin(step_size) carrying
 * sign(c_k).  Per-sample overlap goes into `ct.vals`
 * and through `qd->sw`.
 *
 * Module-private helpers live in the unnamed block
 * below:
 *   - ranct_init / _free      -- carrier for the
 *                                random pool + CDF.
 *   - ranct_calc_cdf          -- thin wrapper over
 *                                prob_cdf_from_array_strided.
 *   - ranct_sample            -- one independent draw,
 *                                filling hm_ran with
 *                                signed unit-weight
 *                                terms.
 *
 * The PRNG is xoshiro256** seeded from `dt.seed` (or a
 * static fallback when zero).  Sign carrying uses
 * `signof` from circ/internal.h.  Error order is
 * stochastic O(1/sqrt(samples)); see doc/phase2.md
 * §5.3.
 */

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

#include "internal.h"

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

static void ranct_calc_cdf(
	struct qdrift_ranct *rct, struct circ_hamil_term *terms)
{
	/* cf is the first field of struct circ_hamil_term, so &terms[0].cf
	 * equals (double *)terms; stride is the term-record size. */
	prob_cdf_from_array_strided(&rct->cdf, &terms[0].cf,
		sizeof terms[0], NULL);
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
	}

	return 0;
}