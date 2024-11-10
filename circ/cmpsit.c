#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2.h"
#include "xoshiro256ss.h"

#include "circ/cmpsit.h"

#define SEED UINT64_C(0xafb424901446f21f)

/*
 * Sort the Hamiltonian by absolute val of coefficients, in descending order.
 */
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
	struct cmpsit_ranct *rct, struct circ_hamil_term *terms)
{
	struct get_vals_data data = { .i = 0, .terms = terms };
	prob_cdf_from_iter(&rct->cdf, get_vals, &data);
}

static int ranct_init(struct cmpsit_ranct *rct, const struct circ_hamil *hm,
	const size_t length)
{
	if (prob_cdf_init(&rct->cdf, hm->len - length) < 0)
		return -1;

	const uint32_t qb = hm->qb;

	struct circ_hamil hm_tmp;
	if (circ_hamil_init(&hm_tmp, qb, hm->len) < 0)
		goto err_hm_tmp;
	for (size_t i = 0; i < hm->len; i++)
		hm_tmp.terms[i] = hm->terms[i];
	qsort(hm_tmp.terms, hm_tmp.len, sizeof(struct circ_hamil_term),
		hamil_term_cmp_abscf_desc);
	ranct_calc_cdf(rct, hm_tmp.terms + length);

	if (circ_hamil_init(&rct->hm_det, qb, length) < 0)
		goto err_hm_det;
	for (size_t i = 0; i < length; i++)
		rct->hm_det.terms[i] = hm_tmp.terms[i];
	circ_hamil_sort_lex(&rct->hm_det);

	if (circ_hamil_init(&rct->hm_ran, qb, hm->len - length) < 0)
		goto err_hm_ran;
	for (size_t i = length; i < hm->len; i++)
		rct->hm_ran.terms[i - length] = hm_tmp.terms[i];
	circ_hamil_sort_lex(&rct->hm_det);

	circ_hamil_free(&hm_tmp);

	return 0;

	// circ_hamil_free(&rct->hm_ran);
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

int cmpsit_init(
	struct cmpsit *cp, const struct cmpsit_data *dt, const data_id fid)
{
	if (circ_init(&cp->ct, fid, dt->samples) < 0)
		goto err_circ_init;

	cp->dt = *dt;

	if (ranct_init(&cp->ranct, &cp->ct.hm, dt->length) < 0)
		goto err_ranct_init;

	/* TODO: seed it with user-supplied seed */
	xoshiro256ss_init(&cp->rng, SEED);

	return 0;

	// ranct_free(&cp->ranct);
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

static int ranct_hmsmpl(struct cmpsit *cp)
{
	int rt = -1;

	const size_t hm_smpl_len_max =
		cp->dt.length + cp->dt.depth * CMPSIT_TRUNC_DIST;
	struct circ_hamil_term *trm = malloc(sizeof *trm * hm_smpl_len_max);
	if (!trm)
		goto trm_alloc;

	size_t i;
	for (i = 0; i < cp->dt.length; i++)
		trm[i] = cp->ranct.hm_det.terms[i];
	for (size_t d = 0; d < cp->dt.depth; d++) {
		if (i == hm_smpl_len_max) {
			log_warn("Max. sampled circuit size reached: %zu", i);
			break;
		}

		/* Sample terms */
	}

	const size_t hm_smpl_len = i;
	if (circ_hamil_init(&cp->ranct.hm_smpl, cp->ct.hm.qb, hm_smpl_len) < 0)
		goto hm_smpl_init;
	for (i = 0; i < hm_smpl_len; i++)
		cp->ranct.hm_smpl.terms[i] = trm[i];

	rt = 0;
hm_smpl_init:
	free(trm);
trm_alloc:
	return rt;
}

static void ranct_hmsmpl_free(struct cmpsit *cp)
{
	circ_hamil_free(&cp->ranct.hm_smpl);
}

int cmpsit_simul(struct cmpsit *cp)
{
	/* Second order Trotter */
	struct circ *ct = &cp->ct;
	struct circ_values *vals = &ct->vals;

	struct circ_prog prog;
	circ_prog_init(&prog, vals->len);
	for (size_t i = 0; i < vals->len; i++) {
		if (ranct_hmsmpl(cp) < 0)
			return -1;

		circ_prepst(ct);
		if (circ_step(&cp->ct, &cp->ranct.hm_smpl, 1.0) < 0)
			return -1;
		vals->z[i] = circ_measure(ct);

		ranct_hmsmpl_free(cp);
		circ_prog_tick(&prog);
	}

	return 0;
}

int cmpsit_write_res(struct cmpsit *cp, data_id fid)
{
	int rt = -1;

	if (data_grp_create(fid, DATA_CIRCCMPSIT) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_DEPTH,
		    cp->dt.depth) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_LENGTH,
		    cp->dt.length) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_STEPSIZE,
		    cp->dt.step_size) < 0)
		goto data_res_write;
	if (data_attr_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_STEPS,
		    cp->dt.steps) < 0)
		goto data_res_write;
	if (data_res_write(fid, DATA_CIRCCMPSIT, DATA_CIRCCMPSIT_VALUES,
		    cp->ct.vals.z, cp->ct.vals.len) < 0)
		goto data_res_write;

	rt = 0;

data_res_write:
	return rt;
}
