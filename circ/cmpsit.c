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

static void cmpsit_hamil_rearrange(struct cmpsit *cp)
{
	struct circ_hamil *hm = &cp->ct.hm;
	qsort(hm->terms, hm->len, sizeof(struct circ_hamil_term),
		hamil_term_cmp_abscf_desc);
}

static int cmpsit_ranct_init(struct cmpsit_ranct *rct, const uint32_t qb,
	const size_t hm_ran_len, const size_t cdf_len)
{
	if (prob_cdf_init(&rct->cdf, cdf_len) < 0)
		return -1;
	if (circ_hamil_init(&rct->hm_ran, qb, hm_ran_len) < 0)
		return -1;

	return 0;
}

static void cmpsit_ranct_free(struct cmpsit_ranct *rct)
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

static void cmpsit_ranct_calc_cdf(
	struct cmpsit_ranct *rct, struct circ_hamil_term *terms)
{
	struct get_vals_data data = { .i = 0, .terms = terms };
	prob_cdf_from(&rct->cdf, get_vals, &data);
}

int cmpsit_init(
	struct cmpsit *cp, const struct cmpsit_data *dt, const data_id fid)
{
	if (circ_init(&cp->ct, fid, dt->samples) < 0)
		goto err_circ_init;

	cp->dt = *dt;

	cmpsit_hamil_rearrange(cp);

	/* Calculate the sampled circuit size (no. of pauli rotations). */
	const size_t hm_ran_len = dt->length * dt->depth * CMPSIT_TRUNC_DIST;
	if (cmpsit_ranct_init(&cp->ranct, cp->ct.hm.qb, hm_ran_len,
		    cp->ct.hm.len - dt->depth) < 0)
		goto err_ranct_init;
	cmpsit_ranct_calc_cdf(&cp->ranct, cp->ct.hm.terms + dt->depth);

	/* TODO: seed it with user-supplied seed */
	xoshiro256ss_init(&cp->rng, SEED);

	return 0;

	// cmpsit_ranct_free(&cp->ranct);
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

static void cmpsit_ranct_sample(struct cmpsit *cp)
{
}

int cmpsit_simul(struct cmpsit *cp)
{
	/* Second order Trotter */
	struct circ *ct = &cp->ct;
	struct circ_values *vals = &ct->vals;

	struct circ_prog prog;
	circ_prog_init(&prog, vals->len);
	for (size_t i = 0; i < vals->len; i++) {
		cmpsit_ranct_sample(cp);

		circ_prepst(ct);
		// if (circ_step(&cp->ct, &cp) < 0)
		//	return -1;
		vals->z[i] = circ_measure(ct);

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
