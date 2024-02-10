/** circ: silk
 *
 * Quantum phase estimation with Hadamard test.
 */

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "algos/silk.h"
#include "circ.h"

struct circ_data {
	const struct silk_data *rd;
	double			prob0_re, prob0_im;
	int			scratch[64];
};

void silk_multidet_init(struct silk_data_multidet *md)
{
	md->num_dets = 0;
	md->dets     = NULL;
}

void silk_multidet_destroy(struct silk_data_multidet *md)
{
	if (md->dets) {
		free(md->dets);
		md->dets = NULL;
	}
	md->num_dets = 0;
}

struct iter_multidet_data {
	size_t			   idx;
	struct silk_data_multidet *md;
};

static int iter_multidet(_Complex double coeff, size_t idx, void *op_data)
{
	struct iter_multidet_data *imd = op_data;

	imd->md->dets[imd->idx].coeff = coeff;
	imd->md->dets[imd->idx].index = idx;
	imd->idx++;

	return 0;
}

int silk_multidet_from_data(struct silk_data_multidet *md, const data2_id fid)
{
	size_t num_qubits, num_dets;
	if (data2_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		return -1;
	md->dets = malloc(sizeof *md->dets * num_dets);
	if (!md->dets)
		return -1;

	struct iter_multidet_data imd;
	imd.idx = 0;
	imd.md	= md;
	if (data2_multidet_foreach(fid, iter_multidet, &imd) < 0)
		goto error;

	md->num_dets = num_dets;

	return 0;
error:
	free(md->dets);
	return -1;
}

void silk_data_init(struct silk_data *rd, size_t num_steps)
{
	circ_hamil_init(&rd->hamil);
	silk_multidet_init(&rd->multidet);
	rd->num_steps = num_steps;
}

void silk_data_destroy(struct silk_data *rd)
{
	silk_multidet_destroy(&rd->multidet);
	circ_hamil_destroy(&rd->hamil);
}

int silk_data_from_data(struct silk_data *rd, data2_id fid)
{
	int rc;

	rc = circ_hamil_from_data2(&rd->hamil, fid);
	rc |= silk_multidet_from_data(&rd->multidet, fid);

	return rc;
}

int silk_prepst(struct circ *c)
{
	struct circ_data *cdat = circ_data(c);

	const struct silk_data_multidet *md = &cdat->rd->multidet;

	circ_ops_blankstate(c);
	for (size_t i = 0; i < md->num_dets; i++) {
		circ_ops_setsysamp(c, md->dets[i].index, md->dets[i].coeff);
	}
	const qbid mea_qb0 = circ_meaqb(c, 0);
	circ_ops_hadamard(c, mea_qb0);

	return 0;
}

static void trotter_step(struct circ *c, double omega)
{
	struct circ_data	*cdat	= circ_data(c);
	const struct circ_hamil *hamil	= &cdat->rd->hamil;
	int			*paulis = cdat->scratch;

	for (size_t i = 0; i < hamil->num_terms; i++) {
		circ_hamil_paulistr(hamil, i, paulis);
		/* *
		 * The minus sign below, together with `sgate` in
		 * silk_measure()` (instead of the Hermitian conjugate
		 * of the S gate) gives the correct sign of the imaginary part
		 * of the expectation value.
		 */
		const double angle = -1.0 * omega * hamil->coeffs[i];
		circ_ops_crotpauli(c, paulis, angle);
	}
}

int silk_effect(struct circ *c)
{
	struct circ_data *cdat = circ_data(c);
	const double	  t    = cdat->rd->time_factor;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	trotter_step(c, t);

	return 0;
}

int silk_measure(struct circ *c)
{
	struct circ_data *d	  = circ_data(c);
	const qbid	  mea_qb0 = circ_meaqb(c, 0);

	circ_ops_hadamard(c, mea_qb0);
	d->prob0_re = circ_ops_prob0(c, mea_qb0);

	/* Revert the H gate */
	circ_ops_hadamard(c, mea_qb0);
	/**
	 * To obtain the correct sign of the imaginary part of the
	 * expectation value, the gate effected here should be
	 * `S^{\dagger}`.  We use `sgate()` function and change the
	 * sign of the angle argument in `silk_effect()` instead.
	 */
	circ_ops_sgate(c, mea_qb0);
	circ_ops_hadamard(c, mea_qb0);
	d->prob0_im = circ_ops_prob0(c, mea_qb0);

	return 0;
}

static void silk_circuit_init(struct circuit *ct, size_t num_sys_qb)
{
	ct->name       = SILK_NAME;
	ct->num_mea_qb = SILK_NUM_MEA_QB;
	ct->num_sys_qb = num_sys_qb;
	ct->num_anc_qb = SILK_NUM_ANC_QB;
	ct->reset      = NULL;
	ct->prepst     = silk_prepst;
	ct->effect     = silk_effect;
	ct->measure    = silk_measure;
}

static void silk_circuit_destroy(const struct circuit *ct)
{
	(void)(ct);
}

int silk_simulate(const struct silk_data *rd)
{
	int ret = 0;

	struct circuit ct;
	silk_circuit_init(&ct, rd->hamil.num_qubits);

	struct circ_data cdat;
	cdat.rd	       = rd;
	struct circ *c = circ_create(&ct, &cdat);
	if (!c)
		goto error;

	for (size_t i = 0; i < rd->num_steps; i++) {
		if (circ_run(c) < 0)
			goto error;

		rd->trotter_steps[i] = 2.0 * cdat.prob0_re - 1.0 +
				       (2.0 * cdat.prob0_im - 1.0) * _Complex_I;
	}

	goto exit;
error:
	ret = -1;
exit:
	circ_destroy(c);
	silk_circuit_destroy(&ct);

	return ret;
}
