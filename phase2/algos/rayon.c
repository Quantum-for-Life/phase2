/** circ: rayon
 *
 * Quantum phase estimation with Hadamard test.
 * Computes expectation value of `exp(i t H)` for a given initial state,
 * Hamiltonian and a sequence of times.
 */

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "QuEST.h"

#include "circ.h"
#include "algos/rayon.h"

static const size_t PAULI_MASK = 3;
static const size_t PAULI_WIDTH = 2;
static const size_t PAULI_PAK_SIZE = sizeof(pauli_pak_t) * 8 / PAULI_WIDTH;

struct circ_data {
	double t;
	double prob0_re, prob0_im;
	enum pauliOpType scratch[64];
};

void rayon_hamil_init(struct rayon_data_hamil *hamil)
{
	hamil->num_qubits = 0;
	hamil->num_terms = 0;
	hamil->coeffs = NULL;
	hamil->pak = NULL;
}

void rayon_hamil_destroy(struct rayon_data_hamil *hamil)
{
	if (hamil->pak) {
		free(hamil->pak);
		hamil->pak = NULL;
	}
	if (hamil->coeffs) {
		free(hamil->coeffs);
		hamil->coeffs = NULL;
	}
	hamil->num_terms = 0;
	hamil->num_qubits = 0;
}

int rayon_hamil_from_data(struct rayon_data_hamil *hamil,
			  const struct data_pauli_hamil *dat_ph)
{
	double *coeffs = malloc(sizeof(*coeffs) * dat_ph->num_terms);
	pauli_pak_t *pak = calloc(
		dat_ph->num_terms * dat_ph->num_qubits / PAULI_PAK_SIZE + 1,
		sizeof(pauli_pak_t));
	if (!(coeffs && pak)) {
		free(coeffs);
		free(pak);
		return -1;
	}

	for (size_t i = 0; i < dat_ph->num_terms; i++) {
		coeffs[i] = dat_ph->coeffs[i] * dat_ph->norm;
		for (size_t j = 0; j < dat_ph->num_qubits; j++) {
			const size_t pauli_idx = i * dat_ph->num_qubits + j;
			const size_t pak_idx = pauli_idx / PAULI_PAK_SIZE;
			const size_t pak_offset =
				PAULI_WIDTH * (pauli_idx % PAULI_PAK_SIZE);
			const pauli_pak_t pauli = dat_ph->paulis[pauli_idx];
			pak[pak_idx] += pauli << pak_offset;
		}
	}
	hamil->num_qubits = dat_ph->num_qubits;
	hamil->num_terms = dat_ph->num_terms;
	hamil->coeffs = coeffs;
	hamil->pak = pak;

	return 0;
}

void rayon_multidet_init(struct rayon_data_multidet *md)
{
	md->num_dets = 0;
	md->dets = NULL;
}

void rayon_multidet_destroy(struct rayon_data_multidet *md)
{
	if (md->dets) {
		free(md->dets);
		md->dets = NULL;
	}
	md->num_dets = 0;
}

int rayon_multidet_from_data(struct rayon_data_multidet *md,
			     const struct data_state_prep_multidet *dat_md)
{
	md->dets = malloc(sizeof(*md->dets) * dat_md->num_terms);
	if (!md->dets)
		return -1;

	for (size_t i = 0; i < dat_md->num_terms; i++) {
		md->dets[i].coeff = dat_md->coeffs[i];
		const unsigned char *det_seq =
			dat_md->dets + dat_md->num_qubits * i;
		long long index = 0;
		for (size_t j = 0; j < dat_md->num_qubits; j++) {
			index += det_seq[j] << j;
		}
		md->dets[i].index = index;
	}
	md->num_dets = dat_md->num_terms;

	return 0;
}

void rayon_times_init(struct rayon_data_times *ts)
{
	ts->num_steps = 0;
	ts->steps = NULL;
}

void rayon_times_destroy(struct rayon_data_times *ts)
{
	if (ts->steps) {
		free(ts->steps);
		ts->steps = NULL;
	}
	ts->num_steps = 0;
}

int rayon_times_from_data(struct rayon_data_times *ts,
			  const struct data_time_series *dat_ts)
{
	ts->steps = malloc(sizeof(*ts->steps) * dat_ts->num_steps);
	if (!(ts->steps))
		return -1;

	for (size_t i = 0; i < dat_ts->num_steps; i++) {
		ts->steps[i].t = dat_ts->times[i];
		ts->steps[i].val = dat_ts->values[i];
	}
	ts->num_steps = dat_ts->num_steps;

	return 0;
}

void rayon_data_write_times(struct data_time_series *dat,
			    const struct rayon_data_times *rt)
{
	for (size_t i = 0; i < rt->num_steps; i++) {
		dat->values[i] = rt->steps[i].val;
	}
}

void rayon_data_init(struct rayon_data *rd)
{
	rayon_hamil_init(&rd->hamil);
	rayon_multidet_init(&rd->multidet);
	rayon_times_init(&rd->times);
}

void rayon_data_destroy(struct rayon_data *rd)
{
	rayon_times_destroy(&rd->times);
	rayon_multidet_destroy(&rd->multidet);
	rayon_hamil_destroy(&rd->hamil);
}

int rayon_data_from_data(struct rayon_data *rd, const struct data *dat)
{
	int rc = 0;
	rc |= rayon_hamil_from_data(&rd->hamil, &dat->pauli_hamil);
	rc |= rayon_multidet_from_data(&rd->multidet,
				       &dat->state_prep.multidet);
	rc |= rayon_times_from_data(&rd->times, &dat->time_series);

	return rc;
}

int rayon_prepst(struct circ *c)
{
	const Qureg *qureg = c->qureg;
	const struct rayon_data_multidet *md =
		&((struct rayon_data *)c->ct->data)->multidet;

	initBlankState(*qureg);
	for (size_t i = 0; i < md->num_dets; i++) {
		const long long start_idx = md->dets[i].index
					    << c->ct->num_mea_qb;
		double real = creal(md->dets[i].coeff);
		double imag = cimag(md->dets[i].coeff);
		setAmps(*qureg, start_idx, &real, &imag, 1);
	}
	hadamard(*qureg, c->mea_qb[0]);

	return 0;
}

static void trotter_step(const struct circ *c, double omega)
{
	const Qureg *qureg = c->qureg;
	const struct rayon_data_hamil *hamil =
		&((struct rayon_data *)c->ct->data)->hamil;
	struct circ_data *cdat = c->data;
	enum pauliOpType *paulis = cdat->scratch;

	for (size_t i = 0; i < hamil->num_terms; i++) {
		for (size_t j = 0; j < hamil->num_qubits; j++) {
			const size_t pauli_idx = i * hamil->num_qubits + j;
			const size_t pak_idx = pauli_idx / PAULI_PAK_SIZE;
			const size_t pak_offset =
				PAULI_WIDTH * (pauli_idx % PAULI_PAK_SIZE);
			paulis[j] = (hamil->pak[pak_idx] >> pak_offset) &
				    PAULI_MASK;
		}
		/* *
		 * multiControlledMultiRotatePauli() below applies minus
		 * to the given angle.  That sign, together with `sGate` in
		 * rayon_measure()` below (instead of the Hermitian conjugate
		 * of the S gate) gives the correct sign of the imaginary part
		 * of the expectation value.
		 */
		const double angle = 2.0 * omega * hamil->coeffs[i];
		multiControlledMultiRotatePauli(*qureg, c->mea_qb,
						c->ct->num_mea_qb, c->sys_qb,
						paulis, hamil->num_qubits,
						angle);
	}
}

int rayon_effect(struct circ *c)
{
	const double t = ((struct circ_data *)c->data)->t;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	const double REPS = t * t;
	for (size_t r = 0; r < (size_t)REPS; r++) {
		trotter_step(c, t / REPS);
	}

	return 0;
}

int rayon_measure(struct circ *c)
{
	struct circ_data *d = c->data;
	const Qureg *qureg = c->qureg;
	const int mea_qb_idx = c->mea_qb[0];

	hadamard(*qureg, mea_qb_idx);
	d->prob0_re = calcProbOfOutcome(*qureg, mea_qb_idx, 0);

	/* Revert the H gate */
	hadamard(*qureg, mea_qb_idx);
	/**
	 * To obtain the correct sign of the imaginary part of the
	 * expectation value, the gate effected here should be
	 * `S^{\dagger}`.  We use `sGate()` function and change the
	 * sign of the angle argument in `rayon_effect()` instead.
	 */
	sGate(*qureg, mea_qb_idx);
	hadamard(*qureg, mea_qb_idx);
	d->prob0_im = calcProbOfOutcome(*qureg, mea_qb_idx, 0);

	return 0;
}

void rayon_circuit_init(struct circuit *ct, const struct rayon_data *ct_dat)
{
	ct->name = RAYON_NAME;
	ct->data = (void *)ct_dat;
	ct->num_mea_qb = RAYON_NUM_MEA_QB;
	ct->num_sys_qb = ct_dat->hamil.num_qubits;
	ct->num_anc_qb = RAYON_NUM_ANC_QB;
	ct->reset = (circ_op)NULL;
	ct->prepst = rayon_prepst;
	ct->effect = rayon_effect;
	ct->measure = rayon_measure;
}

static void rayon_circuit_destroy(const struct circuit *ct)
{
	(void)(ct);
}

int rayon_simulate(struct circ_env *env, const struct rayon_data *rd)
{
	int ret = 0;

	struct circuit ct;
	rayon_circuit_init(&ct, rd);

	struct circ c;
	struct circ_data cdat;
	if (circ_init(&c, env, &ct, &cdat) < 0)
		goto error;
	const struct rayon_data_times *ts = &rd->times;
	for (size_t i = 0; i < ts->num_steps; i++) {
		cdat.t = ts->steps[i].t;
		if (circ_simulate(&c) < 0)
			goto error;

		ts->steps[i].val = 2.0 * cdat.prob0_re - 1.0 +
				   (2.0 * cdat.prob0_im - 1.0) * _Complex_I;
	}

	goto exit;
error:
	ret = -1;
exit:
	circ_destroy(&c);
	rayon_circuit_destroy(&ct);

	return ret;
}
