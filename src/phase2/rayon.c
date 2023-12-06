/** circ: rayon
 * Quantum phase estimation with Hadamard test.
 * Computes expectation value of `exp(i t H)` for a given initial state,
 * Hamiltonian and a sequence of times.
 */

#include <float.h>
#include <stdio.h>

#include "QuEST.h"

#include "circ.h"
#include "rayon.h"

struct rayon_circ_data {
	double time;
	double prob0_re;
	double prob0_im;
};

int rayon_state_prep(struct circ *c)
{
	const Qureg *qureg = c->qureg;
	const struct rayon_data_multidet *md =
		&((struct rayon_data *)c->ct->data)->multidet;

	initBlankState(*qureg);
	for (size_t i = 0; i < md->num_dets; i++) {
		long long start_idx = md->dets[i].det << c->ct->num_mea_qb;
		double real = md->dets[i].coeff_real;
		double imag = md->dets[i].coeff_imag;
		setAmps(*qureg, start_idx, &real, &imag, 1);
	}
	hadamard(*qureg, c->mea_qb[0]);

	return 0;
}

static void trotter_step(struct circ *c, double omega)
{
	const Qureg *qureg = c->qureg;
	const struct rayon_data_hamil *hamil =
		&((struct rayon_data *)c->ct->data)->hamil;
	enum pauliOpType *paulis = (enum pauliOpType *)hamil->paulis;
	int num_mea_qb = (int)c->ct->num_mea_qb;
	int num_sys_qb = (int)c->ct->num_sys_qb;

	for (size_t i = 0; i < hamil->num_terms; i++) {
		/* *
		 * multiControlledMultiRotatePauli() below applies minus
		 * to the given angle.  That sign, together with `sGate` in
		 * rayon_state_post()` below (instead of the Hermitian conjugate
		 * of the S gate) gives the correct sign of the imaginary part
		 * of the expectation value.
		 */
		double angle = 2.0 * omega * hamil->coeffs[i];
		multiControlledMultiRotatePauli(*qureg, c->mea_qb, num_mea_qb,
						c->sys_qb,
						paulis + num_sys_qb * i,
						num_sys_qb, angle);
	}
}

int rayon_routine(struct circ *c)
{
	double t = ((struct rayon_circ_data *)c->data)->time;
	if (isnan(t) || isinf(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	const double REPS = t * t;
	for (size_t r = 0; r < (size_t)REPS; r++) {
		trotter_step(c, t / REPS);
	}

	return 0;
}

int rayon_state_post(struct circ *c)
{
	struct rayon_circ_data *d = c->data;
	const Qureg *qureg = c->qureg;

	hadamard(*qureg, c->mea_qb[0]);
	d->prob0_re = calcProbOfOutcome(*qureg, c->mea_qb[0], 0);

	/* Revert the H gate */
	hadamard(*qureg, c->mea_qb[0]);
	
	/**
	 * To obtain the correct sign of the imaginary part of the
	 * expectation value, the gate effected here should be
	 * `S^{\dagger}`.  We use `sGate()` function and change the
	 * sign of the angle argument in `rayon_routine()` instead.
	 */
	sGate(*qureg, c->mea_qb[0]);
	hadamard(*qureg, c->mea_qb[0]);
	d->prob0_im = calcProbOfOutcome(*qureg, c->mea_qb[0], 0);

	return 0;
}

void rayon_hamil_init(struct rayon_data_hamil *hamil)
{
	hamil->num_qubits = 0;
	hamil->num_terms = 0;
	hamil->coeffs = NULL;
	hamil->paulis = NULL;
}

void rayon_hamil_destroy(struct rayon_data_hamil *hamil)
{
	if (hamil->paulis) {
		free(hamil->paulis);
		hamil->paulis = NULL;
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
	int *paulis = malloc(sizeof(*paulis) * dat_ph->num_terms *
			     dat_ph->num_qubits);
	if (!(coeffs && paulis))
		return -1;

	for (size_t i = 0; i < dat_ph->num_terms; i++) {
		coeffs[i] = dat_ph->coeffs[i] * dat_ph->norm;
		for (size_t j = 0; j < dat_ph->num_qubits; j++) {
			paulis[i * dat_ph->num_qubits + j] =
				dat_ph->paulis[i * dat_ph->num_qubits + j];
		}
	}
	hamil->num_qubits = dat_ph->num_qubits;
	hamil->num_terms = dat_ph->num_terms;
	hamil->coeffs = coeffs;
	hamil->paulis = paulis;

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
		md->dets[i].coeff_real = dat_md->coeffs[2 * i];
		md->dets[i].coeff_imag = dat_md->coeffs[2 * i + 1];
		long long index = 0;
		for (size_t j = 0; j < dat_md->num_qubits; j++) {
			long long bit =
				dat_md->dets[i * dat_md->num_qubits + j] == 0 ?
					0 :
					1;
			index += bit << j;
		}
		md->dets[i].det = index;
	}
	md->num_dets = dat_md->num_terms;

	return 0;
}

void rayon_data_init(struct rayon_data *ct_dat)
{
	rayon_hamil_init(&ct_dat->hamil);
	rayon_multidet_init(&ct_dat->multidet);
}

void rayon_data_destroy(struct rayon_data *ct_dat)
{
	rayon_multidet_destroy(&ct_dat->multidet);
	rayon_hamil_destroy(&ct_dat->hamil);
}

int rayon_data_from_data(struct rayon_data *ct_dat, const struct data *dat)
{
	if (rayon_hamil_from_data(&ct_dat->hamil, &dat->pauli_hamil) < 0)
		return -1;
	if (rayon_multidet_from_data(&ct_dat->multidet,
				     &dat->state_prep.multidet) < 0)
		return -1;

	return 0;
}

void rayon_circuit_init(struct circuit *ct, const struct rayon_data *ct_dat)
{
	ct->name = RAYON_NAME;
	ct->data = (void *)ct_dat;
	ct->num_mea_qb = RAYON_DEFAULT_NUM_MEA_QB;
	ct->num_sys_qb = ct_dat->hamil.num_qubits;
	ct->num_anc_qb = RAYON_DEFAULT_NUM_ANC_QB;
	ct->reset = NULL;
	ct->state_prep = rayon_state_prep;
	ct->routine = rayon_routine;
	ct->state_post = rayon_state_post;
}

static void rayon_circuit_destroy(struct circuit *ct)
{
	(void)(ct);
}

int rayon_compute_expect(struct circ *c, const struct data_time_series *dat_ts)
{
	struct rayon_circ_data *cdat = c->data;
	for (size_t i = 0; i < dat_ts->num_steps; i++) {
		cdat->time = dat_ts->times[i];
		if (circ_simulate(c) < 0)
			return -1;

		dat_ts->values[2 * i] = 2.0 * cdat->prob0_re - 1.0;
		dat_ts->values[2 * i + 1] = 2.0 * cdat->prob0_im - 1.0;
	}

	return 0;
}

int rayon_simulate(struct circ_env *env, const struct rayon_data *ct_dat,
		   const struct data_time_series *dat_ts)
{
	int rc;

	struct circuit ct;
	rayon_circuit_init(&ct, ct_dat);

	struct rayon_circ_data cdat;
	struct circ c;
	rc = circ_init(&c, env, &ct, &cdat);
	if (rc < 0)
		goto cleanup;
	rc = rayon_compute_expect(&c, dat_ts);
cleanup:
	circ_destroy(&c);
	rayon_circuit_destroy(&ct);

	return rc;
}
