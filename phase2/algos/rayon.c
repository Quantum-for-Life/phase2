/** circ: rayon
 *
 * Quantum phase estimation with Hadamard test.
 * Computes expectation value of `exp(i t H)` for a given initial state,
 * Hamiltonian and a sequence of times.
 */

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "circ.h"
#include "algos/rayon.h"

static const size_t PAULI_MASK = 3;
static const size_t PAULI_WIDTH = 2;
static const size_t PAULI_PAK_SIZE = sizeof(pauli_pak_t) * 8 / PAULI_WIDTH;

struct circ_data {
	const struct rayon_data *rd;
	double t;
	double prob0_re, prob0_im;
	int scratch[64];
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
	struct circ_data *cdat = circ_data(c);
	const struct rayon_data_multidet *md =
		&((const struct rayon_data *)cdat->rd)->multidet;

	circ_blankstate(c);
	for (size_t i = 0; i < md->num_dets; i++) {
		circ_setsysamp(c, md->dets[i].index, md->dets[i].coeff);
	}
	const qbid mea_qb0 = circ_mea_qb(c, 0);
	circ_hadamard(c, mea_qb0);

	return 0;
}

static void trotter_step(struct circ *c, double omega)
{
	struct circ_data *cdat = circ_data(c);
	const struct rayon_data_hamil *hamil =
		&((const struct rayon_data *)cdat->rd)->hamil;
	int *paulis = cdat->scratch;

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
		 * The minus sign below, together with `sgate` in
		 * rayon_measure()` (instead of the Hermitian conjugate
		 * of the S gate) gives the correct sign of the imaginary part
		 * of the expectation value.
		 */
		const double angle = -1.0 * omega * hamil->coeffs[i];
		circ_sys_control_rotate_pauli(c, paulis, angle);
	}
}

int rayon_effect(struct circ *c)
{
	struct circ_data *cdat = circ_data(c);
	const double t = ((struct circ_data *)cdat)->t;
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
	struct circ_data *d = circ_data(c);
	const qbid mea_qb0 = circ_mea_qb(c, 0);

	circ_hadamard(c, mea_qb0);
	d->prob0_re = circ_prob0(c, mea_qb0);

	/* Revert the H gate */
	circ_hadamard(c, mea_qb0);
	/**
	 * To obtain the correct sign of the imaginary part of the
	 * expectation value, the gate effected here should be
	 * `S^{\dagger}`.  We use `sgate()` function and change the
	 * sign of the angle argument in `rayon_effect()` instead.
	 */
	circ_sgate(c, mea_qb0);
	circ_hadamard(c, mea_qb0);
	d->prob0_im = circ_prob0(c, mea_qb0);

	return 0;
}

static void rayon_circuit_init(struct circuit *ct, size_t num_sys_qb)
{
	ct->name = RAYON_NAME;
	ct->num_mea_qb = RAYON_NUM_MEA_QB;
	ct->num_sys_qb = num_sys_qb;
	ct->num_anc_qb = RAYON_NUM_ANC_QB;
	ct->reset = NULL;
	ct->prepst = rayon_prepst;
	ct->effect = rayon_effect;
	ct->measure = rayon_measure;
}

static void rayon_circuit_destroy(const struct circuit *ct)
{
	(void)(ct);
}

int rayon_simulate(const struct rayon_data *rd)
{
	int ret = 0;

	struct circuit ct;
	rayon_circuit_init(&ct, rd->hamil.num_qubits);

	struct circ_data cdat;
	cdat.rd = rd;
	struct circ *c = circ_init(&ct, &cdat);
	if (!c)
		goto error;
	const struct rayon_data_times *ts = &rd->times;
	for (size_t i = 0; i < ts->num_steps; i++) {
		cdat.t = ts->steps[i].t;
		if (circ_simulate(c) < 0)
			goto error;

		ts->steps[i].val = 2.0 * cdat.prob0_re - 1.0 +
				   (2.0 * cdat.prob0_im - 1.0) * _Complex_I;
	}

	goto exit;
error:
	ret = -1;
exit:
	circ_destroy(c);
	rayon_circuit_destroy(&ct);

	return ret;
}
