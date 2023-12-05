/** circ: rayon
 *
 *  Quantum phase estimation with Hadamard test.
 */

#include <float.h>

#include "QuEST.h"

#include "log.h"

#include "circ.h"
#include "rayon.h"

struct rayon_circ_data {
	double time;
	int imag_switch;
};

int rayon_state_prep(struct circ *c)
{
	double real, imag;

	const Qureg *qureg = c->qureg;
	const struct rayon_data_multidet *md =
		&((struct rayon_data *)c->ct->data)->multidet;

	initBlankState(*qureg);
	for (size_t i = 0; i < md->num_dets; i++) {
		long long start_idx = md->dets[i].det << c->ct->num_mea_qb;
		real = md->dets[i].coeff_real;
		imag = md->dets[i].coeff_imag;
		setAmps(*qureg, start_idx, &real, &imag, 1);
	}
	hadamard(*qureg, c->mea_qb[0]);

	return 0;
}

int rayon_routine(struct circ *c)
{
	const Qureg *qureg = c->qureg;
	const struct rayon_data_hamil *hamil =
		&((struct rayon_data *)c->ct->data)->hamil;
	enum pauliOpType *paulis = (enum pauliOpType *)hamil->paulis;
	double time = ((struct rayon_circ_data *)c->data)->time;
	int num_mea_qb = (int)c->ct->num_mea_qb;
	int num_sys_qb = (int)c->ct->num_sys_qb;

	if (fabs(time) < DBL_EPSILON) {
		return 0;
	}
	const double REPS = time * time;
	for (size_t r = 0; r < (size_t)REPS; r++) {
		for (size_t i = 0; i < hamil->num_terms; i++) {
			const qreal angle =
				-2.0 * time / REPS * hamil->coeffs[i];
			multiControlledMultiRotatePauli(*qureg, c->mea_qb,
							num_mea_qb, c->sys_qb,
							paulis + num_sys_qb * i,
							num_sys_qb, angle);
		}
	}

	return 0;
}

int rayon_state_post(struct circ *c)
{
	const struct rayon_circ_data *d = c->data;
	const Qureg *qureg = c->qureg;
	if (d->imag_switch == 1) {
		sGate(*qureg, c->mea_qb[0]);
		// pauliZ(*qureg, c->mea_qb[0]);
	}
	hadamard(*qureg, c->mea_qb[0]);

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
	free(hamil->coeffs);
	hamil->coeffs = NULL;
	free(hamil->paulis);
	hamil->paulis = NULL;
	hamil->num_qubits = 0;
	hamil->num_terms = 0;
}

int rayon_hamil_from_data(struct rayon_data_hamil *hamil,
			  const struct data_pauli_hamil *dat_ph)
{
	double *coeffs = malloc(sizeof(*coeffs) * dat_ph->num_terms);
	int *paulis = malloc(sizeof(*paulis) * dat_ph->num_terms *
			     dat_ph->num_qubits);
	if (!(coeffs && paulis)) {
		return -1;
	}

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
	free(md->dets);
	md->dets = NULL;
	md->num_dets = 0;
}

int rayon_multidet_from_data(struct rayon_data_multidet *md,
			     const struct data_state_prep_multidet *dat_md)
{
	md->dets = malloc(sizeof(*md->dets) * dat_md->num_terms);
	if (!md->dets) {
		return -1;
	}
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
	int res = rayon_hamil_from_data(&ct_dat->hamil, &dat->pauli_hamil);
	if (res != 0) {
		return res;
	}
	res = rayon_multidet_from_data(&ct_dat->multidet,
				       &dat->state_prep.multidet);
	if (res != 0) {
		return res;
	}
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
	log_info("Computing expectation values");
	double val[2];

	struct rayon_circ_data *circ_dat = c->data;
	for (size_t i = 0; i < dat_ts->num_steps; i++) {
		double t = dat_ts->times[i];
		val[0] = dat_ts->values[2 * i];
		val[1] = dat_ts->values[2 * i + 1];
		if (!(isnan(val[0]) || isnan(val[1]))) {
			log_trace("time=%f, value already computed; skip", t);
			continue;
		}

		circ_dat->time = t;
		for (int imag_sw = 0; imag_sw <= 1; imag_sw++) {
			circ_dat->imag_switch = imag_sw;
			if (circ_simulate(c) != 0) {
				log_error("Simulation error");
				return -1;
			}
			double prob_0 = c->mea_cl[0] == 0 ?
						c->mea_cl_prob[0] :
						1.0 - c->mea_cl_prob[0];
			val[imag_sw] = 2 * prob_0 - 1;
		}

		dat_ts->values[2 * i] = val[0];
		dat_ts->values[2 * i + 1] = val[1];
		log_trace("t=%f, val_real=%f, val_imag=%f", t, val[0], val[1]);
	}

	return 0;
}

int rayon_simulate(struct circ_env env, const struct rayon_data *ct_dat,
		   const struct data_time_series *dat_ts)
{
	int res;

	log_info("Initialize Pauli Hamiltonian");
	struct circuit ct;
	rayon_circuit_init(&ct, ct_dat);

	log_info("Initialize circ");
	struct rayon_circ_data circ_data = { .imag_switch = 0, .time = 0.0 };
	struct circ c;
	if ((res = circ_init(&c, &env, &ct, &circ_data)) != 0) {
		log_error("Cannot initialize circ");
		goto cleanup;
	}
	res = rayon_compute_expect(&c, dat_ts);
cleanup:
	log_info("Clean up resources");
	circ_destroy(&c);
	rayon_circuit_destroy(&ct);

	return res;
}
