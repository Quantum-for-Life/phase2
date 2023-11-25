/** circ: rayon
 *
 *  Quantum phase estimation with Hadamard test.
 */

#include <float.h>

#include "QuEST.h"

#include "log.h"

#include "circ.h"
#include "circ/rayon.h"


int rayon_state_prep(struct circ *c) {
        Qureg *qureg = c->qureg;
        hadamard(*qureg, c->mea_qb[0]);

        //    pauliX(c.qureg, c.sys_qb[0]);
        //    pauliX(c.qureg, c.sys_qb[2]);
        for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
                hadamard(*qureg, c->sys_qb[i]);
        }

        return CIRC_OK;
}

int rayon_routine(struct circ *c) {
        Qureg *qureg = c->qureg;

        const struct rayon_circuit_data *ct_dat = (struct rayon_circuit_data *)
                c->ct.data;
        const struct rayon_circ_data *dat = (struct rayon_circ_data *) c->data;

        const double time = dat->time;
        if (fabs(time) < DBL_EPSILON) {
                return CIRC_OK;
        }

        const double REPS = time * time;
        for (size_t r = 0; r < (size_t) REPS; r++) {
                for (size_t i = 0; i < ct_dat->hamil.num_terms; i++) {
                        // angle is proportional to time/REPS = 1/time
                        const qreal angle =
                                2.0 / time * ct_dat->hamil.coeffs[i];
                        multiControlledMultiRotatePauli(*qureg, c->mea_qb,
                                                        c->ct.num_mea_qb,
                                                        c->sys_qb,
                                                        (enum pauliOpType *)
                                                                ct_dat->hamil
                                                                        .paulis +
                                                        c->ct.num_sys_qb * i,
                                                        c->ct.num_sys_qb,
                                                        angle);
                }
        }

        return CIRC_OK;
}

int rayon_state_post(struct circ *c) {
        Qureg *qureg = c->qureg;

        //    pauliX(c.qureg, c.sys_qb[0]);
        //    pauliX(c.qureg, c.sys_qb[2]);
        for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
                hadamard(*qureg, c->sys_qb[i]);
        }
        const struct rayon_circ_data *d = (struct rayon_circ_data *) c->data;
        if (d->imag_switch == 1) {
                sGate(*qureg, c->mea_qb[0]);
        }
        hadamard(*qureg, c->mea_qb[0]);

        return CIRC_OK;
}

void rayon_hamil_init(struct rayon_hamil *hamil) {
        hamil->num_qubits = 0;
        hamil->num_terms = 0;
        hamil->coeffs = NULL;
        hamil->paulis = NULL;
}

void rayon_hamil_destroy(struct rayon_hamil *hamil) {
        free(hamil->coeffs);
        hamil->coeffs = NULL;
        free(hamil->paulis);
        hamil->paulis = NULL;
        hamil->num_qubits = 0;
        hamil->num_terms = 0;
}

int
rayon_hamil_from_data(struct rayon_hamil *hamil,
                      const struct data_pauli_hamil *dat_ph) {

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

void rayon_multidet_init(struct rayon_multidet *md) {
        md->num_dets = 0;
        md->dets = NULL;
}

void rayon_multidet_destroy(struct rayon_multidet *md) {
        free(md->dets);
        md->dets = NULL;
        md->num_dets = 0;
}

int rayon_multidet_from_data(struct rayon_multidet *md,
                             const struct data_state_prep_multidet *dat_md) {

        struct rayon_slater_det *dets = malloc(
                sizeof(*dets) * dat_md->num_terms);
        if (!dets) {
                return -1;
        }
        for (size_t i = 0; i < dat_md->num_terms; i++) {
                dets[i].coeff = dat_md->coeffs[i];
                unsigned long long index = 0;
                for (size_t j = 0; j < dat_md->num_qubits; j++) {
                        unsigned long long bit =
                                dat_md->dets[i * dat_md->num_qubits + j] == 0 ?
                                0 : 1;
                        index += bit << j;
                }
                dets[i].det = index;
        }

        md->num_dets = dat_md->num_terms;
        md->dets = dets;

        return 0;
}


void rayon_circuit_data_init(struct rayon_circuit_data *ct_dat) {
        rayon_hamil_init(&ct_dat->hamil);
        rayon_multidet_init(&ct_dat->multidet);
}

void rayon_circuit_data_destroy(struct rayon_circuit_data *ct_dat) {
        rayon_hamil_destroy(&ct_dat->hamil);
        rayon_multidet_destroy(&ct_dat->multidet);
}

int rayon_circuit_data_from_data(struct rayon_circuit_data *ct_dat,
                                 const struct data *dat) {
        int res;
        if ((res = rayon_hamil_from_data(&ct_dat->hamil, &dat->pauli_hamil))
            != 0) {
                return res;
        }
        if ((res = rayon_multidet_from_data(&ct_dat->multidet,
                                            &dat->state_prep.multidet)) != 0) {
                return res;
        }

        return 0;
}

void rayon_circuit_init(struct circuit *ct,
                        const struct rayon_circuit_data *ct_dat) {
        ct->name = RAYON_NAME;
        ct->data = (void *) ct_dat;
        ct->num_mea_qb = RAYON_DEFAULT_NUM_MEA_QB;
        ct->num_sys_qb = ct_dat->hamil.num_qubits;
        ct->num_anc_qb = RAYON_DEFAULT_NUM_ANC_QB;
        ct->reset = NULL;
        ct->state_prep = rayon_state_prep;
        ct->routine = rayon_routine;
        ct->state_post = rayon_state_post;
}

void rayon_circuit_destroy(struct circuit *ct) {
        (void) (ct);
}


int rayon_simulate(struct circ_env env, const struct rayon_circuit_data *ct_dat,
                   struct data_time_series *dat_ts) {
        log_info("Initialize Pauli Hamiltonian");

        struct circuit ct;
        rayon_circuit_init(&ct, ct_dat);

        log_info("Initialize circ");
        struct rayon_circ_data circ_data = {.imag_switch = 0, .time = 0.0};
        struct circ c;
        if (circ_init(&c, env, ct, &circ_data) != CIRC_OK) {
                log_error("Cannot initialize circ");
                return CIRC_ERR;
        }

        log_info("Computing expectation values");
        for (size_t i = 0; i < dat_ts->num_steps; i++) {

                if (!isnan(dat_ts->values[2 * i]) &&
                    !isnan(dat_ts->values[2 * i + 1])) {
                        continue;
                }

                circ_data.imag_switch = 0;
                while (circ_data.imag_switch <= 1) {
                        const size_t offset =
                                circ_data.imag_switch == 0 ? 0 : 1;
                        circ_data.time = dat_ts->times[i];
                        if (circ_simulate(&c) != CIRC_OK) {
                                log_error("Simulation error");
                                rayon_circuit_destroy(&ct);
                                return CIRC_ERR;
                        }
                        const double prob_0 = c.mea_cl[0] == 0
                                              ? c.mea_cl_prob[0]
                                              : 1.0 - c.mea_cl_prob[0];
                        dat_ts->values[2 * i + offset] =
                                2 * prob_0 - 1;
                        circ_data.imag_switch++;
                }

        }
        circ_destroy(&c);
        log_info("End of simulation");
        rayon_circuit_destroy(&ct);

        return CIRC_OK;
}
