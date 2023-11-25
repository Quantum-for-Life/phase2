#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H

#include "data.h"

#define RAYON_NAME "rayon"
#define RAYON_DEFAULT_NUM_MEA_QB 1
#define RAYON_DEFAULT_NUM_SYS_QB 8
#define RAYON_DEFAULT_NUM_ANC_QB 0

struct rayon_hamil {
        size_t num_qubits;
        size_t num_terms;
        double *coeffs;
        int *paulis;
};

void rayon_hamil_init(struct rayon_hamil *hamil);

void rayon_hamil_destroy(struct rayon_hamil *hamil);

int rayon_hamil_from_data(struct rayon_hamil *hamil,
                          const struct data_pauli_hamil *dat_ph);


struct rayon_multidet {
        size_t num_dets;
        struct rayon_slater_det {
                unsigned long long det;
                double coeff;
        } *dets;
};

void rayon_multidet_init(struct rayon_multidet *md);

void rayon_multidet_destroy(struct rayon_multidet *md);

int rayon_multidet_from_data(struct rayon_multidet *md,
                             const struct data_state_prep_multidet *dat_md);


struct rayon_circuit_data {
        struct rayon_hamil hamil;
        struct rayon_multidet multidet;
};

void rayon_circuit_data_init(struct rayon_circuit_data *ct_dat);

void rayon_circuit_data_destroy(struct rayon_circuit_data *ct_dat);

int rayon_circuit_data_from_data(struct rayon_circuit_data *ct_dat,
                                 const struct data_pauli_hamil *dat_ph,
                                 const struct data_state_prep_multidet *dat_md);

struct rayon_circ_data {
        double time;
        int imag_switch;
};

int rayon_state_prep(struct circ *c);

int rayon_routine(struct circ *c);

int rayon_state_post(struct circ *c);

const struct circuit RAYON_CIRCUIT_TEMPLATE = {
        .name = RAYON_NAME,
        .data = NULL,
        .num_mea_qb = RAYON_DEFAULT_NUM_MEA_QB,
        .num_sys_qb = RAYON_DEFAULT_NUM_SYS_QB,
        .num_anc_qb = RAYON_DEFAULT_NUM_ANC_QB,
        .reset = NULL,
        .state_prep = rayon_state_prep,
        .routine = rayon_routine,
        .state_post = rayon_state_post,
};

void rayon_circuit_init(struct circuit *ct,
                        const struct rayon_circuit_data *ct_dat);

void rayon_circuit_destroy(struct circuit *ct) {
        (void) (ct);
}

#endif //PHASE2_RAYON_H
