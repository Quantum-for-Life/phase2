#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H

#include "data.h"

#define RAYON_NAME "rayon"
#define RAYON_DEFAULT_NUM_MEA_QB 1
//#define RAYON_DEFAULT_NUM_SYS_QB 8
#define RAYON_DEFAULT_NUM_ANC_QB 0


struct rayon_hamil {
        size_t num_qubits;
        size_t num_terms;
        double *coeffs;
        int *paulis;
};

struct rayon_slater_det {
        unsigned long long det;
        double coeff;
};

struct rayon_multidet {
        size_t num_dets;
        struct rayon_slater_det *dets;
};

struct rayon_circuit_data {
        struct rayon_hamil hamil;
        struct rayon_multidet multidet;

};

void rayon_circuit_data_init(struct rayon_circuit_data *ct_dat);

void rayon_circuit_data_destroy(struct rayon_circuit_data *ct_dat);

int rayon_circuit_data_from_data(struct rayon_circuit_data *ct_dat,
                                 const struct data *dat);

int rayon_simulate(struct circ_env env, struct rayon_circuit_data *ct_dat,
                   struct data_time_series *dat_ts);

#endif //PHASE2_RAYON_H
