#ifndef PHASE2_RAYON_H
#define PHASE2_RAYON_H

#define RAYON_NAME "rayon"
#define RAYON_DEFAULT_NUM_MEA_QB 1
#define RAYON_DEFAULT_NUM_SYS_QB 8
#define RAYON_DEFAULT_NUM_ANC_QB 0

struct rayon_circuit_data {
        PauliHamil hamil;
        void *data; // state preparation. TBA
};

struct rayon_circ_data {
        double time;
        int imag_switch;
};

int rayon_state_prep(struct circ *c, void *data);

int rayon_routine(struct circ *c, void *data);

int rayon_state_post(struct circ *c, void *data);

struct circuit rayon_circuit = {
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


#endif //PHASE2_RAYON_H
