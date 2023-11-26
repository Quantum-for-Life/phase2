#ifndef PHASE2_LINEN_H
#define PHASE2_LINEN_H

#include "circ.h"

#define LINEN_NAME "linen"
#define LINEN_DEFAULT_NUM_MEA_QB 4
#define LINEN_DEFAULT_NUM_SYS_QB 8
#define LINEN_DEFAULT_NUM_ANC_QB 4

struct linen_circuit_data {
        int state_prep_value;
        int routine_value;
        int state_post_value;
};

struct linen_circ_data {
        int state_prep_value;
        int routine_value;
        int state_post_value;
};

void linen_circuit_init(struct circuit *ct, struct linen_circuit_data *ct_dat);

int linen_simulate(struct circ_env env);

#endif //PHASE2_LINEN_H
