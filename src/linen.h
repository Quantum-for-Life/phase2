#ifndef PHASE2_LINEN_H
#define PHASE2_LINEN_H

#include "circ.h"

#define LINEN_NAME "linen"
#define LINEN_DEFAULT_NUM_MEA_QB 4
#define LINEN_DEFAULT_NUM_SYS_QB 8
#define LINEN_DEFAULT_NUM_ANC_QB 4


int linen_reset(struct circ *c);

int linen_state_prep(struct circ *c, void *data);

int linen_routine(struct circ *c, void *data);

int linen_state_post(struct circ *c, void *data);


struct circuit linen_circuit = {
        .name = LINEN_NAME,
        .data = NULL,
        .num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB,
        .num_sys_qb = LINEN_DEFAULT_NUM_SYS_QB,
        .num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB,
        .reset = linen_reset,
        .state_prep = linen_state_prep,
        .routine = linen_routine,
        .state_post = linen_state_post,
};


#endif //PHASE2_LINEN_H
