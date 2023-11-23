#ifndef MOCK_CIRC_H
#define MOCK_CIRC_H

#include <stdlib.h>

#include "../src/circ.h"

#define MOCK_CIRC_NAME "mock_circ"
#define MOCK_CIRC_NUM_MEA_CL 3
#define MOCK_CIRC_NUM_MEA_QB 3
#define MOCK_CIRC_NUM_SYS_QB 3
#define MOCK_CIRC_NUM_ANC_QB 2


struct mock_circ_sample {
    int input;
    int state;
    int result;
};

circ_result mock_circ_reset(circ *c);

circ_result mock_circ_state_prep(circ *c, void *data);

circ_result mock_circ_state_post(circ *c, void *data);

circ_result mock_circ_routine(circ *c, void *data);

circ_result mock_circ_measure(circ *c, void *data);

static circuit mock_circ_circuit = {
        .name = MOCK_CIRC_NAME,
        .data = NULL,
        .num_mea_cl = MOCK_CIRC_NUM_MEA_CL,
        .num_mea_qb = MOCK_CIRC_NUM_MEA_QB,
        .num_sys_qb = MOCK_CIRC_NUM_SYS_QB,
        .num_anc_qb = MOCK_CIRC_NUM_ANC_QB,
        .reset = mock_circ_reset,
        .state_prep = mock_circ_state_prep,
        .routine = mock_circ_routine,
        .state_post = mock_circ_state_post,
        .measure = mock_circ_measure};

#endif //MOCK_CIRC_H
