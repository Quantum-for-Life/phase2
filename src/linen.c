/** circuit: linen.
 *
 * Simple circuit implementation for testing.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "circ.h"
#include "logger.h"

#define LINEN_NAME "linen"
#define LINEN_DEFAULT_NUM_MEA_CL 3
#define LINEN_DEFAULT_NUM_MEA_QB 3
#define LINEN_DEFAULT_NUM_SYS_QB 3
#define LINEN_DEFAULT_NUM_ANC_QB 3

circ_result linen_reset(circ *c) {
    (void) c;
    logger("debug", "reset circuit");
    return CIRC_OK;
}

circ_result linen_state_prep(circ *c, void *data) {
    (void) c;
    (void) data;
    logger("debug", "state_prep");
    return CIRC_OK;
}

circ_result linen_routine(circ *c, void *data) {
    (void) c;
    (void) data;
    logger("debug", "routine");
    return CIRC_OK;
}

circ_result linen_state_post(circ *c, void *data) {
    (void) c;
    (void) data;
    logger("debug", "state_post");
    return CIRC_OK;
}

circ_result linen_measure(circ *c, void *data) {
    (void) c;
    (void) data;
    logger("debug", "measure");
    return CIRC_OK;
}

circuit linen_circuit = {
        .name = LINEN_NAME,
        .data = NULL,
        .num_mea_cl = LINEN_DEFAULT_NUM_MEA_CL,
        .num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB,
        .num_sys_qb = LINEN_DEFAULT_NUM_SYS_QB,
        .num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB,
        .reset = linen_reset,
        .state_prep = linen_state_prep,
        .routine = linen_routine,
        .state_post = linen_state_post,
        .measure = linen_measure};


int simul_linen(circ_env *env) {

    circuit factory = linen_circuit;

    circ *circ = circ_create(factory, env, NULL);
    if (circ) {
        logger("info", "linen circ created");
    }

    logger("debug", "free circuit instance");
    circ_destroy(circ);

    return EXIT_SUCCESS;
}