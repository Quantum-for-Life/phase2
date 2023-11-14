/** circuit: linen.
 *
 * Simple circuit implementation for testing.
 */

#include "circ.h"
#include "linen.h"
#include "log/log.h"

#define LINEN_NAME "linen"
#define LINEN_DEFAULT_NUM_MEA_QB 3
#define LINEN_DEFAULT_NUM_SYS_QB 3
#define LINEN_DEFAULT_NUM_ANC_QB 3

circ_result linen_reset(circ c) {
    (void) c;
    log_debug(">>> reset");

    return CIRC_OK;
}

circ_result linen_state_prep(circ c, void *data) {
    (void) c;
    (void) data;
    log_debug(">>> state_prep");

    return CIRC_OK;
}

circ_result linen_routine(circ c, void *data) {
    (void) c;
    (void) data;
    log_debug(">>> routine");

    return CIRC_OK;
}

circ_result linen_state_post(circ c, void *data) {
    (void) c;
    (void) data;
    log_debug(">>> state_post");

    return CIRC_OK;
}

void
linen_circuit_init(circuit *ct, void *data) {
    ct->name = LINEN_NAME;
    ct->data = data;
    ct->num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB;
    ct->num_sys_qb = LINEN_DEFAULT_NUM_SYS_QB;
    ct->num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB;
    ct->reset = linen_reset;
    ct->state_prep = linen_state_prep;
    ct->routine = linen_routine;
    ct->state_post = linen_state_post;
}
