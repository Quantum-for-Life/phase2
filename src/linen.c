/** circuit: linen.
 *
 * Simple circuit implementation for testing.
 */

#include "circ.h"
#include "linen.h"
#include "log/log.h"


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

circuit
linen_circuit_factory(void *data) {
    circuit ct = {
            .name = LINEN_NAME,
            .data = data,
            .num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB,
            .num_sys_qb = LINEN_DEFAULT_NUM_SYS_QB,
            .num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB,
            .reset = linen_reset,
            .state_prep = linen_state_prep,
            .routine = linen_routine,
            .state_post = linen_state_post,
    };

    return ct;
}
