/** circuit: linen.
 *
 * Simple circuit implementation for testing.
 */

#include <stdio.h>

#include "circ.h"
#include "log/log.h"

#define LINEN_NAME "linen"
#define LINEN_DEFAULT_NUM_MEA_QB 3
#define LINEN_DEFAULT_NUM_SYS_QB 3
#define LINEN_DEFAULT_NUM_ANC_QB 3

circ_result linen_reset(circ *c) {
    (void) c;
    log_debug(">>> reset");
    return CIRC_OK;
}

circ_result linen_state_prep(circ *c, void *data) {
    (void) c;
    (void) data;
    log_debug(">>> state_prep");
    return CIRC_OK;
}

circ_result linen_routine(circ *c, void *data) {
    (void) c;
    (void) data;
    log_debug(">>> routine");
    return CIRC_OK;
}

circ_result linen_state_post(circ *c, void *data) {
    (void) c;
    (void) data;
    log_debug(">>> state_post");
    return CIRC_OK;
}

circ_result linen_measure(circ *c, void *data) {
    (void) c;
    (void) data;
    log_debug(">>> measure");
    return CIRC_OK;
}

circuit linen_circuit(circuit_data data) {
    circuit ct = {
            .name = LINEN_NAME,
            .data = data,
            .num_mea_qb = LINEN_DEFAULT_NUM_MEA_QB,
            .num_sys_qb = data.hamil.numQubits,
            .num_anc_qb = LINEN_DEFAULT_NUM_ANC_QB,
            .reset = linen_reset,
            .state_prep = linen_state_prep,
            .routine = linen_routine,
            .state_post = linen_state_post,
    };
    return ct;
}

int
linen_simulate(circ_env *env, circuit_data data) {

    log_debug("Report simulation environment");
    circ_report_env(env);

    circuit factory = linen_circuit(data);
    circ *circ = circ_create(factory, env, NULL);
    if (circ == NULL) {
        log_error("Cannot initialize circuit");
        return -1;
    }
    log_debug("\"linen\" circuit created");
    circ_report(circ);
    log_debug("Simulating circuit");
    circ_simulate(circ);
    log_debug("Free circuit instance");
    circ_destroy(circ);

    return 0;
}