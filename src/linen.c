/** circuit: linen.
 *
 * Simple circuit implementation for testing.
 */

#include "circ.h"
#include "linen.h"
#include "log/log.h"


int linen_reset(struct circ *c) {
        (void) c;
        log_debug(">>> reset");

        return CIRC_OK;
}

int linen_state_prep(struct circ *c, void *data) {
        (void) c;
        (void) data;
        log_debug(">>> state_prep");

        return CIRC_OK;
}

int linen_routine(struct circ *c, void *data) {
        (void) c;
        (void) data;
        log_debug(">>> routine");

        return CIRC_OK;
}

int linen_state_post(struct circ *c, void *data) {
        (void) c;
        (void) data;
        log_debug(">>> state_post");

        return CIRC_OK;
}

struct circuit linen_circuit_factory(void *data) {
        struct circuit ct = {
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
