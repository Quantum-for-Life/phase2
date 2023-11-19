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

int linen_state_post(struct circ *c, void *data) {
        (void) c;
        (void) data;
        log_debug(">>> state_post");

        return CIRC_OK;
}

int linen_routine(struct circ *c, void *data) {
        (void) c;
        (void) data;
        log_debug(">>> routine");

        return CIRC_OK;
}
