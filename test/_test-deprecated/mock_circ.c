/*
 * A Mockup of a circuit.
 */

#include "../src/circ.h"

#include "mock_circ.h"


circ_result mock_circ_reset(circ *c) {
    (void) c;
    return CIRC_OK;
}

circ_result mock_circ_state_prep(circ *c, void *data) {
    (void) c;
    struct mock_circ_sample *s = (struct mock_circ_sample *) data;
    s->state = s->input;

    return CIRC_OK;
}

circ_result mock_circ_state_post(circ *c, void *data) {
    (void) c;
    (void) data;
    return CIRC_OK;
}

circ_result mock_circ_routine(circ *c, void *data) {
    (void) c;
    struct mock_circ_sample *s = (struct mock_circ_sample *) data;
    s->state++;

    return CIRC_OK;
}

circ_result mock_circ_measure(circ *c, void *data) {
    (void) c;
    struct mock_circ_sample *s = (struct mock_circ_sample *) data;
    s->result = s->state;

    return CIRC_OK;
}
