/** circ: linen.
 *
 * Simple circ implementation for testing.
 */

#include <stdio.h>
#include "circ.h"
#include "circ/linen.h"


int linen_reset(struct circ *c) {
        (void) c;

        return CIRC_OK;
}

int linen_state_prep(struct circ *c) {
        struct linen_circuit_data *ct_dat = c->ct.data;
        struct linen_circ_data *dat = c->data;

        dat->state_prep_value = ct_dat->state_prep_value;

        return CIRC_OK;
}

int linen_state_post(struct circ *c) {
        struct linen_circuit_data *ct_dat = c->ct.data;
        struct linen_circ_data *dat = c->data;

        dat->state_post_value = ct_dat->state_post_value;

        return CIRC_OK;
}

int linen_routine(struct circ *c) {
        struct linen_circuit_data *ct_dat = c->ct.data;
        struct linen_circ_data *dat = c->data;

        dat->routine_value = ct_dat->routine_value;

        return CIRC_OK;
}
