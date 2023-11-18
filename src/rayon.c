/** circuit: rayon
 *
 *  Quantum phase estimation with Hadamard test.
 */

#include <float.h>
#include "QuEST.h"
#include "circ.h"
#include "rayon.h"

#define RAYON_NAME "rayon"
#define RAYON_DEFAULT_NUM_MEA_QB 1
// #define RAYON_DEFAULT_NUM_SYS_QB 8
#define RAYON_DEFAULT_NUM_ANC_QB 0


int rayon_state_prep(struct circ *c, void *data) {
    (void) data;

    hadamard(*c->qureg, c->mea_qb[0]);

    //    pauliX(c.qureg, c.sys_qb[0]);
    //    pauliX(c.qureg, c.sys_qb[2]);
    for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
        hadamard(*c->qureg, c->sys_qb[i]);
    }

    return CIRC_OK;
}

int rayon_routine(struct circ *c, void *data) {
    struct rayon_circuit_data *ctdat = (struct rayon_circuit_data *) c->ct.data;
    struct rayon_circ_data *dat = (struct rayon_circ_data *) data;

    double time = dat->time;
    if (fabs(time) < DBL_EPSILON) {
        return CIRC_OK;
    }

    const double REPS = time * time;
    for (size_t r = 0; r < (size_t) REPS; r++) {
        for (int i = 0; i < ctdat->hamil.numSumTerms; i++) {
            // angle is proportional to time/REPS = 1/time
            qreal angle = 2.0 / time * ctdat->hamil.termCoeffs[i];
            multiControlledMultiRotatePauli(*c->qureg, c->mea_qb,
                                            c->ct.num_mea_qb,
                                            c->sys_qb,
                                            ctdat->hamil.pauliCodes +
                                            (i * c->ct.num_sys_qb),
                                            c->ct.num_sys_qb, angle);
        }
    }

    return CIRC_OK;
}

int rayon_state_post(struct circ *c, void *data) {
    //    pauliX(c.qureg, c.sys_qb[0]);
    //    pauliX(c.qureg, c.sys_qb[2]);
    for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
        hadamard(*c->qureg, c->sys_qb[i]);
    }
    struct rayon_circ_data *d = (struct rayon_circ_data *) data;
    if (d->imag_switch == 1) {
        sGate(*c->qureg, c->mea_qb[0]);
    }
    hadamard(*c->qureg, c->mea_qb[0]);

    return CIRC_OK;
}

struct circuit rayon_circuit_factory(struct rayon_circuit_data *data) {
    struct circuit ct = {
            .name = RAYON_NAME,
            .data = data,
            .num_mea_qb = RAYON_DEFAULT_NUM_MEA_QB,
            .num_sys_qb = data->hamil.numQubits,
            .num_anc_qb = RAYON_DEFAULT_NUM_ANC_QB,
            .reset = NULL,
            .state_prep = rayon_state_prep,
            .routine = rayon_routine,
            .state_post = rayon_state_post,
    };

    return ct;
}
