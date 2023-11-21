/** circ: rayon
 *
 *  Quantum phase estimation with Hadamard test.
 */

#include <float.h>

#include "QuEST.h"

#include "circ.h"
#include "circ/rayon.h"


int rayon_state_prep(struct circ* c) {
        hadamard(*c->qureg, c->mea_qb[0]);

        //    pauliX(c.qureg, c.sys_qb[0]);
        //    pauliX(c.qureg, c.sys_qb[2]);
        for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
                hadamard(*c->qureg, c->sys_qb[i]);
        }

        return CIRC_OK;
}

int rayon_routine(struct circ* c) {
        const struct rayon_circuit_data* ct_dat = (struct rayon_circuit_data *)
                c->ct.data;
        const struct rayon_circ_data* dat = (struct rayon_circ_data *)c->data;

        const double time = dat->time;
        if (fabs(time) < DBL_EPSILON) {
                return CIRC_OK;
        }

        const double REPS = time * time;
        for (size_t r = 0; r < (size_t)REPS; r++) {
                for (int i = 0; i < ct_dat->hamil.numSumTerms; i++) {
                        // angle is proportional to time/REPS = 1/time
                        const qreal angle =
                                2.0 / time * ct_dat->hamil.termCoeffs[i];
                        multiControlledMultiRotatePauli(*c->qureg, c->mea_qb,
                                c->ct.num_mea_qb,
                                c->sys_qb,
                                ct_dat->hamil
                                .pauliCodes +
                                i * c->ct.num_sys_qb,
                                c->ct.num_sys_qb,
                                angle);
                }
        }

        return CIRC_OK;
}

int rayon_state_post(struct circ* c) {
        //    pauliX(c.qureg, c.sys_qb[0]);
        //    pauliX(c.qureg, c.sys_qb[2]);
        for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
                hadamard(*c->qureg, c->sys_qb[i]);
        }
        const struct rayon_circ_data* d = (struct rayon_circ_data *)c->data;
        if (d->imag_switch == 1) {
                sGate(*c->qureg, c->mea_qb[0]);
        }
        hadamard(*c->qureg, c->mea_qb[0]);

        return CIRC_OK;
}
