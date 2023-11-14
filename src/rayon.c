/** circuit: rayon
 *
 *  Quantum phase estimation with Hadamard test.
 */

#include "QuEST.h"
#include "circ.h"
#include "rayon.h"

#define RAYON_NAME "rayon"
#define RAYON_DEFAULT_NUM_MEA_QB 1
// #define RAYON_DEFAULT_NUM_SYS_QB 8
#define RAYON_DEFAULT_NUM_ANC_QB 0

circ_result rayon_reset(circ c) {
    (void) c;
    return CIRC_OK;
}

circ_result rayon_state_prep(circ c, void *data) {
    (void) data;

    hadamard(c.qureg, c.mea_qb[0]);

//    pauliX(qureg, sys_qb[0]);
//    pauliX(qureg, sys_qb[2]);
    for (size_t i = 0; i < c.ct.num_sys_qb; i++) {
        hadamard(c.qureg, c.sys_qb[i]);
    }

    return CIRC_OK;
}

circ_result rayon_routine(circ c, void *data) {
    PauliHamil *hamil = (PauliHamil *) circ_circuit_data(c);

    rayon_circ_data *d = (rayon_circ_data *) data;
    double time = d->time;
    for (int i = 0; i < hamil->numSumTerms; i++) {
        qreal angle = 2.0 * time * hamil->termCoeffs[i];
        multiControlledMultiRotatePauli(c.qureg, c.mea_qb, c.ct.num_mea_qb,
                                        c.sys_qb,
                                        hamil->pauliCodes +
                                        (i * hamil->numQubits),
                                        hamil->numQubits, angle);
    }

    return CIRC_OK;
}

circ_result rayon_state_post(circ c, void *data) {


//    pauliX(qureg, sys_qb[0]);
//    pauliX(qureg, sys_qb[2]);

    for (size_t i = 0; i < c.ct.num_sys_qb; i++) {
        hadamard(c.qureg, c.sys_qb[i]);
    }

    rayon_circ_data *d = (rayon_circ_data *) data;
    if (d->imag_switch == 1) {
        sGate(c.qureg, c.mea_qb[0]);
    }
    hadamard(c.qureg, c.mea_qb[0]);

    return CIRC_OK;
}

circuit rayon_circuit_factory(rayon_circuit_data *data) {
    circuit ct = {
            .name = RAYON_NAME,
            .data = data,
            .num_mea_qb = RAYON_DEFAULT_NUM_MEA_QB,
            .num_sys_qb = data->hamil.numQubits,
            .num_anc_qb = RAYON_DEFAULT_NUM_ANC_QB,
            .reset = rayon_reset,
            .state_prep = rayon_state_prep,
            .routine = rayon_routine,
            .state_post = rayon_state_post,
    };
    return ct;
}
