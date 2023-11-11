/** circuit: rayon
 *
 *  Quantum phase estimation with Hadamard test.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "QuEST.h"
#include "circ.h"
#include <hdf5/serial/hdf5.h>
#include "logger.h"

#define RAYON_NAME "rayon"
#define RAYON_DEFAULT_NUM_MEA_CL 1
#define RAYON_DEFAULT_NUM_MEA_QB 1
#define RAYON_DEFAULT_NUM_SYS_QB 8
#define RAYON_DEFAULT_NUM_ANC_QB 0

typedef struct {
    double time;
    int imag_switch;
    int outcome;
    double outcome_prob;
} rayon_CircData;

circ_result rayon_reset(circ *c) {
    Qureg qureg = circ_qureg(c);
    int *mea_cl = circ_mea_cl(c);

    initZeroState(qureg);
    for (int i = 0; i < circ_num_mea_cl(c); i++) {
        mea_cl[i] = 0;
    }

    return CIRC_OK;
}

circ_result rayon_state_prep(circ *c, void *data) {
    (void) data;
    Qureg qureg = circ_qureg(c);
    int *mea_qb = circ_mea_qb(c);
    int *sys_qb = circ_sys_qb(c);

    hadamard(qureg, mea_qb[0]);

//    pauliX(qureg, sys_qb[0]);
//    pauliX(qureg, sys_qb[2]);
    for (int i = 0; i < circ_num_sys_qb(c); i++) {
        hadamard(qureg, sys_qb[i]);
    }

    return CIRC_OK;
}

circ_result rayon_routine(circ *c, void *data) {
    PauliHamil *hamil = (PauliHamil *) circ_circuit_data(c);
    rayon_CircData *circ_data = (rayon_CircData *) data;
    Qureg qureg = circ_qureg(c);
    int *mea_qb = circ_mea_qb(c);
    int *sys_qb = circ_sys_qb(c);

    double time = circ_data->time;
    for (int i = 0; i < hamil->numSumTerms; i++) {
        qreal angle = 2.0 * time * hamil->termCoeffs[i];
        multiControlledMultiRotatePauli(qureg, mea_qb, circ_num_mea_qb(c),
                                        sys_qb,
                                        hamil->pauliCodes + (i *
                                                             hamil->numQubits),
                                        hamil->numQubits, angle);
    }

    return CIRC_OK;
}

circ_result rayon_state_post(circ *c, void *data) {
    rayon_CircData *circ_data = (rayon_CircData *) data;
    Qureg qureg = circ_qureg(c);
    int *mea_qb = circ_mea_qb(c);
    int *sys_qb = circ_sys_qb(c);

//    pauliX(qureg, sys_qb[0]);
//    pauliX(qureg, sys_qb[2]);

    for (int i = 0; i < circ_num_sys_qb(c); i++) {
        hadamard(qureg, sys_qb[i]);
    }

    if (circ_data->imag_switch == 1) {
        sGate(qureg, mea_qb[0]);
    }
    hadamard(qureg, mea_qb[0]);

    return CIRC_OK;
}

circ_result rayon_measure(circ *c, void *data) {
    rayon_CircData *circ_data = (rayon_CircData *) data;
    Qureg qureg = circ_qureg(c);
    int *mea_cl = circ_mea_cl(c);
    int *mea_qb = circ_mea_qb(c);

    int outcome = measureWithStats(qureg, mea_qb[0], &circ_data->outcome_prob);
    mea_cl[0] = outcome;
    circ_data->outcome = outcome;

    return CIRC_OK;
}

circuit rayon_circuit = {
        .name = RAYON_NAME,
        .data = NULL,
        .num_mea_cl = RAYON_DEFAULT_NUM_MEA_CL,
        .num_mea_qb = RAYON_DEFAULT_NUM_MEA_QB,
        .num_sys_qb = RAYON_DEFAULT_NUM_SYS_QB,
        .num_anc_qb = RAYON_DEFAULT_NUM_ANC_QB,
        .reset = rayon_reset,
        .state_prep = rayon_state_prep,
        .routine = rayon_routine,
        .state_post = rayon_state_post,
        .measure = rayon_measure};

int simul_rayon(circ_env *env, char *h5file) {

    logger("info", "open data file");
    hid_t file_id = H5Fopen(h5file, H5F_ACC_RDWR, H5P_DEFAULT);

    logger("info", "reading hamiltonian data");
    hid_t group_id = H5Gopen1(file_id, "hamiltonian");

    hid_t num_qubits_attr_id = H5Aopen(group_id, "num_qubits", H5P_DEFAULT);
    int num_qubits;
    herr_t status = H5Aread(num_qubits_attr_id, H5T_NATIVE_INT, &num_qubits);
    hid_t num_sum_terms_attr_id = H5Aopen(group_id, "num_sum_terms",
                                          H5P_DEFAULT);
    int num_sum_terms;
    status = H5Aread(num_sum_terms_attr_id, H5T_NATIVE_INT, &num_sum_terms);
    status = H5Aclose(num_sum_terms_attr_id);
    status = H5Aclose(num_qubits_attr_id);

    char buf[100];
    snprintf(buf, 100, "num qubits: %d, num_sum_terms: %d", num_qubits,
             num_sum_terms);
    logger("info", buf);

    logger("info", "reading coeffs");
    hid_t dset_coeffs_id = H5Dopen2(group_id, "coeffs", H5P_DEFAULT);
    double *coeffs = (double *) malloc(sizeof(double) * num_sum_terms);
    assert(coeffs != NULL);
    status = H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, coeffs);
    status = H5Dclose(dset_coeffs_id);

    logger("info", "computing hamiltonian one-norm");
    double hamil_one_norm = 0.0;
    for (int i = 0; i < num_sum_terms; i++) {
        hamil_one_norm += fabs(coeffs[i]);
    }
    snprintf(buf, 100, "norm: %f", hamil_one_norm);
    logger("debug", buf);
    hid_t dspace_hamil_one_norm_id = H5Screate(H5S_SCALAR);
    hid_t attr_hamil_one_norm_id = H5Acreate1(group_id, "one_norm",
                                              H5T_IEEE_F64LE,
                                              dspace_hamil_one_norm_id,
                                              H5P_DEFAULT);
    status = H5Awrite(attr_hamil_one_norm_id, H5T_NATIVE_DOUBLE,
                      &hamil_one_norm);
    status = H5Sclose(dspace_hamil_one_norm_id);
    status = H5Aclose(attr_hamil_one_norm_id);

    logger("info", "reading paulis");
    hid_t dset_paulis_id = H5Dopen2(group_id, "paulis", H5P_DEFAULT);
    int *paulis = (int *) malloc(sizeof(int) * num_sum_terms * num_qubits);
    assert(paulis != NULL);
    status = H5Dread(dset_paulis_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, paulis);
    status = H5Dclose(dset_paulis_id);
    status = H5Gclose(group_id);

    logger("debug", "rescaling hamiltonian");
    for (int i = 0; i < num_sum_terms; i++) {
        coeffs[i] /= hamil_one_norm;
    }

    logger("info", "initialize Pauli Hamiltonian");
    PauliHamil hamil = createPauliHamil(num_qubits, num_sum_terms);
    initPauliHamil(hamil, coeffs, (enum pauliOpType *) paulis);
    reportPauliHamil(hamil);
    free(coeffs);
    free(paulis);
    circuit factory = rayon_circuit;
    factory.data = &hamil;
    factory.num_sys_qb = hamil.numQubits;

    logger("info", "read time_series/times");
    hid_t time_series_id = H5Gopen1(file_id, "time_series");
    hid_t steps_attr_id = H5Aopen(time_series_id, "steps", H5P_DEFAULT);
    int steps;
    status = H5Aread(steps_attr_id, H5T_NATIVE_INT, &steps);
    status = H5Aclose(steps_attr_id);
    snprintf(buf, 100, "number of steps: %d", steps);
    logger("info", buf);
    hid_t dset_times_id = H5Dopen1(time_series_id, "times");
    double *times = malloc(sizeof(double) * steps);
    assert(times != NULL);
    status = H5Dread(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, times);
    status = H5Dclose(dset_times_id);

    double *values = malloc(sizeof(double) * steps);
    assert(values != NULL);
    rayon_CircData circ_data = {
            .time = times[0],
            .imag_switch = 0,
            .outcome = 0,
            .outcome_prob = 0.0,
    };
    factory.num_sys_qb = num_qubits;
    circ *circ = circ_create(factory, env, &circ_data);

    logger("info", "evaluating real part of expectation value");
    for (size_t i = 0; i < steps; i++) {
        circ_data.time = times[i];
        circ_reset(circ);
        circ_simulate(circ);

        double prob_0;
        if (circ_data.outcome == 0) {
            prob_0 = circ_data.outcome_prob;
        } else {
            prob_0 = 1.0 - circ_data.outcome_prob;
        }
        values[i] = 2 * prob_0 - 1;
    }
    hid_t dspace_values_id = H5Screate_simple(1, (hsize_t[]) {steps}, NULL);
    hid_t dset_values_id = H5Dcreate1(time_series_id, "values_real",
                                      H5T_IEEE_F64LE, dspace_values_id,
                                      H5P_DEFAULT);
    status = H5Dwrite(dset_values_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                      dspace_values_id, H5P_DEFAULT, values);
    status = H5Sclose(dspace_values_id);
    status = H5Dclose(dset_values_id);

    logger("info", "evaluating imag part of expectation value");
    circ_data.imag_switch = 1;
    for (size_t i = 0; i < steps; i++) {
        circ_data.time = times[i];
        circ_reset(circ);
        circ_simulate(circ);

        double prob_0;
        if (circ_data.outcome == 0) {
            prob_0 = circ_data.outcome_prob;
        } else {
            prob_0 = 1.0 - circ_data.outcome_prob;
        }
        values[i] = 2 * prob_0 - 1;
    }

    dspace_values_id = H5Screate_simple(1, (hsize_t[]) {steps}, NULL);
    dset_values_id = H5Dcreate1(time_series_id, "values_imag",
                                H5T_IEEE_F64LE, dspace_values_id,
                                H5P_DEFAULT);
    status = H5Dwrite(dset_values_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                      dspace_values_id, H5P_DEFAULT, values);
    status = H5Sclose(dspace_values_id);
    status = H5Dclose(dset_values_id);

    free(times);
    free(values);
    circ_destroy(circ);

    status = H5Gclose(time_series_id);
    status = H5Fclose(file_id);

    return EXIT_SUCCESS;
}
