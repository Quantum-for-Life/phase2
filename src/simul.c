#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"

#include "log/log.h"
#include "circ.h"
#include "rayon.h"

#define SIMUL_DEFAULT_H5FILE "simul.h5"

circ_result linen_simulate(circ_env *env, void *data);

circuit linen_circuit(void *data);

circ_result rayon_simulate(circ_env *env, const char *hamil_file);

circuit rayon_circuit(void *data);


void exit_failure() {
    log_error("Failed.");
    exit(EXIT_FAILURE);
}

void print_help_page(int argc, char **argv) {
    (void) argc;
    fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5]\n\n", argv[0]);
    fprintf(stderr,
            "If no simulation input file (HDF5) is specified,\n"
            "the default is " SIMUL_DEFAULT_H5FILE " in the current directory.\n"
    );
}

int main(int argc, char **argv) {

    log_info("*** Init ***");
    log_info("Open log file: simul.log");
    FILE *log_file = fopen("simul.log", "a");
    if (log_file == NULL) {
        log_error("Cannot open log file.");
        exit_failure();
    }
    log_add_fp(log_file, LOG_DEBUG);


    log_info("Initialize simulation environment");
    circ_env *env = circ_create_env();
    if (env == NULL) {
        exit_failure();
    }

    log_debug("Parsing command line arguments");
    if (argc < 2) {
        log_error("No circuit name specified");
        print_help_page(argc, argv);
        exit_failure();
    }
    const char *h5filename = SIMUL_DEFAULT_H5FILE;
    if (argc < 3) {
        log_debug("No simulation input file specified; "
                  "using default: %s",
                  SIMUL_DEFAULT_H5FILE);
    } else {
        h5filename = argv[2];
    }

    log_info("Open simulation input file: %s", h5filename);


    log_info("*** Circuit ***");
    // TODO: Initialize hamiltonian with data from file
    if (strncmp(argv[1], "linen", 5) == 0) {
        log_info("Circuit: linen");
        if (linen_simulate(env, NULL) != 0) {
            exit_failure();
        }
    }

    if (strncmp(argv[1], "rayon", 5) == 0) {
        log_info("Circuit: rayon");
        if (rayon_simulate(env, h5filename) != 0) {
            exit_failure();
        }
    }

    log_info("*** Cleanup ***");
    log_info("Shut down simulation environment");
    circ_destroy_env(env);

    log_info("Done");
    fclose(log_file);
    return EXIT_SUCCESS;
}


circ_result
linen_simulate(circ_env *env, void *data) {

    log_debug("Report simulation environment");
    circ_report_env(env);
    circuit factory = linen_circuit(data);
    circ *circ = circ_create(factory, env, NULL);
    if (circ == NULL) {
        log_error("Cannot initialize circuit");
        return CIRC_ERR;
    }
    log_debug("\"linen\" circuit created");
    circ_report(circ);
    log_debug("Simulating circuit");
    circ_simulate(circ);
    log_debug("Free circuit instance");
    circ_destroy(circ);

    return CIRC_OK;
}


circ_result rayon_simulate(circ_env *env, const char *hamil_file) {

    log_info("open data file");
    hid_t file_id = H5Fopen(hamil_file, H5F_ACC_RDWR, H5P_DEFAULT);

    log_info("reading hamiltonian data");
    hid_t group_id = H5Gopen1(file_id, "hamiltonian");

    hid_t num_qubits_attr_id = H5Aopen(group_id, "num_qubits", H5P_DEFAULT);
    int num_qubits;
    herr_t status;
    if ((status = H5Aread(num_qubits_attr_id, H5T_NATIVE_INT, &num_qubits)) <
        0) {
        log_error("no field: num_qubits");
        return CIRC_ERR;
    };

    hid_t num_sum_terms_attr_id = H5Aopen(group_id, "num_sum_terms",
                                          H5P_DEFAULT);
    int num_sum_terms;
    status = H5Aread(num_sum_terms_attr_id, H5T_NATIVE_INT, &num_sum_terms);
    status = H5Aclose(num_sum_terms_attr_id);
    status = H5Aclose(num_qubits_attr_id);
    assert(status >= 0);

    log_info("num qubits: %d, num_sum_terms: %d", num_qubits,
             num_sum_terms);

    log_info("reading coeffs");
    hid_t dset_coeffs_id = H5Dopen2(group_id, "coeffs", H5P_DEFAULT);
    double *coeffs = (double *) malloc(sizeof(double) * num_sum_terms);
    assert(coeffs != NULL);
    status = H5Dread(dset_coeffs_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, coeffs);
    status = H5Dclose(dset_coeffs_id);

    log_info("computing hamiltonian one-norm");
    double hamil_one_norm = 0.0;
    for (int i = 0; i < num_sum_terms; i++) {
        hamil_one_norm += fabs(coeffs[i]);
    }
    log_info("norm: %f", hamil_one_norm);
    hid_t dspace_hamil_one_norm_id = H5Screate(H5S_SCALAR);
    hid_t attr_hamil_one_norm_id = H5Acreate1(group_id, "one_norm",
                                              H5T_IEEE_F64LE,
                                              dspace_hamil_one_norm_id,
                                              H5P_DEFAULT);
    status = H5Awrite(attr_hamil_one_norm_id, H5T_NATIVE_DOUBLE,
                      &hamil_one_norm);
    status = H5Sclose(dspace_hamil_one_norm_id);
    status = H5Aclose(attr_hamil_one_norm_id);

    log_info("reading paulis");
    hid_t dset_paulis_id = H5Dopen2(group_id, "paulis", H5P_DEFAULT);
    int *paulis = (int *) malloc(sizeof(int) * num_sum_terms * num_qubits);
    assert(paulis != NULL);
    status = H5Dread(dset_paulis_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, paulis);
    status = H5Dclose(dset_paulis_id);
    status = H5Gclose(group_id);

    log_debug("rescaling hamiltonian");
    for (int i = 0; i < num_sum_terms; i++) {
        coeffs[i] /= hamil_one_norm;
    }

    log_info("initialize Pauli Hamiltonian");
    PauliHamil hamil = createPauliHamil(num_qubits, num_sum_terms);
    initPauliHamil(hamil, coeffs, (enum pauliOpType *) paulis);
    reportPauliHamil(hamil);
    free(coeffs);
    free(paulis);
    circuit factory = rayon_circuit(&hamil);
    factory.num_sys_qb = hamil.numQubits;

    log_info("read time_series/times");
    hid_t time_series_id = H5Gopen1(file_id, "time_series");
    hid_t steps_attr_id = H5Aopen(time_series_id, "steps", H5P_DEFAULT);
    int steps_signed;
    status = H5Aread(steps_attr_id, H5T_NATIVE_INT, &steps_signed);
    size_t steps;
    assert(steps_signed > 0);
    steps = steps_signed;
    status = H5Aclose(steps_attr_id);
    log_info("number of steps: %d", steps);
    hid_t dset_times_id = H5Dopen1(time_series_id, "times");
    double *times = malloc(sizeof(double) * steps);
    assert(times != NULL);
    status = H5Dread(dset_times_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, times);
    status = H5Dclose(dset_times_id);

    double *values = malloc(sizeof(double) * steps);
    assert(values != NULL);

    rayon_circ_data circ_data = {
            .time = times[0],
            .imag_switch = 0,
    };
    factory.num_sys_qb = num_qubits;
    circ *circ = circ_create(factory, env, &circ_data);

    log_info("evaluating real part of expectation value");
    for (size_t i = 0; i < steps; i++) {
        circ_data.time = times[i];
        circ_reset(circ);
        circ_simulate(circ);

        double prob_0;
        int *mea_cl = circ_mea_cl(circ);
        double *mea_cl_prob = circ_mea_cl_prob(circ);
        if (mea_cl[0] == 0) {
            prob_0 = mea_cl_prob[0];
        } else {
            prob_0 = 1.0 - mea_cl_prob[0];
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

    log_info("evaluating imag part of expectation value");
    circ_data.imag_switch = 1;
    for (size_t i = 0; i < steps; i++) {
        circ_data.time = times[i];
        circ_reset(circ);
        circ_simulate(circ);

        double prob_0;
        int *mea_cl = circ_mea_cl(circ);
        double *mea_cl_prob = circ_mea_cl_prob(circ);
        if (mea_cl[0] == 0) {
            prob_0 = mea_cl_prob[0];
        } else {
            prob_0 = 1.0 - mea_cl_prob[0];
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

    return CIRC_OK;
}
