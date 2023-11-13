#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"

#include "circ.h"
#include "rayon.h"
#include "linen.h"
#include "sdat.h"
#include "log/log.h"

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"
#define PHASE2_LOG_FILE "simul.log"

#define SIMUL_DEFAULT_H5FILE "simul.h5"

circ_result linen_simulate(circ_env *env, const char *filename);

circ_result rayon_simulate(circ_env *env, const char *filename);


void exit_failure() {
    log_error("Failed.");
    exit(EXIT_FAILURE);
}

void print_help_page(int argc, char **argv) {
    (void) argc;
    fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5]\n\n", argv[0]);
    fprintf(stderr,
            "where CIRCUIT is must be one of:\n"
            "\n"
            "    linen\n"
            "    rayon\n"
            "\n"
            "If no simulation input file (HDF5) is specified,\n"
            "the default is " SIMUL_DEFAULT_H5FILE " in the current directory.\n"
    );
    fprintf(stderr,
            "\n"
            "To set logging level, set environment variable:\n"
            "\n    "  PHASE2_LOG_ENVVAR
            "={trace, debug, info, warn, error, fatal}"
            "\n\n"
    );
}

void simul_set_log_level() {
    char *log_level = getenv(PHASE2_LOG_ENVVAR);
    if (log_level == NULL) {
        log_set_level(LOG_ERROR);
        return;
    }
    if (strncmp(log_level, "trace", 5) == 0) {
        log_set_level(LOG_TRACE);
    }
    if (strncmp(log_level, "debug", 5) == 0) {
        log_set_level(LOG_DEBUG);
    }
    if (strncmp(log_level, "info", 4) == 0) {
        log_set_level(LOG_INFO);
    }
    if (strncmp(log_level, "warn", 4) == 0) {
        log_set_level(LOG_WARN);
    }
    if (strncmp(log_level, "error", 5) == 0) {
        log_set_level(LOG_ERROR);
    }
    if (strncmp(log_level, "fatal", 5) == 0) {
        log_set_level(LOG_FATAL);
    }
}

int main(int argc, char **argv) {

    simul_set_log_level();
    log_info("*** Init ***");
    log_info("Open log file: " PHASE2_LOG_FILE);
    FILE *log_file = fopen(PHASE2_LOG_FILE, "a");
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
        print_help_page(argc, argv);
        return EXIT_FAILURE;
    }
    const char *h5filename = SIMUL_DEFAULT_H5FILE;
    if (argc < 3) {
        log_debug("No simulation input file specified; "
                  "using default: %s",
                  SIMUL_DEFAULT_H5FILE);
    } else {
        h5filename = argv[2];
    }

    log_info("*** Circuit ***");
    if (strncmp(argv[1], "linen", 5) == 0) {
        log_info("Circuit: linen");
        if (linen_simulate(env, h5filename) != CIRC_OK) {
            exit_failure();
        }
    } else if (strncmp(argv[1], "rayon", 5) == 0) {
        log_info("Circuit: rayon");
        if (rayon_simulate(env, h5filename) != CIRC_OK) {
            exit_failure();
        }
    } else {
        log_error("No circuit named %s", argv[1]);
        exit_failure();
    }

    log_info("*** Cleanup ***");
    log_info("Shut down simulation environment");
    circ_destroy_env(env);

    log_info("Done");
    fclose(log_file);
    return EXIT_SUCCESS;
}


circ_result
linen_simulate(circ_env *env, const char *h5filename) {
    (void) h5filename;

    log_debug("Read simulation input file: %s", h5filename);
    hid_t file_id = H5Fopen(h5filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    sdat_pauli_hamil ph;
    sdat_pauli_hamil_init(&ph);
    sdat_pauli_hamil_read(&ph, file_id);
    log_debug("Hamiltonian: num_qubits=%zu, num_sum_terms=%zu",
              ph.num_qubits, ph.num_sum_terms);
    sdat_pauli_hamil_drop(ph);

    sdat_time_series ts;
    sdat_time_series_init(&ts);
    sdat_time_series_read(&ts, file_id);
    log_debug("Time series: num_steps=%zu", ts.num_steps);

    H5Fclose(file_id);

    log_debug("Report simulation environment");
    circ_report_env(env);
    circuit factory = linen_circuit_factory(NULL);
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


circ_result rayon_simulate(circ_env *env, const char *h5filename) {

    log_debug("Read simulation input file: %s", h5filename);
    hid_t file_id = H5Fopen(h5filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    sdat_pauli_hamil ph;
    sdat_pauli_hamil_init(&ph);
    sdat_pauli_hamil_read(&ph, file_id);
    log_debug("Hamiltonian: num_qubits=%zu, num_sum_terms=%zu",
              ph.num_qubits, ph.num_sum_terms);

    sdat_time_series ts;
    sdat_time_series_init(&ts);
    sdat_time_series_read(&ts, file_id);
    log_debug("Time series: num_steps=%zu", ts.num_steps);

    H5Fclose(file_id);

    log_info("initialize Pauli Hamiltonian");
    for (size_t i = 0; i < ph.num_sum_terms; i++) {
        printf("%f ", ph.coeffs[i]);
        for (size_t j = 0; j < ph.num_qubits; j++) {
            printf("%u ", ph.paulis[i * ph.num_qubits + j]);
        }
        printf("\n");
    }
    PauliHamil hamil = createPauliHamil(ph.num_qubits, ph.num_sum_terms);

    int *paulis = malloc(sizeof(int) * ph.num_qubits * ph.num_sum_terms);
    if (paulis == NULL) {
        return CIRC_ERR;
    }
    for (size_t i = 0; i < ph.num_qubits * ph.num_sum_terms; i++) {
        paulis[i] = ph.paulis[i];
    }
    initPauliHamil(hamil, ph.coeffs, (enum pauliOpType *) paulis);
    free(paulis);

    reportPauliHamil(hamil);

    rayon_circuit_data ct_data = {.hamil = hamil, .data = NULL};
    circuit factory = rayon_circuit_factory(&ct_data);
    factory.num_sys_qb = hamil.numQubits;

    rayon_circ_data circ_data = {
            .time = ts.times[0],
            .imag_switch = 0,
    };
    factory.num_sys_qb = hamil.numQubits;
    circ *circ = circ_create(factory, env, &circ_data);

    log_info("evaluating real part of expectation value");
    for (size_t i = 0; i < ts.num_steps; i++) {
        circ_data.time = ts.times[i];

        double prob_0;
        int *mea_cl = circ_mea_cl(circ);
        double *mea_cl_prob = circ_mea_cl_prob(circ);

        circ_simulate(circ);
        if (mea_cl[0] == 0) {
            prob_0 = mea_cl_prob[0];
        } else {
            prob_0 = 1.0 - mea_cl_prob[0];
        }
        ts.values[2 * i] = 2 * prob_0 - 1;
        
    }

    log_info("evaluating imag part of expectation value");
    circ_data.imag_switch = 1;
    for (size_t i = 0; i < ts.num_steps; i++) {
        circ_data.time = ts.times[i];
        circ_simulate(circ);

        double prob_0;
        int *mea_cl = circ_mea_cl(circ);
        double *mea_cl_prob = circ_mea_cl_prob(circ);
        if (mea_cl[0] == 0) {
            prob_0 = mea_cl_prob[0];
        } else {
            prob_0 = 1.0 - mea_cl_prob[0];
        }
        ts.values[2 * i + 1] = 2 * prob_0 - 1;
    }

    file_id = H5Fopen(h5filename, H5F_ACC_RDWR, H5P_DEFAULT);
    sdat_time_series_write(ts, file_id);
    H5Fclose(file_id);

    sdat_pauli_hamil_drop(ph);
    sdat_time_series_drop(ts);

    return CIRC_OK;
}
