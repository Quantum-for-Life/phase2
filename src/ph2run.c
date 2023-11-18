#include <stdlib.h>
#include <string.h>

#ifdef DISTRIBUTED

#include "mpi.h"

#endif

#include "hdf5.h"

#include "circ.h"
#include "rayon.h"
#include "linen.h"
#include "sdat.h"
#include "log/log.h"

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"
#define PHASE2_LOG_FILE "simul.log"

#define PHASE2_DEFAULT_H5FILE "simul.h5"

int linen_simulate(struct circ_env env);

int rayon_simulate(struct circ_env env, sdat_pauli_hamil ph,
                   sdat_time_series *ts);

void exit_failure(const char *msg) {
    log_error("Failure: %s", msg);
    exit(EXIT_FAILURE);
}

void help_page(int argc, char **argv) {
    (void) argc;
    fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5]\n\n", argv[0]);
    fprintf(stderr,
            "where CIRCUIT is must be one of:\n"
            "\n"
            "    linen\n"
            "    rayon\n"
            "\n"
            "If no simulation input file (HDF5) is specified,\n"
            "the default is " PHASE2_DEFAULT_H5FILE " in the current directory.\n"
    );
    fprintf(stderr,
            "\n"
            "To set logging level, set environment variable:\n"
            "\n    "  PHASE2_LOG_ENVVAR
            "={trace, debug, info, warn, error, fatal}"
            "\n\n"
    );
}

void set_log_level() {
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

    struct circ_env env;
    if (circ_env_init(&env) != circ_ok) {
        exit_failure("initialize environment");
    }

#ifdef DISTRIBUTED

    int rank, num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    log_info("*** Init ***");
    log_info("Initialize MPI environment");
    log_info("MPI num_ranks: %d", num_ranks);
    log_info("This is rank no. %d", rank);
#else
    log_info("*** Init ***");
    log_info("MPI mode not enabled.");
    log_info("To enable distributed mode, set "
             "-DDISTRIBUTED "
             "flag during compilation");
#endif

    set_log_level();
    log_info("Open log file: " PHASE2_LOG_FILE);
    FILE *log_file = fopen(PHASE2_LOG_FILE, "a");
    if (!log_file) {
        exit_failure("open log file");
    }
    log_add_fp(log_file, LOG_DEBUG);

    log_debug("Parsing command line arguments");
    if (argc < 2) {
        help_page(argc, argv);
        return EXIT_FAILURE;
    }
    const char *h5filename = PHASE2_DEFAULT_H5FILE;
    if (argc < 3) {
        log_debug("No simulation input file specified; "
                  "using default: %s",
                  PHASE2_DEFAULT_H5FILE);
    } else {
        h5filename = argv[2];
    }

    log_debug("Read simulation input file: %s", h5filename);
    hid_t file_id, access_plist;
#ifdef DISTRIBUTED
    access_plist = H5Pcreate(H5P_FILE_ACCESS);
    log_debug("H5P_FILE_ACCESS: %" PRId64, access_plist);
    H5Pset_fapl_mpio(access_plist, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    access_plist = H5P_DEFAULT;
#endif
    file_id = H5Fopen(h5filename, H5F_ACC_RDONLY, access_plist);

    sdat_pauli_hamil dat_ph;
    sdat_pauli_hamil_init(&dat_ph);
    if (sdat_pauli_hamil_read(&dat_ph, file_id) != sdat_ok) {
        exit_failure("read Hamiltonian data");
    }
    log_debug("Hamiltonian: num_qubits=%zu, num_sum_terms=%zu",
              dat_ph.num_qubits, dat_ph.num_sum_terms);

    sdat_time_series dat_ts;
    sdat_time_series_init(&dat_ts);
    if (sdat_time_series_read(&dat_ts, file_id) != sdat_ok) {
        exit_failure("read time series data");
    }
    log_debug("Time series: num_steps=%zu", dat_ts.num_steps);
    H5Fclose(file_id);

    log_info("*** Circuit ***");
    int sucess;
    if (strncmp(argv[1], "linen", 5) == 0) {
        log_info("Circuit: linen");
        sucess = linen_simulate(env) == circ_ok;

    } else if (strncmp(argv[1], "rayon", 5) == 0) {
        log_info("Circuit: rayon");
        sucess = rayon_simulate(env, dat_ph, &dat_ts) != circ_ok;
    } else {
        log_error("No circuit named %s", argv[1]);
        sucess = 0;
    }
    if (!sucess) {
        exit_failure("simulation error");
    }

    log_info("Saving data");
#ifdef DISTRIBUTED
    access_plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access_plist, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    access_plist = H5P_DEFAULT;
#endif
    file_id = H5Fopen(h5filename, H5F_ACC_RDWR, access_plist);
    sdat_time_series_write(dat_ts, file_id);
    H5Fclose(file_id);

    log_info("*** Cleanup ***");
    sdat_pauli_hamil_drop(dat_ph);
    sdat_time_series_drop(dat_ts);

    log_info("Shut down simulation environment");
    circ_env_destroy(&env);

    log_info("Done");
    fclose(log_file);
    return EXIT_SUCCESS;
}


int
linen_simulate(struct circ_env env) {

    log_debug("Report simulation environment");
    circ_env_report(env);

    struct circuit factory = linen_circuit_factory(NULL);
    struct circ c;
    if (circ_init(&c, factory, env, NULL) != circ_ok) {
        log_error("Cannot initialize circuit");
        return circ_err;
    }
    log_debug("\"linen\" circuit created");
    circ_report(c);
    log_debug("Simulating circuit");
    circ_simulate(&c);
    log_debug("Free circuit instance");
    circ_destroy(&c);

    return circ_ok;
}


void
rayon_simulate_cleanup(PauliHamil hamil) {
    destroyPauliHamil(hamil);
}

int
rayon_simulate(struct circ_env env, sdat_pauli_hamil ph,
               sdat_time_series *ts) {
    log_info("Initialize Pauli Hamiltonian");
    PauliHamil hamil = createPauliHamil(ph.num_qubits, ph.num_sum_terms);
    for (size_t i = 0; i < ph.num_sum_terms; i++) {
        hamil.termCoeffs[i] = ph.coeffs[i];
        for (size_t j = 0; j < ph.num_qubits; j++) {
            size_t pauli_idx = i * ph.num_qubits + j;
            hamil.pauliCodes[pauli_idx] = (enum pauliOpType) ph.paulis[pauli_idx];
        }
    }

    log_info("Initialize circuit");
    struct rayon_circuit_data ct_data = {.hamil = hamil, .data = NULL};
    struct circuit factory = rayon_circuit_factory(&ct_data);
    struct rayon_circ_data circ_data = {.imag_switch = 0, .time = 0.0};
    struct circ c;
    if (circ_init(&c, factory, env, &circ_data) != circ_ok) {
        log_error("Cannot initialize circuit");
        return circ_err;
    }

    log_info("Computing expectation value");
    double prob_0;
    circ_data.imag_switch = 0;
    while (circ_data.imag_switch <= 1) {
        size_t offset = circ_data.imag_switch == 0 ? 0 : 1;
        for (size_t i = 0; i < ts->num_steps; i++) {
            circ_data.time = ts->times[i];
            if (circ_simulate(&c) != circ_ok) {
                log_error("Simulation error");
                rayon_simulate_cleanup(hamil);
                return circ_err;
            }
            prob_0 = c.mea_cl[0] == 0 ? c.mea_cl_prob[0] : 1.0 -
                                                           c.mea_cl_prob[0];
            ts->values[2 * i + offset] = 2 * prob_0 - 1;
        }
        circ_data.imag_switch++;
    }
    circ_destroy(&c);
    log_info("End of simulation");
    rayon_simulate_cleanup(hamil);

    return circ_ok;
}
