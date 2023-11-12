#include <stdlib.h>
#include <string.h>
#include "hdf5.h"

#include "log/log.h"
#include "circ.h"

#define SIMUL_DEFAULT_H5FILE "simul.h5"

int linen_simulate(circ_env *env, circuit_data data);

int rayon_simulate(circ_env *env, char *hamil_file);


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
    circuit_data ct_data = {.data = NULL, .hamil = createPauliHamil(1,1)};
    if (strncmp(argv[1], "linen", 5) == 0) {
        log_info("Circuit: linen");
        if (linen_simulate(env, ct_data) != 0) {
            exit_failure();
        }
    }

    if (strncmp(argv[1], "rayon", 5) == 0) {
        log_info("Circuit: rayon");
        if (rayon_simulate(env, argv[2]) != 0) {
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
