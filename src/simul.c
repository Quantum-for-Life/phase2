#include <stdlib.h>
#include <string.h>

#include "log/log.h"
#include "circ.h"


int simul_linen(circ_env *env);

int simul_rayon(circ_env *env, char *hamil_file);


void exit_failure() {
    log_error("Failed.");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

    log_info("Start");
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
        exit_failure();
    }

    if (strncmp(argv[1], "linen",5) == 0) {
        log_info("Circuit: linen");
        simul_linen(env);
    }

    if (strncmp(argv[1], "rayon",5) == 0) {
        log_info("Circuit: rayon");
        if (argc < 3) {
            log_error("No simul.h5 file provided");
            exit_failure();
        }
        simul_rayon(env, argv[2]);
    }

    log_info("Shut down simulation environment");
    circ_destroy_env(env);

    log_info("Done");
    fclose(log_file);
    return EXIT_SUCCESS;
}
