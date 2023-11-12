#include <stdlib.h>
#include <string.h>

#include "log/log.h"
#include "circ.h"


int simul_linen(circ_env *env);

int simul_rayon(circ_env *env, char *hamil_file);


int main(int argc, char **argv) {

    FILE *f_log = fopen("simul.log", "a");
    log_add_fp(f_log, LOG_DEBUG);
    log_debug("init");

    log_info("initialize circ_env");
    circ_env *env = circ_create_env();
    if (env == NULL) {
        log_error("cannot initialize circ_env");
        return EXIT_FAILURE;
    }
    log_info("circ_env initialized");

    log_debug("parsing command line arguments");
    if (argc < 2) {
        log_error("no circuit specified");
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "linen") == 0) {
        log_info("running circuit: linen");
        simul_linen(env);
    }

    if (strcmp(argv[1], "rayon") == 0) {
        log_info("running circuit: rayon");
        if (argc < 3) {
            log_error("no simul.h5 file provided");
            return EXIT_FAILURE;
        }
        simul_rayon(env, argv[2]);
    }

    log_info("deactivate circ_env");
    circ_destroy_env(env);

    log_info("done");
    fclose(f_log);
    return EXIT_SUCCESS;
}
