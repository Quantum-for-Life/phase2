#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "logger.h"
#include "circ.h"


int simul_linen(circ_env *env);

int simul_rayon(circ_env *env, char *hamil_file);


int main(int argc, char **argv) {

    logger("debug", "init");

    logger("info", "initialize circ_env");
    circ_env *env = circ_create_env();
    if (env == NULL) {
        logger("error", "cannot initialize circ_env");
        return EXIT_FAILURE;
    }
    logger("info", "circ_env initialized");

    logger("debug", "parsing command line arguments");
    if (argc < 2) {
        logger("error", "no circuit specified");
        return EXIT_FAILURE;
    }

    if (strcmp(argv[1], "linen") == 0) {
        logger("info", "running circuit: linen");
        simul_linen(env);
    }

    if (strcmp(argv[1], "rayon") == 0) {
        logger("info", "running circuit: rayon");
        if (argc < 3) {
            logger("error", "no simul.h5 file provided");
            return EXIT_FAILURE;
        }
        simul_rayon(env, argv[2]);
    }

    logger("info", "deactivate circ_env");
    circ_destroy_env(env);

    logger("info", "done");
    return EXIT_SUCCESS;
}
