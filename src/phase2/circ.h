/** Circuit interface.
 */

#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

#include "QuEST.h"

enum {
        CIRC_OK,
        CIRC_ERR,
};


/** circ environment for a run on MPI cluster.
 *
 * To be created and destroyed _only once_
 * during the experiment.
 */

struct circ_env {
        QuESTEnv *quest_env;
};

struct circ;

/** circ specification.
 */
struct circuit {
        const char *name;
        void *data;

        size_t num_mea_qb;
        size_t num_sys_qb;
        size_t num_anc_qb;

        /** Reset the circ.
         *
         * On top of standard resetting the Qureg and mea_cl.
         * This can be NULL.
         */
        int (*reset)(struct circ *);

        /** circ specification divided into 4 steps.
         *
         * The order of execution of these steps in circ_simulate() is always:
         *   1. state_prep()
         *   2. routine()
         *   3. state_post()
         *
         * The cicuit is _not_ reset before the call to circ_simulate().
         *
         * The data pointer passed to these function is the same stored in
         * circ.data, and the same as passed to circ_create_circuit().
         *
         * There pointers can be NULL, meaning the step is not specified.
         */
        int (*state_prep)(struct circ *, void *);

        int (*routine)(struct circ *, void *);

        int (*state_post)(struct circ *, void *);
};


/** circ instance.
 */
struct circ {
        struct circ_env env;
        struct circuit ct;
        void *data;

        Qureg *qureg;

        int *mea_cl;
        double *mea_cl_prob;
        int *mea_qb;
        int *sys_qb;
        int *anc_qb;
};


int circ_env_init(struct circ_env *env);

void circ_env_destroy(struct circ_env *env);

void circ_env_report(struct circ_env const *env);

int circ_init(struct circ *c,
              struct circ_env env, struct circuit ct, void *data);

void circ_destroy(struct circ *c);

void circ_report(struct circ const *c);

int circ_reset(struct circ *c);

int circ_simulate(struct circ *c);

#endif //PHASE2_CIRC_H