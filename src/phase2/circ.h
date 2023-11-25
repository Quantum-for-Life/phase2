#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdlib.h>

#include "QuEST.h"

enum {
        CIRC_OK,
        CIRC_ERR,
};

struct circ_env {
        QuESTEnv *quest_env;
};

struct circ;

struct circuit {
        const char *name;
        void *data;

        size_t num_mea_qb;
        size_t num_sys_qb;
        size_t num_anc_qb;

        int (*reset)(struct circ *);

        int (*state_prep)(struct circ *);

        int (*routine)(struct circ *);

        int (*state_post)(struct circ *);
};

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