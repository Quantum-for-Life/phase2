#include <stdio.h>
#include <stdlib.h>

#include "QuEST.h"

#include "circ.h"

size_t circ_circuit_num_tot_qb(const struct circuit ct) {
        return ct.num_mea_qb + ct.num_sys_qb + ct.num_anc_qb;
}

int circ_env_init(struct circ_env* env) {
        env->quest_env = malloc(sizeof(QuESTEnv));
        if (!env->quest_env) {
                return CIRC_ERR;
        }
        *env->quest_env = createQuESTEnv();

        return CIRC_OK;
}

void circ_env_destroy(struct circ_env* env) {
        if (!env->quest_env) {
                return;
        }
        destroyQuESTEnv(*env->quest_env);
        free(env->quest_env);
        env->quest_env = NULL;
}

void circ_env_report(struct circ_env const* env) {
        reportQuESTEnv(*env->quest_env);
}

void zero_mea_cl(const struct circ* c) {
        for (size_t i = 0; i < c->ct.num_mea_qb; i++) {
                c->mea_cl[i] = 0;
        }
}

void init_mea_qb(const struct circ* c) {
        for (size_t i = 0; i < c->ct.num_mea_qb; i++) {
                c->mea_qb[i] = i;
        }
}

void init_sys_qb(const struct circ* c) {
        for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
                c->sys_qb[i] = c->ct.num_mea_qb + i;
        }
}

void init_anc_qb(const struct circ* c) {
        for (size_t i = 0; i < c->ct.num_anc_qb; i++) {
                c->anc_qb[i] = c->ct.num_mea_qb + c->ct.num_sys_qb + i;
        }
}

int
circ_init(struct circ* c,
          const struct circ_env env, struct circuit const ct, void* data) {
        c->env = env;
        c->ct = ct;
        c->data = data;

        Qureg* qureg = malloc(sizeof(Qureg));
        if (!qureg) {
                return CIRC_ERR;
        }
        *qureg = createQureg(
                circ_circuit_num_tot_qb(ct),
                *env.quest_env);
        c->qureg = qureg;

        int* mea_cl = malloc(sizeof(int) * ct.num_mea_qb);
        double* mea_cl_prob = malloc(sizeof(double) * ct.num_mea_qb);
        int* mea_qb = malloc(sizeof(int) * ct.num_mea_qb);
        int* sys_qb = malloc(sizeof(int) * ct.num_sys_qb);
        int* anc_qb = malloc(sizeof(int) * ct.num_anc_qb);
        if (!(mea_cl && mea_cl_prob && mea_qb && sys_qb && anc_qb)) {
                free(anc_qb);
                free(sys_qb);
                free(mea_qb);
                free(mea_cl_prob);
                free(mea_cl);
                return CIRC_ERR;
        }
        c->mea_cl = mea_cl;
        c->mea_cl_prob = mea_cl_prob;
        c->mea_qb = mea_qb;
        c->sys_qb = sys_qb;
        c->anc_qb = anc_qb;

        zero_mea_cl(c);
        init_mea_qb(c);
        init_sys_qb(c);
        init_anc_qb(c);

        return CIRC_OK;
}

void circ_destroy(struct circ* c) {
        destroyQureg(*c->qureg, *c->env.quest_env);
        free(c->qureg);
        c->qureg = NULL;
        free(c->mea_qb);
        c->mea_qb = NULL;
        free(c->sys_qb);
        c->sys_qb = NULL;
        free(c->anc_qb);
        c->anc_qb = NULL;
        free(c->mea_cl_prob);
        c->mea_cl_prob = NULL;
        free(c->mea_cl);
        c->mea_cl = NULL;
}

void circ_report(struct circ const* c) {
        printf("----------------\n");
        printf("CIRCUIT: %s\n", c->ct.name);
        reportQuregParams(*c->qureg);

        printf("mea_cl register: [");
        for (size_t i = 0; i < c->ct.num_mea_qb; i++) {
                printf("%d", c->mea_cl[i]);
        }
        printf("]\n");

        printf("mea_qb indices: { ");
        for (size_t i = 0; i < c->ct.num_mea_qb; i++) {
                printf("%d ", c->mea_qb[i]);
        }
        printf("}\n");

        printf("sys_qb indices: { ");
        for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
                printf("%d ", c->sys_qb[i]);
        }
        printf("}\n");

        printf("anc_qb indices: { ");
        for (size_t i = 0; i < c->ct.num_anc_qb; i++) {
                printf("%d ", c->anc_qb[i]);
        }
        printf("}\n");
        printf("----------------\n");
}

int circ_reset(struct circ* c) {
        initZeroState(*c->qureg);
        zero_mea_cl(c);
        if (c->ct.reset) {
                return c->ct.reset(c);
        }
        return CIRC_OK;
}

int circ_simulate(struct circ* c) {
        circ_reset(c);

        int result;
        if (c->ct.state_prep) {
                result = c->ct.state_prep(c);
                if (result != CIRC_OK) {
                        return result;
                }
        }
        if (c->ct.routine) {
                result = c->ct.routine(c);
                if (result != CIRC_OK) {
                        return result;
                }
        }
        if (c->ct.state_post) {
                result = c->ct.state_post(c);
                if (result != CIRC_OK) {
                        return result;
                }
        }

        /* Measure qubits */
        for (size_t i = 0; i < c->ct.num_mea_qb; i++) {
                c->mea_cl[i] = measureWithStats(*c->qureg, c->mea_qb[i],
                                                &c->mea_cl_prob[i]);
        }

        return CIRC_OK;
}
