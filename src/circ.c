#include <stdio.h>
#include <stdlib.h>

#include "QuEST.h"

#include "circ.h"
#include "log/log.h"

#define CIRC_LOG_TAG "[circ] "

struct circ_env_ {
    QuESTEnv quest_env;
};

struct circ_ {
    circuit ct;
    /** Arbitrary data passed during initialization.
     *
     * The user is responsible to manage the memory pointed to.
     */
    void *data;

    circ_env *env;
    Qureg qureg;

    /** Subregisters:
     *
     *    mea (measurement)
     *       mea_cl - classical register
     *       mea_qb - quantum register
     *    sys (system)
     *    anc (ancilla)
     */
    int *mea_cl;
    double *mea_cl_prob;

    int *mea_qb;
    int *sys_qb;
    int *anc_qb;

    size_t simul_counter;
};


size_t circ_circuit_num_tot_qb(circuit ct) {
    return ct.num_mea_qb + ct.num_sys_qb + ct.num_anc_qb;
}

circ_env *circ_create_env() {
    log_debug(CIRC_LOG_TAG "Create circ_env");
    QuESTEnv quest_env = createQuESTEnv();

    circ_env *env = (circ_env *) malloc(sizeof(circ_env));
    if (env == NULL) {
        destroyQuESTEnv(quest_env);
        return NULL;
    }
    env->quest_env = quest_env;

    return env;
}

void circ_destroy_env(circ_env *const env) {
    log_debug(CIRC_LOG_TAG "Destroy circ_env");
    if (env == NULL) {
        return;
    }
    destroyQuESTEnv(env->quest_env);
    free(env);
}

void circ_report_env(circ_env *const env) {
    reportQuESTEnv(env->quest_env);
}

void zero_mea_cl(circ *const c) {
    for (size_t i = 0; i < c->ct.num_mea_qb; i++) {
        c->mea_cl[i] = 0;
    }
}

void init_mea_qb(circ *const c) {
    for (size_t i = 0; i < c->ct.num_mea_qb; i++) {
        c->mea_qb[i] = i;
    }
}

void init_sys_qb(circ *const c) {
    for (size_t i = 0; i < c->ct.num_sys_qb; i++) {
        c->sys_qb[i] = c->ct.num_mea_qb + i;
    }
}

void init_anc_qb(circ *const c) {
    for (size_t i = 0; i < c->ct.num_anc_qb; i++) {
        c->anc_qb[i] = c->ct.num_mea_qb + c->ct.num_sys_qb + i;
    }
}


circ *circ_create(circuit const ct, circ_env *const env, void *data) {
    log_debug(CIRC_LOG_TAG "Create circ");
    if (env == NULL) {
        return NULL;
    }
    circ *c = malloc(sizeof(circ));
    int *mea_cl = malloc(sizeof(int) * ct.num_mea_qb);
    double *mea_cl_prob = malloc(sizeof(double) * ct.num_mea_qb);
    int *mea_qb = malloc(sizeof(int) * circ_circuit_num_tot_qb(ct));
    if (!(c && mea_cl && mea_cl_prob && mea_qb)) {
        free(mea_qb);
        free(mea_cl_prob);
        free(mea_cl);
        free(c);
        return NULL;
    }
    c->mea_cl = mea_cl;
    c->mea_cl_prob = mea_cl_prob;
    c->mea_qb = mea_qb;
    c->sys_qb = c->mea_qb + ct.num_mea_qb;
    c->anc_qb = c->sys_qb + ct.num_sys_qb;

    c->ct = ct;
    c->data = data;
    c->env = env;

    c->qureg = createQureg(
            circ_circuit_num_tot_qb(ct),
            env->quest_env);
    c->simul_counter = 0;

    zero_mea_cl(c);
    init_mea_qb(c);
    init_sys_qb(c);
    init_anc_qb(c);

    return c;
}

void circ_destroy(circ *c) {
    log_debug(CIRC_LOG_TAG "Destroy circ");
    if (c == NULL) {
        return;
    }
    destroyQureg(c->qureg, c->env->quest_env);
    free(c->mea_qb);
    free(c->mea_cl_prob);
    free(c->mea_cl);
    free(c);
}

int *circ_mea_cl(circ *const c) {
    return c->mea_cl;
}

double *circ_mea_cl_prob(circ *const c) {
    return c->mea_cl_prob;
}

int *circ_mea_qb(circ *const c) {
    return c->mea_qb;
}

int *circ_sys_qb(circ *const c) {
    return c->sys_qb;
}

int *circ_anc_qb(circ *const c) {
    return c->anc_qb;
}

size_t circ_num_mea_qb(circ *const c) {
    return c->ct.num_mea_qb;
}

size_t circ_num_sys_qb(circ *const c) {
    return c->ct.num_sys_qb;
}

size_t circ_num_anc_qb(circ *const c) {
    return c->ct.num_anc_qb;
}


size_t circ_num_tot_qb(circ *const c) {
    return circ_circuit_num_tot_qb(c->ct);
}

const char *circ_name(circ *const c) {
    return c->ct.name;
}

void circ_report(circ *const c) {
    printf("----------------\n");
    printf("CIRCUIT: %s\n", circ_name(c));
    reportQuregParams(c->qureg);

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

void *circ_circuit_data(circ *const c) {
    return c->ct.data;
}

Qureg circ_qureg(circ *const c) {
    return c->qureg;
}

circ_result circ_reset(circ *const c) {
    log_trace(CIRC_LOG_TAG "reset");
    initZeroState(circ_qureg(c));
    zero_mea_cl(c);
    if (c->ct.reset != NULL) {
        return c->ct.reset(c);
    }
    return CIRC_OK;
}

circ_result circ_simulate(circ *const c) {
    if (c == NULL) {
        return CIRC_ERR;
    }
    c->simul_counter++;
    log_trace(CIRC_LOG_TAG "simulate (%zu)", c->simul_counter);

    circ_reset(c);

    circ_result result;
    if (c->ct.state_prep != NULL) {
        log_trace(CIRC_LOG_TAG "state_prep");
        result = c->ct.state_prep(c, c->data);
        if (result != CIRC_OK) {
            return result;
        }
    }
    if (c->ct.routine != NULL) {
        log_trace(CIRC_LOG_TAG "routine");
        result = c->ct.routine(c, c->data);
        if (result != CIRC_OK) {
            return result;
        }
    }
    if (c->ct.state_post != NULL) {
        log_trace(CIRC_LOG_TAG "state_post");
        result = c->ct.state_post(c, c->data);
        if (result != CIRC_OK) {
            return result;
        }
    }

    /* Measure qubits */
    log_trace(CIRC_LOG_TAG "measure");
    int *mea_cl = circ_mea_cl(c);
    double *mea_cl_prob = circ_mea_cl_prob(c);
    int *mea_qb = circ_mea_qb(c);
    Qureg qureg = circ_qureg(c);
    for (size_t i = 0; i < circ_num_mea_qb(c); i++) {
        mea_cl[i] = measureWithStats(qureg, mea_qb[i], &mea_cl_prob[i]);
    }

    return CIRC_OK;
}
