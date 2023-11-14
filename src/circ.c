#include <stdio.h>
#include <stdlib.h>

#include "QuEST.h"

#include "circ.h"
#include "log/log.h"

#define CIRC_LOG_TAG "[circ] "

size_t circ_circuit_num_tot_qb(circuit ct) {
    return ct.num_mea_qb + ct.num_sys_qb + ct.num_anc_qb;
}

circ_result circ_env_init(circ_env *env) {
    log_debug(CIRC_LOG_TAG "Init circ_env");
    QuESTEnv quest_env = createQuESTEnv();
    env->quest_env = quest_env;

    return CIRC_OK;
}

void circ_env_drop(circ_env env) {
    log_debug(CIRC_LOG_TAG "Destroy circ_env");
    destroyQuESTEnv(env.quest_env);
}

void circ_env_report(circ_env env) {
    reportQuESTEnv(env.quest_env);
}

void zero_mea_cl(circ c) {
    for (size_t i = 0; i < c.ct.num_mea_qb; i++) {
        c.mea_cl[i] = 0;
    }
}

void init_mea_qb(circ c) {
    for (size_t i = 0; i < c.ct.num_mea_qb; i++) {
        c.mea_qb[i] = i;
    }
}

void init_sys_qb(circ c) {
    for (size_t i = 0; i < c.ct.num_sys_qb; i++) {
        c.sys_qb[i] = c.ct.num_mea_qb + i;
    }
}

void init_anc_qb(circ c) {
    for (size_t i = 0; i < c.ct.num_anc_qb; i++) {
        c.anc_qb[i] = c.ct.num_mea_qb + c.ct.num_sys_qb + i;
    }
}

circ_result
circ_init(circ *c, circuit const ct, circ_env env, void *data) {
    log_debug(CIRC_LOG_TAG "Init circ");
    log_debug("num_tot: %zu", circ_circuit_num_tot_qb(ct));

    int *mea_cl = malloc(sizeof(int) * ct.num_mea_qb);
    double *mea_cl_prob = malloc(sizeof(double) * ct.num_mea_qb);
    int *mea_qb = malloc(sizeof(int) * circ_circuit_num_tot_qb(ct));
    if (!(mea_cl && mea_cl_prob && mea_qb)) {
        free(mea_qb);
        free(mea_cl_prob);
        free(mea_cl);
        return CIRC_ERR;
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
            (int) circ_circuit_num_tot_qb(ct),
            env.quest_env);
    c->simul_counter = 0;

    zero_mea_cl(*c);
    init_mea_qb(*c);
    init_sys_qb(*c);
    init_anc_qb(*c);

    return CIRC_OK;
}

void circ_drop(circ c) {
    log_debug(CIRC_LOG_TAG "Destroy circ");
    destroyQureg(c.qureg, c.env.quest_env);
    free(c.mea_qb);
    free(c.mea_cl_prob);
    free(c.mea_cl);
}

size_t circ_num_tot_qb(circ c) {
    return circ_circuit_num_tot_qb(c.ct);
}

const char *circ_name(circ c) {
    return c.ct.name;
}

void circ_report(circ c) {
    printf("----------------\n");
    printf("CIRCUIT: %s\n", circ_name(c));
    reportQuregParams(c.qureg);

    printf("mea_cl register: [");
    for (size_t i = 0; i < c.ct.num_mea_qb; i++) {
        printf("%d", c.mea_cl[i]);
    }
    printf("]\n");

    printf("mea_qb indices: { ");
    for (size_t i = 0; i < c.ct.num_mea_qb; i++) {
        printf("%d ", c.mea_qb[i]);
    }
    printf("}\n");

    printf("sys_qb indices: { ");
    for (size_t i = 0; i < c.ct.num_sys_qb; i++) {
        printf("%d ", c.sys_qb[i]);
    }
    printf("}\n");

    printf("anc_qb indices: { ");
    for (size_t i = 0; i < c.ct.num_anc_qb; i++) {
        printf("%d ", c.anc_qb[i]);
    }
    printf("}\n");
    printf("----------------\n");
}

void *circ_circuit_data(circ c) {
    return c.ct.data;
}


circ_result circ_reset(circ c) {
    log_trace(CIRC_LOG_TAG "reset");
    initZeroState(c.qureg);
    zero_mea_cl(c);
    if (c.ct.reset != NULL) {
        return c.ct.reset(c);
    }
    return CIRC_OK;
}

circ_result circ_simulate(circ c) {
    c.simul_counter++;
    log_trace(CIRC_LOG_TAG "simulate (%zu)", c.simul_counter);

    circ_reset(c);

    circ_result result;
    if (c.ct.state_prep != NULL) {
        log_trace(CIRC_LOG_TAG "state_prep");
        result = c.ct.state_prep(c, c.data);
        if (result != CIRC_OK) {
            return result;
        }
    }
    if (c.ct.routine != NULL) {
        log_trace(CIRC_LOG_TAG "routine");
        result = c.ct.routine(c, c.data);
        if (result != CIRC_OK) {
            return result;
        }
    }
    if (c.ct.state_post != NULL) {
        log_trace(CIRC_LOG_TAG "state_post");
        result = c.ct.state_post(c, c.data);
        if (result != CIRC_OK) {
            return result;
        }
    }

    /* Measure qubits */
    log_trace(CIRC_LOG_TAG "measure");
    for (size_t i = 0; i < c.ct.num_mea_qb; i++) {
        c.mea_cl[i] = measureWithStats(c.qureg, c.mea_qb[i],
                                       &c.mea_cl_prob[i]);
    }

    return CIRC_OK;
}
