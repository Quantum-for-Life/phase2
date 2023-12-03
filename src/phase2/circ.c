#include <stdio.h>
#include <stdlib.h>

#include "QuEST.h"

#include "circ.h"

static size_t circ_circuit_num_tot_qb(const struct circuit *ct)
{
	return ct->num_mea_qb + ct->num_sys_qb + ct->num_anc_qb;
}

int circ_env_init(struct circ_env *env)
{
	QuESTEnv *quest_env = malloc(sizeof(QuESTEnv));
	if (!quest_env)
		return -1;
	*quest_env = createQuESTEnv();
	env->quest_env = quest_env;

	return 0;
}

void circ_env_destroy(struct circ_env *env)
{
	if (!env->quest_env)
		return;

	const QuESTEnv *quest_env = env->quest_env;
	destroyQuESTEnv(*quest_env);
	free(env->quest_env);
	env->quest_env = NULL;
}

void circ_env_report(struct circ_env const *env)
{
	const QuESTEnv *quest_env = env->quest_env;
	reportQuESTEnv(*quest_env);
}

static void init_circ_indices(struct circ *c)
{
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		c->mea_cl[i] = 0;
		c->mea_qb[i] = i;
	}
	for (size_t i = 0; i < c->ct->num_sys_qb; i++) {
		c->sys_qb[i] = c->ct->num_mea_qb + i;
	}
	for (size_t i = 0; i < c->ct->num_anc_qb; i++) {
		c->anc_qb[i] = c->ct->num_mea_qb + c->ct->num_sys_qb + i;
	}
}

int circ_init(struct circ *c, struct circ_env *env, struct circuit *ct,
	      void *data)
{
	const QuESTEnv *quest_env = env->quest_env;
	c->env = env;
	c->ct = ct;
	c->data = data;

	Qureg *qureg = malloc(sizeof(Qureg));
	if (!qureg)
		return -1;

	*qureg = createQureg(circ_circuit_num_tot_qb(ct), *quest_env);
	c->qureg = qureg;

	int *mea_cl = malloc(sizeof(int) * ct->num_mea_qb);
	double *mea_cl_prob = malloc(sizeof(double) * ct->num_mea_qb);
	int *mea_qb = malloc(sizeof(int) * ct->num_mea_qb);
	int *sys_qb = malloc(sizeof(int) * ct->num_sys_qb);
	int *anc_qb = malloc(sizeof(int) * ct->num_anc_qb);
	if (!(mea_cl && mea_cl_prob && mea_qb && sys_qb && anc_qb)) {
		free(anc_qb);
		free(sys_qb);
		free(mea_qb);
		free(mea_cl_prob);
		free(mea_cl);
		return -1;
	}
	c->mea_cl = mea_cl;
	c->mea_cl_prob = mea_cl_prob;
	c->mea_qb = mea_qb;
	c->sys_qb = sys_qb;
	c->anc_qb = anc_qb;

	init_circ_indices(c);

	return 0;
}

void circ_destroy(struct circ *c)
{
	const QuESTEnv *quest_env = c->env->quest_env;
	const Qureg *qureg = c->qureg;
	destroyQureg(*qureg, *quest_env);
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

void circ_report(struct circ const *c)
{
	const Qureg *qureg = c->qureg;
	printf("----------------\n");
	printf("CIRCUIT: %s\n", c->ct->name);
	reportQuregParams(*qureg);

	printf("mea_cl register: [");
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		printf("%d", c->mea_cl[i]);
	}
	printf("]\n");

	printf("mea_qb indices: { ");
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		printf("%d ", c->mea_qb[i]);
	}
	printf("}\n");

	printf("sys_qb indices: { ");
	for (size_t i = 0; i < c->ct->num_sys_qb; i++) {
		printf("%d ", c->sys_qb[i]);
	}
	printf("}\n");

	printf("anc_qb indices: { ");
	for (size_t i = 0; i < c->ct->num_anc_qb; i++) {
		printf("%d ", c->anc_qb[i]);
	}
	printf("}\n");
	printf("----------------\n");
}

int circ_reset(struct circ *c)
{
	const Qureg *qureg = c->qureg;
	initZeroState(*qureg);
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		c->mea_cl[i] = 0;
	}
	if (c->ct->reset)
		return c->ct->reset(c);

	return 0;
}

int circ_simulate(struct circ *c)
{
	if (circ_reset(c) < 0)
		return -1;

	const struct circuit *ct = c->ct;
	int (*ops[3])(struct circ *) = {
		[0] = ct->state_prep, [1] = ct->routine, [2] = ct->state_post
	};
	for (int i = 0; i < 3; i++) {
		if (ops[i]) {
			if (ops[i](c) < 0)
				return -1;
		}
	}

	/* Measure qubits */
	const Qureg *qureg = c->qureg;
	for (size_t i = 0; i < ct->num_mea_qb; i++) {
		c->mea_cl[i] = measureWithStats(*qureg, c->mea_qb[i],
						&c->mea_cl_prob[i]);
	}

	return 0;
}
