#include <stdio.h>
#include <stdlib.h>

#include "QuEST.h"

#include "circ.h"

static size_t ct_num_tot_qb(const struct circuit *ct)
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

int circ_init(struct circ *c, struct circ_env *env, struct circuit *ct,
	      void *data)
{
	Qureg *qureg = malloc(sizeof(Qureg));
	if (!qureg)
		goto qureg_alloc_fail;
	const QuESTEnv *quest_env = env->quest_env;
	*qureg = createQureg(ct_num_tot_qb(ct), *quest_env);

	int *mea_cl = malloc(sizeof(int) * ct->num_mea_qb);
	if (!mea_cl)
		goto mea_alloc_fail;
	int *qb = malloc(sizeof(int) * ct_num_tot_qb(ct));
	if (!qb)
		goto qb_alloc_fail;

	c->env = env;
	c->ct = ct;
	c->data = data;
	c->qb = qureg;
	c->mea_cl = mea_cl;
	c->mea_qb = qb;
	c->sys_qb = c->mea_qb + ct->num_mea_qb;
	c->anc_qb = c->sys_qb + ct->num_sys_qb;
	for (size_t i = 0; i < ct_num_tot_qb(c->ct); i++) {
		c->mea_qb[i] = i;
	}

	return 0;

qb_alloc_fail:
	free(mea_cl);
mea_alloc_fail:
	destroyQureg(*qureg, *quest_env);
	free(qureg);
qureg_alloc_fail:
	return -1;
}

void circ_destroy(struct circ *c)
{
	if (c->qb) {
		const QuESTEnv *quest_env = c->env->quest_env;
		const Qureg *qureg = c->qb;
		destroyQureg(*qureg, *quest_env);
		free(c->qb);
		c->qb = NULL;
	}
	if (c->mea_qb) {
		free(c->mea_qb);
		c->mea_qb = NULL;
		c->sys_qb = NULL;
		c->anc_qb = NULL;
	}
	if (c->mea_cl) {
		free(c->mea_cl);
		c->mea_cl = NULL;
	}
	c->data = NULL;
	c->ct = NULL;
	c->env = NULL;
}

void circ_report(struct circ const *c)
{
	const Qureg *qureg = c->qb;
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

	const circ_op ops[3] = {
		[0] = c->ct->prepst, [1] = c->ct->effect, [2] = c->ct->measure
	};
	for (int i = 0; i < 3; i++) {
		if (ops[i]) {
			if (ops[i](c) < 0)
				return -1;
		}
	}

	return 0;
}
