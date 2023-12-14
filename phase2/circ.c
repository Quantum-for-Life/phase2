#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "QuEST.h"

#include "circ.h"

static struct {
	_Bool init;
	QuESTEnv quest_env;
} CIRC_ENV;

static void env_init()
{
	CIRC_ENV.quest_env = createQuESTEnv();
	CIRC_ENV.init = true;
}

static void env_destroy()
{
	if (CIRC_ENV.init) {
		destroyQuESTEnv(CIRC_ENV.quest_env);
		CIRC_ENV.init = false;
	}
}

static QuESTEnv env_get_questenv(void)
{
	if (!CIRC_ENV.init) {
		env_init(&CIRC_ENV);
	}
	return CIRC_ENV.quest_env;
}

static void env_report(void)
{
	reportQuESTEnv(env_get_questenv());
}

static size_t ct_num_tot_qb(const struct circuit *ct)
{
	return ct->num_mea_qb + ct->num_sys_qb + ct->num_anc_qb;
}

int circ_init(struct circ *c, struct circuit *ct, void *data)
{
	Qureg *qureg = malloc(sizeof(Qureg));
	if (!qureg)
		goto qureg_alloc_fail;

	*qureg = createQureg(ct_num_tot_qb(ct), env_get_questenv());

	int *mea_cl = malloc(sizeof(int) * ct->num_mea_qb);
	if (!mea_cl)
		goto mea_alloc_fail;
	int *qb = malloc(sizeof(int) * ct_num_tot_qb(ct));
	if (!qb)
		goto qb_alloc_fail;

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
	destroyQureg(*qureg, env_get_questenv());
	free(qureg);
qureg_alloc_fail:
	return -1;
}

void circ_destroy(struct circ *c)
{
	if (c->qb) {
		const Qureg *qureg = c->qb;
		destroyQureg(*qureg, env_get_questenv());
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
}

void circ_report(struct circ const *c)
{
	env_report();

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
