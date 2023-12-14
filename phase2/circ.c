#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>

#include "QuEST.h"

#include "circ.h"

struct circ {
	struct circuit *ct;
	void *data;

	/* Qubit register */
	void *reg;

	int *mea_cl;
	int *mea_qb;
	int *sys_qb;
	int *anc_qb;
};

static struct {
	_Bool init;
	QuESTEnv quest_env;
} circ_env = { .init = false };

static void env_init()
{
	circ_env.quest_env = createQuESTEnv();
	circ_env.init = true;
}

static void env_destroy()
{
	if (circ_env.init) {
		destroyQuESTEnv(circ_env.quest_env);
		circ_env.init = false;
	}
}

static QuESTEnv get_questenv(void)
{
	if (!circ_env.init) {
		env_init(&circ_env);
	}
	return circ_env.quest_env;
}

static void env_report(void)
{
	reportQuESTEnv(get_questenv());
}

static size_t ct_num_tot_qb(const struct circuit *ct)
{
	return ct->num_mea_qb + ct->num_sys_qb + ct->num_anc_qb;
}

struct circ *circ_init(struct circuit *ct, void *data)
{
	struct circ *c = malloc(sizeof(*c));
	if (!c)
		return NULL;

	Qureg *qureg = malloc(sizeof(Qureg));
	if (!qureg)
		goto qureg_alloc_fail;

	*qureg = createQureg(ct_num_tot_qb(ct), get_questenv());

	int *mea_cl = malloc(sizeof(int) * ct->num_mea_qb);
	if (!mea_cl)
		goto mea_alloc_fail;
	int *qb = malloc(sizeof(int) * ct_num_tot_qb(ct));
	if (!qb)
		goto qb_alloc_fail;

	c->ct = ct;
	c->data = data;
	c->reg = qureg;
	c->mea_cl = mea_cl;
	c->mea_qb = qb;
	c->sys_qb = c->mea_qb + ct->num_mea_qb;
	c->anc_qb = c->sys_qb + ct->num_sys_qb;
	for (size_t i = 0; i < ct_num_tot_qb(c->ct); i++) {
		c->mea_qb[i] = i;
	}

	return c;

qb_alloc_fail:
	free(mea_cl);
mea_alloc_fail:
	destroyQureg(*qureg, get_questenv());
	free(qureg);
qureg_alloc_fail:
	return NULL;
}

void circ_destroy(struct circ *c)
{
	if (c->reg) {
		const Qureg *qureg = c->reg;
		destroyQureg(*qureg, get_questenv());
		free(c->reg);
	}
	if (c->mea_qb) {
		free(c->mea_qb);
	}
	if (c->mea_cl) {
		free(c->mea_cl);
	}
	free(c);
}

void *circ_data(struct circ *c)
{
	return c->data;
}

void circ_report(struct circ const *c)
{
	env_report();

	const Qureg *qureg = c->reg;
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
	Qureg *qureg = c->reg;
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

	int (*ops[3])(struct circ *) = {
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

qbid circ_mea_qb(struct circ *c, size_t idx)
{
	return c->mea_qb[idx];
}

qbid circ_sys_qb(struct circ *c, size_t idx)
{
	return c->sys_qb[idx];
}

qbid circ_anc_qb(struct circ *c, size_t idx)
{
	return c->anc_qb[idx];
}

void circ_hadamard(struct circ *c, qbid qb)
{
	Qureg *qureg = c->reg;
	if (qureg)
		hadamard(*qureg, qb);
}

void circ_sgate(struct circ *c, qbid qb)
{
	Qureg *qureg = c->reg;
	if (qureg)
		sGate(*qureg, qb);
}

double circ_prob0(struct circ *c, qbid qb)
{
	Qureg *qureg = c->reg;
	return calcProbOfOutcome(*qureg, qb, 0);
}

void circ_blankstate(struct circ *c)
{
	Qureg *qureg = c->reg;
	initBlankState(*qureg);
}

void circ_setsysamp(struct circ *c, size_t idx, _Complex double amp)
{
	Qureg *qureg = c->reg;
	double amp_re = creal(amp);
	double amp_im = cimag(amp);
	setAmps(*qureg, idx << c->ct->num_mea_qb, &amp_re, &amp_im, 1);
}

void circ_sys_control_rotate_pauli(struct circ *c, int *paulis, double angle)
{
	Qureg *qureg = c->reg;
	multiControlledMultiRotatePauli(*qureg, c->mea_qb, c->ct->num_mea_qb,
					c->sys_qb, (enum pauliOpType *)paulis,
					c->ct->num_sys_qb, -2.0 * angle);
}