#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <threads.h>

#include "QuEST.h"

#include "circ.h"

static struct {
	_Atomic _Bool init;
	atomic_flag lock;
	QuESTEnv quest_env;
} circ_env = {
	.init = false,
	.lock = ATOMIC_FLAG_INIT,
};

struct circ {
	struct circuit *ct;
	void *data;

	/* Qubit register */
	Qureg quest_qureg;

	int *cl;
	int *qb;
	int *sys_qb;
	int *anc_qb;
};

/** Initialize circuit environment.
 *
 * This function can be called multiple times from different threads.
 *
 * Return value:
 *   0 - if the environment was successfully initialized
 *       by this function call
 *   1 - if the environment has already been initialized
 *       by another call to this function
 *  -1 - in case of failure
 */
int circ_env_initialize()
{
	int rc = 1;

	/**
	 * If circ_env.init flag is set, someone else has already successfully
	 * initialized the environment, or is in the process of doing it.
	 */
	if (atomic_load_explicit(&circ_env.init, memory_order_relaxed))
		return rc;

	/**
	 * Try to obtain a spinlock. Enter critical section until the lock is released.
	 */
	while (atomic_flag_test_and_set(&circ_env.lock))
		thrd_yield();
	/**
	 * Check if someone hasn't set everything up for us while we were waiting.
	 * If the flag is not set, we need to be sure other members of circ_env
	 * haven't been touched yet either.  Hence `acquire` memory order.
	 */
	if (!atomic_load_explicit(&circ_env.init, memory_order_acquire)) {
		/**
		 * A call to QuEST always succeeds. If there's something else
		 * here to do that might fail, conditionally set the return code
		 * rc=-1 and leave the flag off.
		 */
		circ_env.quest_env = createQuESTEnv();
		if (true) {
			atomic_store_explicit(&circ_env.init, true,
					      memory_order_release);
			rc = 0;
		} else {
			rc = -1;
		}
	}
	atomic_flag_clear(&circ_env.lock);

	return rc;
}

int circ_env_shutdown()
{
	int rc = 1;
	if (!atomic_load_explicit(&circ_env.init, memory_order_relaxed))
		return rc;

	while (atomic_flag_test_and_set(&circ_env.lock))
		thrd_yield();
	if (atomic_load_explicit(&circ_env.init, memory_order_acquire)) {
		destroyQuESTEnv(circ_env.quest_env);
		atomic_store_explicit(&circ_env.init, false,
				      memory_order_release);
		rc = 0;
	}
	atomic_flag_clear(&circ_env.lock);

	return rc;
}

static QuESTEnv *env_get_questenv(void)
{
	if (!atomic_load_explicit(&circ_env.init, memory_order_acquire)) {
		if (circ_env_initialize() < 0)
			return NULL;
	}
	return &circ_env.quest_env;
}

static int env_report(void)
{
	QuESTEnv *quest_env = env_get_questenv();
	if (!quest_env) {
		fprintf(stderr, "Error: circuit environment not initialized");
		return -1;
	}
	reportQuESTEnv(*quest_env);

	return 0;
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

	int *mea_cl = malloc(sizeof(int) * ct->num_mea_qb);
	if (!mea_cl)
		goto mea_alloc_fail;
	int *qb = malloc(sizeof(int) * ct_num_tot_qb(ct));
	if (!qb)
		goto qb_alloc_fail;
	Qureg qureg = createQureg(ct_num_tot_qb(ct), *env_get_questenv());

	c->ct = ct;
	c->data = data;
	c->quest_qureg = qureg;
	c->cl = mea_cl;
	c->qb = qb;
	c->sys_qb = c->qb + ct->num_mea_qb;
	c->anc_qb = c->sys_qb + ct->num_sys_qb;
	for (size_t i = 0; i < ct_num_tot_qb(c->ct); i++) {
		c->qb[i] = i;
	}

	return c;

qb_alloc_fail:
	free(mea_cl);
mea_alloc_fail:
	destroyQureg(qureg, *env_get_questenv());
	return NULL;
}

void circ_destroy(struct circ *c)
{
	destroyQureg(c->quest_qureg, *env_get_questenv());
	if (c->qb) {
		free(c->qb);
	}
	if (c->cl) {
		free(c->cl);
	}
	free(c);
}

void *circ_data(struct circ *c)
{
	return c->data;
}

int circ_report(struct circ const *c)
{
	if (env_report() < 0)
		return -1;

	printf("----------------\n");
	printf("CIRCUIT: %s\n", c->ct->name);
	reportQuregParams(c->quest_qureg);

	printf("cl register: [");
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		printf("%d", c->cl[i]);
	}
	printf("]\n");

	printf("qb indices: { ");
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		printf("%d ", c->qb[i]);
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

	return 0;
}

int circ_reset(struct circ *c)
{
	initZeroState(c->quest_qureg);
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		c->cl[i] = 0;
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
	return c->qb[idx];
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
	hadamard(c->quest_qureg, qb);
}

void circ_sgate(struct circ *c, qbid qb)
{
	sGate(c->quest_qureg, qb);
}

double circ_prob0(struct circ *c, qbid qb)
{
	return calcProbOfOutcome(c->quest_qureg, qb, 0);
}

void circ_blankstate(struct circ *c)
{
	initBlankState(c->quest_qureg);
}

void circ_setsysamp(struct circ *c, size_t idx, _Complex double amp)
{
	double amp_re = creal(amp);
	double amp_im = cimag(amp);
	setAmps(c->quest_qureg, idx << c->ct->num_mea_qb, &amp_re, &amp_im, 1);
}

void circ_sys_control_rotate_pauli(struct circ *c, int *paulis, double angle)
{
	multiControlledMultiRotatePauli(c->quest_qureg, c->qb,
					c->ct->num_mea_qb, c->sys_qb,
					(enum pauliOpType *)paulis,
					c->ct->num_sys_qb, -2.0 * angle);
}