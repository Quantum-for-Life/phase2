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

	int *cl, *qb;

	/* Qubit register */
	Qureg quest_qureg;
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
int circ_initialize()
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

int circ_shutdown()
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
		if (circ_initialize() < 0)
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

struct circ *circ_create(struct circuit *ct, void *data)
{
	const size_t num_qb_tot =
		ct->num_mea_qb + ct->num_sys_qb + ct->num_anc_qb;
	struct circ *c = malloc(sizeof(*c));
	if (!c)
		goto circ_fail;
	int *cl = malloc(sizeof(*cl) * ct->num_mea_qb);
	if (!cl)
		goto cl_fail;
	int *qb = malloc(sizeof(*qb) * num_qb_tot);
	if (!qb)
		goto qb_fail;
	QuESTEnv *quest_env = env_get_questenv();
	if (!quest_env)
		goto quest_env_fail;
	Qureg qureg = createQureg(num_qb_tot, *quest_env);

	c->ct = ct;
	c->data = data;
	c->quest_qureg = qureg;
	c->cl = cl;
	c->qb = qb;

	return c;

quest_env_fail:
	free(qb);
qb_fail:
	free(cl);
cl_fail:
	free(c);
circ_fail:
	return NULL;
}

void circ_destroy(struct circ *c)
{
	QuESTEnv *quest_env = env_get_questenv();
	if (!quest_env) {
		return;
	}
	destroyQureg(c->quest_qureg, *quest_env);
	if (c->qb) {
		free(c->qb);
	}
	if (c->cl) {
		free(c->cl);
	}
	free(c);
}

void *circ_data(const struct circ *c)
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

	printf("num_mea_qb: %zu\n", circ_num_mea_qb(c));
	printf("num_sys_qb: %zu\n", circ_num_sys_qb(c));
	printf("num_anc_qb: %zu\n", circ_num_anc_qb(c));

	printf("cl register: [");
	for (size_t i = 0; i < c->ct->num_mea_qb; i++) {
		printf("%d", c->cl[i]);
	}
	printf("]\n");

	printf("----------------\n");

	return 0;
}

int circ_reset(struct circ *c)
{
	initZeroState(c->quest_qureg);
	for (size_t i = 0; i < circ_num_mea_qb(c); i++) {
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

size_t circ_num_mea_qb(const struct circ *c)
{
	return c->ct->num_mea_qb;
}

size_t circ_num_sys_qb(const struct circ *c)
{
	return c->ct->num_sys_qb;
}

size_t circ_num_anc_qb(const struct circ *c)
{
	return c->ct->num_anc_qb;
}

qbid circ_mea_qb(const struct circ *c, size_t idx)
{
	(void)c;
	return idx;
}

qbid circ_sys_qb(const struct circ *c, size_t idx)
{
	return idx + circ_num_mea_qb(c);
}

qbid circ_anc_qb(const struct circ *c, size_t idx)
{
	return idx + circ_num_mea_qb(c) + circ_num_sys_qb(c);
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
	long long start_ind = idx << circ_num_mea_qb(c);
	setAmps(c->quest_qureg, start_ind, &amp_re, &amp_im, 1);
}

void circ_sys_control_rotate_pauli(struct circ *c, int *paulis, double angle)
{
	const size_t num_mea_qb = circ_num_mea_qb(c);
	const size_t num_sys_qb = circ_num_sys_qb(c);
	for (size_t i = 0; i < num_mea_qb + num_sys_qb; i++) {
		c->qb[i] = i;
	}
	int *qb_ctl = c->qb;
	int *qb_trg = c->qb + num_mea_qb;

	multiControlledMultiRotatePauli(c->quest_qureg, qb_ctl, num_mea_qb,
					qb_trg, (enum pauliOpType *)paulis,
					num_sys_qb, -2.0 * angle);
}