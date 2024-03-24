#include <stdatomic.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "ev.h"

#include "circ.h"
#include "circ_private.h"

/* ----------------------------------------------------------------------------
 * Initialize and control global circ environment.
 *
 * This part was meant to provide safe initialization of a global environment in
 * a multithreaded setup.  Hence the use of atomic flags and memory ordering.
 *
 * Unfortunatelly, our main cluster setup (ETHZ Euler) does not support C11
 * threads. Hence, we put on hold the development of a concurrent simulation
 * environment using bare C threads, and focus insead on OpenMP/MPI concurrency.
 *
 * The code below should work just as well in a single-threaded situation, so
 * there's no reason to remove it.
 * -------------------------------------------------------------------------- */
static struct {
	_Atomic _Bool init;
	atomic_flag   lock;
	struct ev      ev;
} circ_env = {
	.init = false,
	.lock = ATOMIC_FLAG_INIT,
};

#ifdef __STDC_NO_THREADS__
void thrd_yield(void)
{
}
#else
#include <threads.h>
#endif

int circ_initialize(void)
{
	int rc = 1;

	/* If circ_env.init flag is set, someone else has already successfully
	 * initialized the environment, or is in the process of doing it. */
	if (atomic_load_explicit(&circ_env.init, memory_order_relaxed))
		return rc;

	/* Try to obtain a spinlock. Enter critical section until the lock is
	released. */
	while (atomic_flag_test_and_set(&circ_env.lock))
		thrd_yield();
	/* Check if someone hasn't set everything up for us while we were
	waiting.  If the flag is not set, we need to be sure other members
	of circ_env haven't been touched yet either.  Hence `acquire`
	memory order. */
	if (!atomic_load_explicit(&circ_env.init, memory_order_acquire)) {
		/* The call to QuEST always succeeds. If there's something else
		here to do that might fail, conditionally set the return code
		rc=-1 and leave the flag off. */
		ev_init(&circ_env.ev);
		int atex_rc	   = atexit(circ_shutdown);
		if (atex_rc == 0) {
			atomic_store_explicit(
				&circ_env.init, true, memory_order_release);
			rc = 0;
		} else {
			rc = -1;
		}
	}
	atomic_flag_clear(&circ_env.lock);

	return rc;
}

void circ_shutdown(void)
{
	if (!atomic_load_explicit(&circ_env.init, memory_order_relaxed))
		return;

	while (atomic_flag_test_and_set(&circ_env.lock))
		thrd_yield();
	if (atomic_load_explicit(&circ_env.init, memory_order_acquire)) {
		ev_destroy(&circ_env.ev);
		atomic_store_explicit(
			&circ_env.init, false, memory_order_release);
	}
	atomic_flag_clear(&circ_env.lock);
}

static struct ev *env_get_questenv(void)
{
	if (!atomic_load_explicit(&circ_env.init, memory_order_acquire)) {
		if (circ_initialize() < 0)
			return NULL;
	}

	return &circ_env.ev;
}

static int env_report(void)
{
	const struct ev *quest_env = env_get_questenv();
	if (!quest_env) {
		fprintf(stderr, "Error: circuit environment not initialized\n");
		return -1;
	}
	//reportQuESTEnv(*quest_env);

	return 0;
}

struct circ *circ_create(struct circuit *ct, void *data)
{
	struct circ *c = malloc(sizeof(*c));
	if (!c)
		goto circ_fail;
	int *cl = malloc(sizeof(*cl) * ct->num_mea_qb);
	if (!cl)
		goto cl_fail;
	const size_t num_qb_tot =
		ct->num_mea_qb + ct->num_sys_qb + ct->num_anc_qb;
	int *qb = malloc(sizeof(*qb) * num_qb_tot);
	if (!qb)
		goto qb_fail;
	struct ev *ev = env_get_questenv();
	if (!ev)
		goto ev_fail;

	struct qreg reg;
	qreg_init(&reg, num_qb_tot,ev);



	c->ct	       = ct;
	c->data	       = data;
	c->quest_qureg = reg;
	c->cl	       = cl;
	c->qb	       = qb;

	return c;

ev_fail:
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
	const struct ev *quest_env = env_get_questenv();
	if (!quest_env)
		return;
	qreg_destroy(&c->quest_qureg);
	if (c->qb)
		free(c->qb);
	if (c->cl)
		free(c->cl);
	free(c);
	c = NULL;
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

	//reportQuregParams(c->quest_qureg);

	printf("num_mea_qb: %zu\n", circ_num_meaqb(c));
	printf("num_sys_qb: %zu\n", circ_num_sysqb(c));
	printf("num_anc_qb: %zu\n", circ_num_ancqb(c));

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
	// initZeroState(c->quest_qureg);
	for (size_t i = 0; i < circ_num_meaqb(c); i++) {
		c->cl[i] = 0;
	}
	if (c->ct->reset)
		return c->ct->reset(c);

	return 0;
}

int circ_run(struct circ *c)
{
	if (circ_reset(c) < 0)
		return -1;

	int (*ops[3])(struct circ *) = { c->ct->prepst, c->ct->effect,
		c->ct->measure };
	for (int i = 0; i < 3; i++) {
		if (ops[i] && ops[i](c) < 0)
			return -1;
	}

	return 0;
}

size_t circ_num_meaqb(const struct circ *c)
{
	return c->ct->num_mea_qb;
}

size_t circ_num_sysqb(const struct circ *c)
{
	return c->ct->num_sys_qb;
}

size_t circ_num_ancqb(const struct circ *c)
{
	return c->ct->num_anc_qb;
}

qbid circ_meaqb(const struct circ *c, size_t idx)
{
	(void)c;
	return idx;
}

qbid circ_sysqb(const struct circ *c, size_t idx)
{
	return idx + circ_num_meaqb(c);
}

qbid circ_ancqb(const struct circ *c, size_t idx)
{
	return idx + circ_num_meaqb(c) + circ_num_sysqb(c);
}
