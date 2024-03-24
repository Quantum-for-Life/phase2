#include <complex.h>
#include <float.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "qreg.h"
#include "circ.h"

#include <math.h>

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
	struct ev     ev;
} circ_env = {
	.init = false,
	.lock = ATOMIC_FLAG_INIT,
};

struct circ {
	struct circuit *ct;
	void	       *data;
	int	       *cl, *qb;

	/* Qubit register */
	struct qreg quest_qureg;
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
		int atex_rc = atexit(circ_shutdown);
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
	// reportQuESTEnv(*quest_env);

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
	qreg_init(&reg, num_qb_tot, ev);

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

	// reportQuregParams(c->quest_qureg);

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

static const size_t PAULI_MASK	   = 3;
static const size_t PAULI_WIDTH	   = 2;
static const size_t PAULI_PAK_SIZE = sizeof(pauli_pak_t) * 8 / PAULI_WIDTH;

void circ_hamil_init(struct circ_hamil *h)
{
	h->num_qubits = 0;
	h->num_terms  = 0;
	h->coeffs     = NULL;
	h->pak	      = NULL;
}

void circ_hamil_destroy(struct circ_hamil *h)
{
	if (h->pak) {
		free(h->pak);
		h->pak = NULL;
	}
	if (h->coeffs) {
		free(h->coeffs);
		h->coeffs = NULL;
	}
	h->num_terms  = 0;
	h->num_qubits = 0;
}

struct hamil_iter_data {
	size_t	     idx;
	size_t	     num_qubits;
	double	     norm;
	double	    *coeffs;
	pauli_pak_t *pak;
};

static int hamil_iter(double coeff, unsigned char *paulis, void *iter_data)
{
	struct hamil_iter_data *idat = iter_data;
	size_t			i = idat->idx++, num_qubits = idat->num_qubits;

	idat->coeffs[i] = coeff * idat->norm;
	for (size_t j = 0; j < num_qubits; j++) {
		ldiv_t		  dv = ldiv(i * num_qubits + j, PAULI_PAK_SIZE);
		const pauli_pak_t pauli = paulis[j];
		idat->pak[dv.quot] += pauli << (dv.rem * PAULI_WIDTH);
	}

	return 0;
}

int circ_hamil_from_data2(struct circ_hamil *h, data2_id fid)
{
	size_t num_qubits, num_terms;
	double norm;

	if (data2_hamil_getnums(fid, &num_qubits, &num_terms) < 0)
		return -1;
	if (data2_hamil_getnorm(fid, &norm) < 0)
		return -1;

	double	    *coeffs = malloc(sizeof *coeffs * num_terms);
	pauli_pak_t *pak = calloc(num_terms * num_qubits / PAULI_PAK_SIZE + 1,
		sizeof(pauli_pak_t));
	if (!(coeffs && pak))
		goto err;

	struct hamil_iter_data idat = { .idx = 0,
		.num_qubits		     = num_qubits,
		.norm			     = norm,
		.coeffs			     = coeffs,
		.pak			     = pak };
	if (data2_hamil_foreach(fid, hamil_iter, &idat) != 0)
		goto err;

	h->num_qubits = num_qubits;
	h->num_terms  = num_terms;
	h->coeffs     = coeffs;
	h->pak	      = pak;

	return 0;
err:
	free(coeffs);
	free(pak);
	return -1;
}

void circ_hamil_paulistr(const struct circ_hamil *h, size_t n, int *paulis)
{
	for (size_t j = 0; j < h->num_qubits; j++) {
		const size_t pauli_idx = n * h->num_qubits + j;
		ldiv_t	     dv	       = ldiv(pauli_idx, PAULI_PAK_SIZE);
		paulis[j] = (h->pak[dv.quot] >> (dv.rem * PAULI_WIDTH)) &
			    PAULI_MASK;
	}
}

double circ_ops_prob0(struct circ *c, qbid qb)
{
	return 0.0; // calcProbOfOutcome(c->quest_qureg, qb, 0);
}

void circ_ops_blankstate(struct circ *c)
{
	qreg_blankstate(&c->quest_qureg);
}

void circ_ops_setsysamp(struct circ *c, size_t idx, _Complex double amp)
{
	double	  amps[2]   = { creal(amp), cimag(amp) };
	long long start_ind = idx << circ_num_meaqb(c);
	qreg_setamp(&c->quest_qureg, start_ind, amps);
}

void circ_ops_getsysamp(struct circ *c, size_t idx, _Complex double *amp)
{
	double	  amps[2];
	long long start_ind = idx << circ_num_meaqb(c);
	qreg_getamp(&c->quest_qureg, start_ind, &amps);

	*amp = amps[0] + _Complex_I * amps[1];
}

void circ_ops_paulirot(struct circ *c, const struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, const size_t num_codes)
{
	qreg_paulirot(&c->quest_qureg, code_hi, codes_lo, angles, num_codes);
}

#define MAX_CODES (1024)

struct code_cache {
	struct paulis code_hi;
	struct paulis codes_lo[MAX_CODES];
	fl	      angles[MAX_CODES];
	size_t	      num_codes;
};

struct circ_data {
	const struct silk_data *rd;
	_Complex double		prod;
	int			scratch[64];

	struct code_cache cache;
};

void silk_multidet_init(struct silk_data_multidet *md)
{
	md->num_dets = 0;
	md->dets     = NULL;
}

void silk_multidet_destroy(struct silk_data_multidet *md)
{
	if (md->dets) {
		free(md->dets);
		md->dets = NULL;
	}
	md->num_dets = 0;
}

struct iter_multidet_data {
	size_t			   idx;
	struct silk_data_multidet *md;
};

static int iter_multidet(_Complex double coeff, size_t idx, void *op_data)
{
	struct iter_multidet_data *imd = op_data;

	imd->md->dets[imd->idx].coeff = coeff;
	imd->md->dets[imd->idx].index = idx;
	imd->idx++;

	return 0;
}

int silk_multidet_from_data(struct silk_data_multidet *md, const data2_id fid)
{
	size_t num_qubits, num_dets;
	if (data2_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		return -1;
	md->dets = malloc(sizeof *md->dets * num_dets);
	if (!md->dets)
		return -1;

	struct iter_multidet_data imd;
	imd.idx = 0;
	imd.md	= md;
	if (data2_multidet_foreach(fid, iter_multidet, &imd) < 0)
		goto error;

	md->num_dets = num_dets;

	return 0;
error:
	free(md->dets);
	return -1;
}

int silk_data_init(struct silk_data *rd, size_t num_steps)
{
	circ_hamil_init(&rd->hamil);
	silk_multidet_init(&rd->multidet);
	rd->num_steps = num_steps;

	rd->trotter_steps = malloc(sizeof *rd->trotter_steps * num_steps);
	if (!rd->trotter_steps)
		return -1;

	return 0;
}

void silk_data_destroy(struct silk_data *rd)
{
	silk_multidet_destroy(&rd->multidet);
	circ_hamil_destroy(&rd->hamil);

	free(rd->trotter_steps);
}

int silk_data_from_data(struct silk_data *rd, data2_id fid)
{
	int rc;

	rc = circ_hamil_from_data2(&rd->hamil, fid);
	rc |= silk_multidet_from_data(&rd->multidet, fid);
	data2_trotter_get_factor(fid, &rd->time_factor);

	return rc;
}

int silk_prepst(struct circ *c)
{
	struct circ_data *cdat = circ_data(c);

	const struct silk_data_multidet *md = &cdat->rd->multidet;

	circ_ops_blankstate(c);
	for (size_t i = 0; i < md->num_dets; i++) {
		circ_ops_setsysamp(c, md->dets[i].index, md->dets[i].coeff);
	}

	return 0;
}

static void trotter_step(struct circ *c, double omega)
{
	struct circ_data	*cdat	= circ_data(c);
	const struct circ_hamil *hamil	= &cdat->rd->hamil;
	int			*paulis = cdat->scratch;

	struct code_cache cache = cdat->cache;
	cache.num_codes		= 0;

	for (size_t i = 0; i < hamil->num_terms; i++) {
		circ_hamil_paulistr(hamil, i, paulis);
		const double  angle = omega * hamil->coeffs[i];
		struct paulis code  = paulis_new();
		for (u32 k = 0; k < circ_num_sysqb(c); k++)
			paulis_set(&code, paulis[k], k);

		struct paulis code_hi, code_lo;
		paulis_split(code, c->quest_qureg.qb_lo, c->quest_qureg.qb_hi,
			&code_lo, &code_hi);
		paulis_shr(&code_hi, c->quest_qureg.qb_lo);

		if (cache.num_codes == 0) {
			cache.code_hi	  = code_hi;
			cache.codes_lo[0] = code_lo;
			cache.angles[0]	  = angle;
			cache.num_codes++;
			continue;
		}

		if (paulis_eq(cache.code_hi, code_hi) &&
			cache.num_codes < MAX_CODES) {
			const size_t k	  = cache.num_codes++;
			cache.codes_lo[k] = code_lo;
			cache.angles[k]	  = angle;
			continue;
		}

		circ_ops_paulirot(c, cache.code_hi, cache.codes_lo,
			cache.angles, cache.num_codes);

		cache.num_codes	  = 1;
		cache.code_hi	  = code_hi;
		cache.codes_lo[0] = code_lo;
		cache.angles[0]	  = angle;
	}

	if (cache.num_codes > 0)
		circ_ops_paulirot(c, cache.code_hi, cache.codes_lo,
			cache.angles, cache.num_codes);
}

int silk_effect(struct circ *c)
{
	struct circ_data *cdat = circ_data(c);
	const double	  t    = cdat->rd->time_factor;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	trotter_step(c, t);

	return 0;
}

int silk_measure(struct circ *c)
{
	struct circ_data *cdat = circ_data(c);

	const struct silk_data_multidet *md = &cdat->rd->multidet;

	_Complex double prod = 0;
	for (size_t i = 0; i < md->num_dets; i++) {
		_Complex double amp;
		circ_ops_getsysamp(c, md->dets[i].index, &amp);
		prod += amp * conj(md->dets[i].coeff);
	}

	cdat->prod = prod;
	return 0;
}

static void silk_circuit_init(struct circuit *ct, size_t num_sys_qb)
{
	ct->name       = SILK_NAME;
	ct->num_mea_qb = SILK_NUM_MEA_QB;
	ct->num_sys_qb = num_sys_qb;
	ct->num_anc_qb = SILK_NUM_ANC_QB;
	ct->reset      = NULL;
	ct->prepst     = NULL;
	ct->effect     = silk_effect;
	ct->measure    = silk_measure;
}

static void silk_circuit_destroy(const struct circuit *ct)
{
	(void)ct;
}

int silk_simulate(const struct silk_data *rd)
{
	int ret = 0;

	struct circuit ct;
	silk_circuit_init(&ct, rd->hamil.num_qubits);

	struct circ_data cdat;
	cdat.rd	       = rd;
	struct circ *c = circ_create(&ct, &cdat);
	if (!c)
		goto error;
	silk_prepst(c);

	for (size_t i = 0; i < rd->num_steps; i++) {
		if (circ_run(c) < 0)
			goto error;
		rd->trotter_steps[i] = cdat.prod;
	}

	goto exit;
error:
	ret = -1;
exit:
	circ_destroy(c);
	silk_circuit_destroy(&ct);

	return ret;
}
