#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "circ.h"
#include "common.h"
#include "qreg.h"

struct circ;

struct circuit {
	size_t num_mea_qb;
	size_t num_sys_qb;
	size_t num_anc_qb;

	int (*reset)(struct circ *);

	int (*prepst)(struct circ *);
	int (*effect)(struct circ *);
	int (*measure)(struct circ *);
};

/*
 * Create instance of a circuit.
 *
 * Using the circuit description `ct`, allocate memory for a new circuit
 * instance.  The pointer `data` will be stored and will be available during
 * call to `circ_simulate()`.
 */
struct circ *circ_create(struct circuit *ct, void *data);

/*
 * Destroy circuit instance.
 *
 * Free allocated memory.  After the function call, the pointer `c` will no
 * longer point to a valid structure and must be discarded.
 */
void circ_destroy(struct circ *c);

/*
 * Retrieve generic pointer to user data.
 */
void *circ_data(const struct circ *c);

/*
 * Reset circuit.
 *
 * Set classical register to zero and call `reset()` function specified by
 * the circuit template (if any).
 *
 * Return value:	 0	if successful
 *			-1	in case of failure
 */
int circ_reset(struct circ *c);

/*
 * Run circuit simulation.
 *
 * This will reset the circuit with circ_reset() and call functions specified by
 * the circuit template in the following order:
 *	prepst(),
 *	effect(),
 *	measure()
 *
 * Return value:	 0	if successful
 *			-1	in case of failure
 */
int circ_run(struct circ *c);

/*
 * Return number of qubits in the "measurement" register.
 */
size_t circ_num_meaqb(const struct circ *c);

/*
 * Return number of qubits in the "system" register.
 */
size_t circ_num_sysqb(const struct circ *c);

void circ_ops_blank(struct circ *c);

void circ_ops_setsysamp(struct circ *c, size_t idx, _Complex double amp);

void circ_ops_getsysamp(struct circ *c, size_t idx, _Complex double *amp);

void circ_ops_paulirot(struct circ *c, struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, size_t num_codes);

/*
 * Hamiltonian module.
 */

/*
 * Initialize Hamiltonian.
 *
 * This function must be called first before any other operation on the
 * structure is performed.
 */
void circ_hamil_init(struct circ_hamil *h);

/*
 * Destroy Hamiltonian.
 *
 * Free allocated memory.
 */
void circ_hamil_destroy(struct circ_hamil *h);

/*
 * Parse data from the open file represented by `fid` descriptor.
 *
 * Return value:  	 0	Data was parsed successfully
 *			-1	Error while reading data
 */
int circ_hamil_from_data2(struct circ_hamil *h, data2_id fid);

/*
 * Retrieve a single Pauli string corresponding to the index `n` in the sum of
 * terms.
 *
 * The value of `n` must be smaller than `h->num_qubits`.  The array `paulis`
 * must be of size at least `n`.
 *
 * After the call to this function, the value of `paulis` will be overwritten
 * with numbers specifying operators in the Pauli string:
 *	0	- I
 *	1	- X
 *	2	- Y
 *	3	- Z
 */
void circ_hamil_paulistr(const struct circ_hamil *h, size_t n, int *paulis);

struct circ {
	struct circuit *ct;
	void	       *data;
	int	       *cl, *qb;

	/* Qubit register */
	struct qreg reg;
};

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

	struct qreg reg;
	qreg_init(&reg, num_qb_tot);

	c->ct	= ct;
	c->data = data;
	c->reg	= reg;
	c->cl	= cl;
	c->qb	= qb;

	return c;

	// free(qb);
qb_fail:
	free(cl);
cl_fail:
	free(c);
circ_fail:
	return NULL;
}

void circ_destroy(struct circ *c)
{
	qreg_destroy(&c->reg);
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

int circ_reset(struct circ *c)
{
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

static const size_t PAULI_MASK	   = 3;
static const size_t PAULI_WIDTH	   = 2;
static const size_t PAULI_PAK_SIZE = sizeof(u64) * 8 / PAULI_WIDTH;

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
	size_t	idx;
	size_t	num_qubits;
	double	norm;
	double *coeffs;
	u64    *pak;
};

static int hamil_iter(double coeff, unsigned char *paulis, void *iter_data)
{
	struct hamil_iter_data *idat = iter_data;
	size_t			i = idat->idx++, num_qubits = idat->num_qubits;

	idat->coeffs[i] = coeff * idat->norm;
	for (size_t j = 0; j < num_qubits; j++) {
		ldiv_t	  dv	= ldiv(i * num_qubits + j, PAULI_PAK_SIZE);
		const u64 pauli = paulis[j];
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

	double *coeffs = malloc(sizeof *coeffs * num_terms);
	u64    *pak    = calloc(
		      num_terms * num_qubits / PAULI_PAK_SIZE + 1, sizeof(u64));
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

void circ_ops_blank(struct circ *c)
{
	qreg_zero(&c->reg);
}

void circ_ops_setsysamp(struct circ *c, size_t idx, _Complex double amp)
{
	double	  amps[2]   = { creal(amp), cimag(amp) };
	long long start_ind = idx << circ_num_meaqb(c);
	qreg_setamp(&c->reg, start_ind, amps);
}

void circ_ops_getsysamp(struct circ *c, size_t idx, _Complex double *amp)
{
	double	  amps[2];
	long long start_ind = idx << circ_num_meaqb(c);
	qreg_getamp(&c->reg, start_ind, &amps);

	*amp = amps[0] + _Complex_I * amps[1];
}

void circ_ops_paulirot(struct circ *c, const struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, const size_t num_codes)
{
	qreg_paulirot(&c->reg, code_hi, codes_lo, angles, num_codes);
}

#define MAX_CODES (1024)

struct code_cache {
	struct paulis code_hi;
	struct paulis codes_lo[MAX_CODES];
	fl	      angles[MAX_CODES];
	size_t	      num_codes;
};

struct circ_data {
	const struct circuit_data *rd;
	_Complex double		   prod;
	int			   scratch[64];

	struct code_cache cache;
};

void silk_multidet_init(struct circuit_multidet *md)
{
	md->num_dets = 0;
	md->dets     = NULL;
}

void silk_multidet_destroy(struct circuit_multidet *md)
{
	if (md->dets) {
		free(md->dets);
		md->dets = NULL;
	}
	md->num_dets = 0;
}

struct iter_multidet_data {
	size_t			 idx;
	struct circuit_multidet *md;
};

static int iter_multidet(_Complex double coeff, size_t idx, void *op_data)
{
	struct iter_multidet_data *imd = op_data;

	imd->md->dets[imd->idx].coeff = coeff;
	imd->md->dets[imd->idx].index = idx;
	imd->idx++;

	return 0;
}

int silk_multidet_from_data(struct circuit_multidet *md, const data2_id fid)
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

int circuit_data_init(struct circuit_data *rd, size_t num_steps)
{
	circ_hamil_init(&rd->hamil);
	silk_multidet_init(&rd->multidet);
	rd->num_steps = num_steps;

	rd->trotter_steps = malloc(sizeof *rd->trotter_steps * num_steps);
	if (!rd->trotter_steps)
		return -1;

	return 0;
}

void circuit_data_destroy(struct circuit_data *rd)
{
	silk_multidet_destroy(&rd->multidet);
	circ_hamil_destroy(&rd->hamil);

	free(rd->trotter_steps);
}

int circuit_data_from_file(struct circuit_data *rd, data2_id fid)
{
	int rc;

	rc = circ_hamil_from_data2(&rd->hamil, fid);
	rc |= silk_multidet_from_data(&rd->multidet, fid);
	data2_trotter_get_factor(fid, &rd->time_factor);

	return rc;
}

int silk_prepst(struct circ *c)
{
	const struct circ_data *cdat = circ_data(c);

	const struct circuit_multidet *md = &cdat->rd->multidet;

	circ_ops_blank(c);
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
		paulis_split(
			code, c->reg.qb_lo, c->reg.qb_hi, &code_lo, &code_hi);
		paulis_shr(&code_hi, c->reg.qb_lo);

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
	const struct circ_data *cdat = circ_data(c);
	const double		t    = cdat->rd->time_factor;
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

	const struct circuit_multidet *md = &cdat->rd->multidet;

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
	ct->num_mea_qb = 0;
	ct->num_sys_qb = num_sys_qb;
	ct->num_anc_qb = 0;
	ct->reset      = NULL;
	ct->prepst     = NULL;
	ct->effect     = silk_effect;
	ct->measure    = silk_measure;
}

static void silk_circuit_destroy(const struct circuit *ct)
{
	(void)ct;
}

int circuit_simulate(const struct circuit_data *rd)
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
