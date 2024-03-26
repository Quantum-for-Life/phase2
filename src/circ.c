#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "circ.h"
#include "qreg.h"

#define MAX_CACHE_CODES (1024)

#define HAMIL_PAULI_MASK (3)
#define HAMIL_PAULI_WIDTH (2)
#define HAMIL_PAULI_PAK_SIZE (32)

struct circ {
	size_t	    num_qb;
	struct qreg reg;

	const struct circ_data *data;
	_Complex double		prod;
	int			scratch[64];

	struct code_cache {
		struct paulis code_hi;
		struct paulis codes_lo[MAX_CACHE_CODES];
		fl	      angles[MAX_CACHE_CODES];
		size_t	      num_codes;
	} cache;
};

int circ_create(
	struct circ *c, const struct circ_data *data, const size_t num_qubits)
{
	struct qreg reg;
	if (qreg_init(&reg, num_qubits) < 0)
		return -1;

	c->num_qb = num_qubits;
	c->data	  = data;
	c->reg	  = reg;

	return 0;
}

void circ_destroy(struct circ *c)
{
	qreg_destroy(&c->reg);
}

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

static int hamil_iter(
	const double coeff, const unsigned char *paulis, void *iter_data)
{
	struct hamil_iter_data *idat	   = iter_data;
	const size_t		i	   = idat->idx++;
	const size_t		num_qubits = idat->num_qubits;

	idat->coeffs[i] = coeff * idat->norm;
	for (size_t j = 0; j < num_qubits; j++) {
		const ldiv_t dv =
			ldiv(i * num_qubits + j, HAMIL_PAULI_PAK_SIZE);
		const u64 pauli = paulis[j];
		idat->pak[dv.quot] += pauli << (dv.rem * HAMIL_PAULI_WIDTH);
	}

	return 0;
}

int circ_hamil_from_data2(struct circ_hamil *h, const data2_id fid)
{
	size_t num_qubits, num_terms;
	double norm;

	if (data2_hamil_getnums(fid, &num_qubits, &num_terms) < 0)
		return -1;
	if (data2_hamil_getnorm(fid, &norm) < 0)
		return -1;

	double *coeffs = malloc(sizeof *coeffs * num_terms);
	u64    *pak    = calloc(
		      num_terms * num_qubits / HAMIL_PAULI_PAK_SIZE + 1, sizeof(u64));
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

void circ_hamil_paulistr(
	const struct circ_hamil *h, const size_t n, int *paulis)
{
	for (size_t j = 0; j < h->num_qubits; j++) {
		const size_t pauli_idx = n * h->num_qubits + j;
		const ldiv_t dv	       = ldiv(pauli_idx, HAMIL_PAULI_PAK_SIZE);
		paulis[j] = (h->pak[dv.quot] >> (dv.rem * HAMIL_PAULI_WIDTH)) &
			    HAMIL_PAULI_MASK;
	}
}

void circ_reg_blank(struct circ *c)
{
	qreg_zero(&c->reg);
}

void circ_reg_setamp(
	struct circ *c, const size_t idx, const _Complex double amp)
{
	const double amps[2] = { creal(amp), cimag(amp) };
	qreg_setamp(&c->reg, idx, amps);
}

void circ_reg_getamp(struct circ *c, const size_t idx, _Complex double *amp)
{
	double amps[2];
	qreg_getamp(&c->reg, idx, &amps);

	*amp = amps[0] + _Complex_I * amps[1];
}

void circ_reg_paulirot(struct circ *c, const struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, const size_t num_codes)
{
	qreg_paulirot(&c->reg, code_hi, codes_lo, angles, num_codes);
}

void silk_multidet_init(struct circ_multidet *md)
{
	md->num_dets = 0;
	md->dets     = NULL;
}

void silk_multidet_destroy(struct circ_multidet *md)
{
	if (md->dets) {
		free(md->dets);
		md->dets = NULL;
	}
	md->num_dets = 0;
}

struct iter_multidet_data {
	uint64_t	      idx;
	struct circ_multidet *md;
};

static int iter_multidet(
	const _Complex double coeff, const uint64_t idx, void *op_data)
{
	struct iter_multidet_data *imd = op_data;

	imd->md->dets[imd->idx].coeff = coeff;
	imd->md->dets[imd->idx].idx = idx;
	imd->idx++;

	return 0;
}

int circuit_multidet_from_data(struct circ_multidet *md, const data2_id fid)
{
	size_t num_qubits, num_dets;
	if (data2_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		return -1;
	md->dets = malloc(sizeof *md->dets * num_dets);
	if (md->dets == NULL)
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

int circ_data_init(struct circ_data *cd, const size_t num_steps)
{
	circ_hamil_init(&cd->hamil);
	silk_multidet_init(&cd->multidet);
	cd->num_steps = num_steps;

	cd->trotter_steps = malloc(sizeof *cd->trotter_steps * num_steps);
	if (cd->trotter_steps == NULL)
		return -1;

	return 0;
}

void circ_data_destroy(struct circ_data *cd)
{
	silk_multidet_destroy(&cd->multidet);
	circ_hamil_destroy(&cd->hamil);

	free(cd->trotter_steps);
}

int circ_data_from_file(struct circ_data *cd, const data2_id fid)
{
	int rc = circ_hamil_from_data2(&cd->hamil, fid);
	rc |= circuit_multidet_from_data(&cd->multidet, fid);
	data2_trotter_get_factor(fid, &cd->time_factor);

	return rc;
}

int circuit_prepst(struct circ *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	circ_reg_blank(c);
	for (size_t i = 0; i < md->num_dets; i++) {
		circ_reg_setamp(c, md->dets[i].idx, md->dets[i].coeff);
	}

	return 0;
}

static void trotter_step(struct circ *c, const double omega)
{
	const struct circ_hamil *hamil	= &c->data->hamil;
	int			*paulis = c->scratch;

	struct code_cache cache = c->cache;
	cache.num_codes		= 0;

	for (size_t i = 0; i < hamil->num_terms; i++) {
		circ_hamil_paulistr(hamil, i, paulis);
		const double  angle = omega * hamil->coeffs[i];
		struct paulis code  = paulis_new();
		for (u32 k = 0; k < c->num_qb; k++)
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
			cache.num_codes < MAX_CACHE_CODES) {
			const size_t k	  = cache.num_codes++;
			cache.codes_lo[k] = code_lo;
			cache.angles[k]	  = angle;
			continue;
		}

		circ_reg_paulirot(c, cache.code_hi, cache.codes_lo,
			cache.angles, cache.num_codes);

		cache.num_codes	  = 1;
		cache.code_hi	  = code_hi;
		cache.codes_lo[0] = code_lo;
		cache.angles[0]	  = angle;
	}

	if (cache.num_codes > 0)
		circ_reg_paulirot(c, cache.code_hi, cache.codes_lo,
			cache.angles, cache.num_codes);
}

static int circ_effect(struct circ *c)
{
	const double t = c->data->time_factor;
	if (isnan(t))
		return -1;
	if (fabs(t) < DBL_EPSILON)
		return 0;

	trotter_step(c, t);

	return 0;
}

static int circ_measure(struct circ *c)
{
	const struct circ_multidet *md = &c->data->multidet;

	_Complex double prod = 0;
	for (size_t i = 0; i < md->num_dets; i++) {
		_Complex double amp;
		circ_reg_getamp(c, md->dets[i].idx, &amp);
		prod += amp * conj(md->dets[i].coeff);
	}

	c->prod = prod;
	return 0;
}

int circ_simulate(const struct circ_data *cd)
{
	int ret = 0;

	const size_t num_qb = cd->hamil.num_qubits;

	struct circ c;
	if (circ_create(&c, cd, num_qb) < 0)
		goto error;
	circuit_prepst(&c);

	for (size_t i = 0; i < cd->num_steps; i++) {
		if (circ_effect(&c) < 0)
			goto error;
		circ_measure(&c);
		cd->trotter_steps[i] = c.prod;
	}

	goto exit;
error:
	ret = -1;
exit:
	circ_destroy(&c);

	return ret;
}
