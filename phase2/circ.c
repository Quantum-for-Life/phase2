#include "c23_compat.h"
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/paulis.h"

int circ_cache_init(struct circ_cache *ch, size_t qb_lo, size_t qb_hi)
{
	struct paulis *lo = malloc(sizeof(struct paulis) * MAX_CACHE_CODES);
	if (!lo)
		goto err_lo;
	double *angles = malloc(sizeof(double) * MAX_CACHE_CODES);
	if (!angles)
		goto err_angles;

	ch->codes_lo = lo;
	ch->phis = angles;
	ch->qb_lo = qb_lo;
	ch->qb_hi = qb_hi;
	ch->n = 0;

	return 0;

	// free(phs);
err_angles:
	free(lo);
err_lo:
	return -1;
}

void circ_cache_destroy(struct circ_cache *ch)
{
	free(ch->codes_lo);
	free(ch->phis);
}

int circ_cache_insert(struct circ_cache *ch, struct paulis code, double phi)
{
	struct paulis lo, hi;
	paulis_split(code, ch->qb_lo, ch->qb_hi, &lo, &hi);
	if (ch->n == 0) {
		ch->code_hi = hi;
		ch->codes_lo[0] = lo;
		ch->phis[0] = phi;
		ch->n = 1;
		return 0;
	}

	if (ch->n < MAX_CACHE_CODES && paulis_eq(ch->code_hi, hi)) {
		const size_t k = ch->n++;
		ch->codes_lo[k] = lo;
		ch->phis[k] = phi;
		return 0;
	}

	return -1;
}

void circ_cache_flush(struct circ_cache *ch,
	void (*op)(struct paulis, struct paulis *, double *, size_t, void *),
	void *data)
{
	if (ch->n > 0 && op)
		op(ch->code_hi, ch->codes_lo, ch->phis, ch->n, data);
	ch->n = 0;
}

static int circ_hamil_init(struct circ_hamil *h, uint32_t nqb, size_t nterms)
{
	h->terms = malloc(sizeof *h->terms * nterms);
	if (!h->terms)
		return -1;
	h->nterms = nterms;
	h->nqb = nqb;

	return 0;
}

static void circ_hamil_destroy(struct circ_hamil *h)
{
	if (h->terms != nullptr)
		free(h->terms);
}

struct hamil_iter_data {
	size_t i;
	double norm;
	struct circ_hamil *hamil;
};

static int hamil_iter(double cf, unsigned char *ops, void *iter_data)
{
	struct hamil_iter_data *id = iter_data;
	struct circ_hamil *h = id->hamil;
	const size_t i = id->i++;
	const uint32_t nqb = h->nqb;

	h->terms[i].cf = cf * id->norm;
	struct paulis op = paulis_new();
	for (uint32_t j = 0; j < nqb; j++)
		paulis_set(&op, ops[j], j);
	h->terms[i].op = op;

	return 0;
}

static int circ_hamil_from_file(struct circ_hamil *h, const data_id fid)
{
	uint32_t nqb;
	size_t nterms;
	double norm;

	if (data_hamil_getnums(fid, &nqb, &nterms) < 0)
		return -1;
	if (data_hamil_getnorm(fid, &norm) < 0)
		return -1;
	if (circ_hamil_init(h, nqb, nterms) < 0)
		return -1;

	struct hamil_iter_data id = { .i = 0, .norm = norm, .hamil = h };
	if (data_hamil_foreach(fid, hamil_iter, &id) != 0)
		return -1;

	return 0;
}

static int circ_muldet_init(struct circ_muldet *m, size_t ndets)
{
	m->dets = malloc(sizeof *m->dets * ndets);
	if (!m->dets)
		return -1;
	m->ndets = ndets;

	return 0;
}

static void circ_muldet_destroy(struct circ_muldet *m)
{
	if (m->dets != nullptr)
		free(m->dets);
}

struct iter_muldet_data {
	size_t i;
	struct circ_muldet *muldet;
};

static int iter_muldet(_Complex double cf, const uint64_t idx, void *iter_data)
{
	struct iter_muldet_data *id = iter_data;
	struct circ_muldet *m = id->muldet;
	const size_t i = id->i++;

	m->dets[i].cf = cf;
	m->dets[i].idx = idx;

	return 0;
}

static int circ_muldet_from_file(struct circ_muldet *m, const data_id fid)
{
	uint32_t nqb;
	size_t ndets;

	if (data_multidet_getnums(fid, &nqb, &ndets) < 0)
		return -1;
	if (circ_muldet_init(m, ndets) < 0)
		return -1;

	struct iter_muldet_data id = { .i = 0, .muldet = m };
	if (data_multidet_foreach(fid, iter_muldet, &id) < 0)
		return -1;

	return 0;
}

int circ_init(
	struct circ *ct, const data_id fid, int (*simul)(struct circ *))
{
	if (circ_hamil_from_file(&ct->hamil, fid) < 0)
		goto err_hamil_init;
	if (circ_muldet_from_file(&ct->muldet, fid) < 0)
		goto err_muldet_init;
	if (qreg_init(&ct->reg, ct->hamil.nqb) < 0)
		goto err_qreg_init;
	if (circ_cache_init(&ct->cache, ct->reg.qb_lo, ct->reg.qb_hi) < 0)
		goto err_cache_init;

	ct->simul = simul;

	return 0;

	// circ_cache_destroy(&c->cache);
err_cache_init:
	qreg_destroy(&ct->reg);
err_qreg_init:
	circ_muldet_destroy(&ct->muldet);
err_muldet_init:
	circ_hamil_destroy(&ct->hamil);
err_hamil_init:
	return -1;
}

void circ_destroy(struct circ *ct)
{
	circ_hamil_destroy(&ct->hamil);
	circ_muldet_destroy(&ct->muldet);
	circ_cache_destroy(&ct->cache);
	qreg_destroy(&ct->reg);
}

int circ_prepst(struct circ *ct)
{
	const struct circ_muldet *md = &ct->muldet;

	qreg_zero(&ct->reg);
	for (size_t i = 0; i < md->ndets; i++)
		qreg_setamp(&ct->reg, md->dets[i].idx, md->dets[i].cf);

	return 0;
}

static void circ_flush(struct paulis code_hi, struct paulis *codes_lo,
	double *phis, size_t ncodes, void *data)
{
	struct qreg *reg = data;
	qreg_paulirot(reg, code_hi, codes_lo, phis, ncodes);
}

int circ_step(struct circ *ct, const struct circ_hamil *hm, const double omega)
{
	struct circ_cache *cache = &ct->cache;

	for (size_t i = 0; i < hm->nterms; i++) {
		const double phi = omega * hm->terms[i].cf;
		const struct paulis code = hm->terms[i].op;

		if (circ_cache_insert(cache, code, phi) == 0)
			continue;

		log_trace("paulirot, term: %zu, num_codes: %zu", i, cache->n);
		circ_cache_flush(cache, circ_flush, &ct->reg);
		if (circ_cache_insert(cache, code, phi) < 0)
			return -1;
	}
	log_trace("paulirot, last term group, num_codes: %zu", cache->n);
	circ_cache_flush(cache, circ_flush, &ct->reg);

	return 0;
}

_Complex double circ_measure(struct circ *ct)
{
	const struct circ_muldet *md = &ct->muldet;

	_Complex double pr = 0.0;
	for (size_t i = 0; i < md->ndets; i++) {
		_Complex double a;
		qreg_getamp(&ct->reg, md->dets[i].idx, &a);
		pr += a * conj(md->dets[i].cf);
	}

	return pr;
}

inline int circ_simul(struct circ *ct)
{
	return ct->simul(ct);
}

static int hamil_term_cmp_lex(const void *a, const void *b)
{
	const struct paulis x = ((const struct circ_hamil_term *)a)->op;
	const struct paulis y = ((const struct circ_hamil_term *)b)->op;

	return paulis_cmp(x, y);
}

void circ_hamil_sort_lex(struct circ_hamil *hm)
{
	qsort(hm->terms, hm->nterms, sizeof(struct circ_hamil_term),
		hamil_term_cmp_lex);
}
