#include "c23_compat.h"
#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/paulis.h"

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

int circ_init(struct circ *c, data_id fid, void *data)
{
	if (circ_hamil_from_file(&c->hamil, fid) < 0)
		return -1;
	if (circ_muldet_from_file(&c->muldet, fid) < 0)
		return -1;
	c->data = data;
	if (circ_res_init(c) < 0)
		return -1;

	return 0;
}

void circ_destroy(struct circ *c)
{
	circ_hamil_destroy(&c->hamil);
	circ_muldet_destroy(&c->muldet);
	circ_res_destroy(c);
}


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
