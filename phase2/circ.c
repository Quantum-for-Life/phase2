#include "c23_compat.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/paulis.h"

int circ_hamil_init(struct circ_hamil *h, size_t nterms)
{
	double *cfs = malloc(sizeof *cfs * nterms);
	if (cfs == nullptr)
		goto err_malloc_cfs;
	struct paulis *ops = malloc(sizeof *ops * nterms);
	if (ops == nullptr)
		goto err_malloc_ops;

	h->nterms = nterms;
	h->cfs = cfs;
	h->ops = ops;

	return 0;

	free(ops);
err_malloc_ops:
	free(cfs);
err_malloc_cfs:

	return -1;
}

void circ_hamil_destroy(struct circ_hamil *h)
{
	if (h->ops != nullptr)
		free(h->ops);
	if (h->cfs != nullptr)
		free(h->cfs);
}

struct hamil_iter_data {
	size_t idx;
	size_t nqb;
	double norm;
	double *cfs;
	struct paulis *ops;
};

static int hamil_iter(double cf, unsigned char *ops, void *iter_data)
{
	struct hamil_iter_data *id = iter_data;
	const size_t i = id->idx++;
	const size_t nqb = id->nqb;

	id->cfs[i] = cf * id->norm;
	struct paulis p = paulis_new();
	for (size_t j = 0; j < nqb; j++) {
		paulis_set(&p, ops[j], j);
	}
	id->ops[i] = p;

	return 0;
}

int circ_hamil_init_from_file(struct circ_hamil *h, const data_id fid)
{
	size_t nqb, nterms;
	double norm;

	if (data_hamil_getnums(fid, &nqb, &nterms) < 0)
		return -1;
	if (data_hamil_getnorm(fid, &norm) < 0)
		return -1;

	if (circ_hamil_init(h, nterms) < 0)
		return -1;

	struct hamil_iter_data id = {
		.idx = 0, .nqb = nqb, .norm = norm, .cfs = h->cfs, .ops = h->ops
	};
	if (data_hamil_foreach(fid, hamil_iter, &id) != 0)
		return -1;

	h->nqb = nqb;

	return 0;
}

int circ_multidet_init(struct circ_multidet *md, size_t ndets)
{
	md->dets = malloc(sizeof *md->dets * ndets);
	if (md->dets == nullptr)
		return -1;

	md->ndets = ndets;

	return 0;
}

void circ_multidet_destroy(struct circ_multidet *md)
{
	if (md->dets != nullptr)
		free(md->dets);
}

struct iter_multidet_data {
	size_t i;
	struct circ_multidet *md;
};

static int iter_multidet(double cf[2], const uint64_t idx, void *iter_data)
{
	struct iter_multidet_data *id = iter_data;

	id->md->dets[id->i].cf[0] = cf[0];
	id->md->dets[id->i].cf[1] = cf[1];
	id->md->dets[id->i].idx = idx;
	id->i++;

	return 0;
}

int circ_multidet_init_from_file(struct circ_multidet *md, const data_id fid)
{
	size_t nqb, ndets;
	if (data_multidet_getnums(fid, &nqb, &ndets) < 0)
		return -1;

	if (circ_multidet_init(md, ndets) < 0)
		return -1;

	struct iter_multidet_data id = { .i = 0, .md = md };
	if (data_multidet_foreach(fid, iter_multidet, &id) < 0)
		return -1;

	return 0;
}
