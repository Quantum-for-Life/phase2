#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "circ.h"

void circ_hamil_destroy(struct circ_hamil *h)
{
	if (h->paulis)
		free(h->paulis);
	if (h->coeffs)
		free(h->coeffs);
}

struct hamil_iter_data {
	size_t	       idx;
	size_t	       num_qubits;
	double	       norm;
	double	      *coeffs;
	struct paulis *paulis;
};

static int hamil_iter(double coeff, unsigned char *paulis, void *iter_data)
{
	struct hamil_iter_data *idat	   = iter_data;
	const size_t		i	   = idat->idx++;
	const size_t		num_qubits = idat->num_qubits;

	idat->coeffs[i] = coeff * idat->norm;
	struct paulis p = paulis_new();
	for (size_t j = 0; j < num_qubits; j++) {
		paulis_set(&p, paulis[j], j);
	}
	idat->paulis[i] = p;

	return 0;
}

int circ_hamil_from_file(struct circ_hamil *h, const data_id fid)
{
	size_t num_qubits, num_terms;
	double norm;

	if (data_hamil_getnums(fid, &num_qubits, &num_terms) < 0)
		return -1;
	if (data_hamil_getnorm(fid, &norm) < 0)
		return -1;

	double	      *coeffs = malloc(sizeof *coeffs * num_terms);
	struct paulis *paulis = malloc(sizeof(struct paulis) * num_terms);
	if (!(coeffs && paulis))
		goto err;

	struct hamil_iter_data idat = { .idx = 0,
		.num_qubits		     = num_qubits,
		.norm			     = norm,
		.coeffs			     = coeffs,
		.paulis			     = paulis };
	if (data_hamil_foreach(fid, hamil_iter, &idat) != 0)
		goto err;

	h->num_qubits = num_qubits;
	h->num_terms  = num_terms;
	h->coeffs     = coeffs;
	h->paulis     = paulis;

	return 0;
err:
	free(coeffs);
	free(paulis);
	return -1;
}

void circ_multidet_destroy(struct circ_multidet *md)
{
	if (md->dets)
		free(md->dets);
}

struct iter_multidet_data {
	size_t		      i;
	struct circ_multidet *md;
};

static int iter_multidet(double coeff[2], const uint64_t idx, void *op_data)
{
	struct iter_multidet_data *imd = op_data;

	imd->md->dets[imd->i].coeff[0] = coeff[0];
	imd->md->dets[imd->i].coeff[1] = coeff[1];
	imd->md->dets[imd->i].idx      = idx;
	imd->i++;

	return 0;
}

int circ_multidet_from_file(struct circ_multidet *md, const data_id fid)
{
	size_t num_qubits, num_dets;
	if (data_multidet_getnums(fid, &num_qubits, &num_dets) < 0)
		return -1;
	md->dets = malloc(sizeof *md->dets * num_dets);
	if (md->dets == NULL)
		return -1;

	struct iter_multidet_data imd;
	imd.i  = 0;
	imd.md = md;
	if (data_multidet_foreach(fid, iter_multidet, &imd) < 0)
		goto error;

	md->num_dets = num_dets;

	return 0;
error:
	free(md->dets);
	return -1;
}
