#include "circ.h"
#include "data.h"

static const size_t PAULI_MASK	   = 3;
static const size_t PAULI_WIDTH	   = 2;
static const size_t PAULI_PAK_SIZE = sizeof(pauli_pak_t) * 8 / PAULI_WIDTH;

void
circ_hamil_init(struct circ_hamil *h)
{
	h->num_qubits = 0;
	h->num_terms  = 0;
	h->coeffs     = NULL;
	h->pak	      = NULL;
}

void
circ_hamil_destroy(struct circ_hamil *h)
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

static int
hamil_iter(double coeff, unsigned char *paulis, void *iter_data)
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

int
circ_hamil_from_data2(struct circ_hamil *h, data2_id fid)
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

void
circ_hamil_paulistr(const struct circ_hamil *h, size_t n, int *paulis)
{
	for (size_t j = 0; j < h->num_qubits; j++) {
		const size_t pauli_idx = n * h->num_qubits + j;
		ldiv_t	     dv	       = ldiv(pauli_idx, PAULI_PAK_SIZE);
		paulis[j] = (h->pak[dv.quot] >> (dv.rem * PAULI_WIDTH)) &
			    PAULI_MASK;
	}
}