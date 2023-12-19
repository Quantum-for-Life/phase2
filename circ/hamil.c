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

int
circ_hamil_from_data(struct circ_hamil *h, const struct data_pauli_hamil *dat)
{
	double	    *coeffs = malloc(sizeof(*coeffs) * dat->num_terms);
	pauli_pak_t *pak =
		calloc(dat->num_terms * dat->num_qubits / PAULI_PAK_SIZE + 1,
			sizeof(pauli_pak_t));
	if (!(coeffs && pak)) {
		free(coeffs);
		free(pak);
		return -1;
	}

	for (size_t i = 0; i < dat->num_terms; i++) {
		coeffs[i] = dat->coeffs[i] * dat->norm;
		for (size_t j = 0; j < dat->num_qubits; j++) {
			const size_t pauli_idx = i * dat->num_qubits + j;
			const size_t pak_idx   = pauli_idx / PAULI_PAK_SIZE;
			const size_t pak_offset =
				PAULI_WIDTH * (pauli_idx % PAULI_PAK_SIZE);
			const pauli_pak_t pauli = dat->paulis[pauli_idx];
			pak[pak_idx] += pauli << pak_offset;
		}
	}
	h->num_qubits = dat->num_qubits;
	h->num_terms  = dat->num_terms;
	h->coeffs     = coeffs;
	h->pak	      = pak;

	return 0;
}

void
circ_hamil_paulis(const struct circ_hamil *h, size_t n, int *paulis)
{
	for (size_t j = 0; j < h->num_qubits; j++) {
		const size_t pauli_idx = n * h->num_qubits + j;
		const size_t pak_idx   = pauli_idx / PAULI_PAK_SIZE;
		const size_t pak_offset =
			PAULI_WIDTH * (pauli_idx % PAULI_PAK_SIZE);
		paulis[j] = (h->pak[pak_idx] >> pak_offset) & PAULI_MASK;
	}
}