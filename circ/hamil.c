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
circ_hamil_from_paulis(struct circ_hamil *h, size_t num_qubits,
	size_t num_terms, unsigned char *paulis, double *coeffs)
{
	double	    *hamil_coeffs = malloc(sizeof(*hamil_coeffs) * num_terms);
	pauli_pak_t *pak = calloc(num_terms * num_qubits / PAULI_PAK_SIZE + 1,
		sizeof(pauli_pak_t));
	if (!(hamil_coeffs && pak)) {
		free(coeffs);
		free(pak);
		return -1;
	}

	for (size_t i = 0; i < num_terms; i++) {
		hamil_coeffs[i] = coeffs[i];
		for (size_t j = 0; j < num_qubits; j++) {
			const size_t	  pauli_idx = i * num_qubits + j;
			const pauli_pak_t pauli	    = paulis[pauli_idx];
			ldiv_t		  dv = ldiv(pauli_idx, PAULI_PAK_SIZE);
			pak[dv.quot] += pauli << (dv.rem * PAULI_WIDTH);
		}
	}
	h->num_qubits = num_qubits;
	h->num_terms  = num_terms;
	h->coeffs     = hamil_coeffs;
	h->pak	      = pak;

	return 0;
}

int
circ_hamil_from_data(struct circ_hamil *h, const struct data_pauli_hamil *dat)
{
	int rc;

	double *norm_coeffs = malloc(sizeof(*norm_coeffs) * dat->num_terms);
	if (!norm_coeffs)
		return -1;
	for (size_t i = 0; i < dat->num_terms; i++) {
		norm_coeffs[i] = dat->coeffs[i] * dat->norm;
	}
	rc = circ_hamil_from_paulis(
		h, dat->num_qubits, dat->num_terms, dat->paulis, norm_coeffs);
	free(norm_coeffs);

	return rc;
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