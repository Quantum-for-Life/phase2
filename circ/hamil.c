#include "circ.h"
#include "data.h"

static const size_t PAULI_MASK	   = 3;
static const size_t PAULI_WIDTH	   = 2;
static const size_t PAULI_PAK_SIZE = sizeof(pauli_pak_t) * 8 / PAULI_WIDTH;

void
circ_hamil_init(struct circ_hamil *hamil)
{
	hamil->num_qubits = 0;
	hamil->num_terms  = 0;
	hamil->coeffs	  = NULL;
	hamil->pak	  = NULL;
}

void
circ_hamil_destroy(struct circ_hamil *hamil)
{
	if (hamil->pak) {
		free(hamil->pak);
		hamil->pak = NULL;
	}
	if (hamil->coeffs) {
		free(hamil->coeffs);
		hamil->coeffs = NULL;
	}
	hamil->num_terms  = 0;
	hamil->num_qubits = 0;
}

int
circ_hamil_from_data(
	struct circ_hamil *hamil, const struct data_pauli_hamil *dat_ph)
{
	double	    *coeffs = malloc(sizeof(*coeffs) * dat_ph->num_terms);
	pauli_pak_t *pak    = calloc(
		   dat_ph->num_terms * dat_ph->num_qubits / PAULI_PAK_SIZE + 1,
		   sizeof(pauli_pak_t));
	if (!(coeffs && pak)) {
		free(coeffs);
		free(pak);
		return -1;
	}

	for (size_t i = 0; i < dat_ph->num_terms; i++) {
		coeffs[i] = dat_ph->coeffs[i] * dat_ph->norm;
		for (size_t j = 0; j < dat_ph->num_qubits; j++) {
			const size_t pauli_idx = i * dat_ph->num_qubits + j;
			const size_t pak_idx   = pauli_idx / PAULI_PAK_SIZE;
			const size_t pak_offset =
				PAULI_WIDTH * (pauli_idx % PAULI_PAK_SIZE);
			const pauli_pak_t pauli = dat_ph->paulis[pauli_idx];
			pak[pak_idx] += pauli << pak_offset;
		}
	}
	hamil->num_qubits = dat_ph->num_qubits;
	hamil->num_terms  = dat_ph->num_terms;
	hamil->coeffs	  = coeffs;
	hamil->pak	  = pak;

	return 0;
}

void
circ_hamil_paulis(const struct circ_hamil *hamil, size_t term, int *paulis)
{
	for (size_t j = 0; j < hamil->num_qubits; j++) {
		const size_t pauli_idx = term * hamil->num_qubits + j;
		const size_t pak_idx   = pauli_idx / PAULI_PAK_SIZE;
		const size_t pak_offset =
			PAULI_WIDTH * (pauli_idx % PAULI_PAK_SIZE);
		paulis[j] = (hamil->pak[pak_idx] >> pak_offset) & PAULI_MASK;
	}
}