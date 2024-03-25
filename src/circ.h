/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include "common.h"
#include "data2.h"

#include "qreg.h"

struct circuit_data_multidet {
	size_t num_dets;
	struct {
		long long	index;
		_Complex double coeff;
	} *dets;
};

/*
 * Structure representing a Hamiltonian as a real linear combination of Pauli
 * strings (codes).  The Pauli operators are packed as a continuous string
 * of pairs of bits.
 */
struct circ_hamil {
	size_t	num_qubits; /* number of qubits */
	size_t	num_terms; /* number of terms in the sum */
	double *coeffs; /* array of coefficients */
	u64    *pak; /* array of Pauli operators */
};

struct circuit_data {
	struct circ_hamil	     hamil;
	struct circuit_data_multidet multidet;
	double			     time_factor;
	size_t			     num_steps;
	_Complex double		    *trotter_steps;
};

int circuit_data_init(struct circuit_data *rd, size_t num_steps);

void circuit_data_destroy(struct circuit_data *rd);

int circuit_data_from_data(struct circuit_data *rd, data2_id fid);

int circuit_simulate(const struct circuit_data *rd);

#endif // PHASE2_CIRC_H
