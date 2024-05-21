#ifndef PHASE2_CIRC_TROTTER_H
#define PHASE2_CIRC_TROTTER_H

#include <stdint.h>

#include "data.h"

struct circ_trotter_hamil {
	size_t	       num_qubits;
	size_t	       num_terms;
	double	      *coeffs;
	struct paulis *paulis;
};

struct circ_trotter_multidet {
	struct {
		uint64_t idx;
		double	 coeff[2];
	}     *dets;
	size_t num_dets;
};

struct circ_trotter_data {
	struct circ_trotter_hamil    hamil;
	struct circ_trotter_multidet multidet;

	double time_factor;

	double *trott_steps[2];
	size_t	num_trott_steps;
};

int
circ_trotter_data_init(struct circ_trotter_data *cd, size_t num_steps);

void
circ_trotter_data_destroy(struct circ_trotter_data *cd);

int
circ_trotter_data_from_file(struct circ_trotter_data *cd, data_id fid);

int
circ_trotter_simulate(const struct circ_trotter_data *cd);

#endif // PHASE2_CIRC_TROTTER_H
