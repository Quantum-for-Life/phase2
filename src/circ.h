#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include <stdint.h>

#include "data2.h"

struct circ_hamil {
	size_t	       num_qubits;
	size_t	       num_terms;
	double	      *coeffs;
	struct paulis *paulis;
};

struct circ_multidet {
	struct {
		uint64_t	idx;
		_Complex double coeff;
	}     *dets;
	size_t num_dets;
};

struct circ_data {
	struct circ_hamil    hamil;
	struct circ_multidet multidet;

	double time_factor;

	_Complex double *trott_steps;
	size_t		 num_trott_steps;
};

int circ_data_init(struct circ_data *cd, size_t num_steps);

void circ_data_destroy(struct circ_data *cd);

int circ_data_from_file(struct circ_data *cd, data2_id fid);

int circ_simulate(const struct circ_data *cd);

#endif // PHASE2_CIRC_H
