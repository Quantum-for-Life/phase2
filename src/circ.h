/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include "common.h"
#include "data2.h"

struct circ_multidet {
	size_t num_dets;
	struct {
		long long	index;
		_Complex double coeff;
	} *dets;
};

struct circ_hamil {
	size_t	num_qubits;
	size_t	num_terms;
	double *coeffs;
	u64    *pak;
};

struct circ_data {
	struct circ_hamil    hamil;
	struct circ_multidet multidet;
	double		     time_factor;
	size_t		     num_steps;
	_Complex double	    *trotter_steps;
};

int circ_data_init(struct circ_data *cd, size_t num_steps);

void circ_data_destroy(struct circ_data *cd);

int circ_data_from_file(struct circ_data *cd, data2_id fid);

int circ_simulate(const struct circ_data *cd);

#endif // PHASE2_CIRC_H
