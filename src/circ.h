/**
 * Circuit interface.
 */
#ifndef PHASE2_CIRC_H
#define PHASE2_CIRC_H

#include "common.h"
#include "data2.h"

#include "qreg.h"

struct circuit_multidet {
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

struct circuit_data {
	struct circ_hamil	hamil;
	struct circuit_multidet multidet;
	double			time_factor;
	size_t			num_steps;
	_Complex double	       *trotter_steps;
};

int circuit_data_init(struct circuit_data *rd, size_t num_steps);

void circuit_data_destroy(struct circuit_data *rd);

int circuit_data_from_file(struct circuit_data *rd, data2_id fid);

int circuit_simulate(const struct circuit_data *rd);

#endif // PHASE2_CIRC_H
