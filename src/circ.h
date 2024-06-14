#ifndef CIRC_H
#define CIRC_H

#include <stddef.h>
#include <stdint.h>

#include "data.h"
#include "qreg.h"

struct circ_hamil {
	size_t	       num_qubits;
	size_t	       num_terms;
	double	      *coeffs;
	struct paulis *paulis;
};

struct circ_multidet {
	struct {
		uint64_t idx;
		double	 coeff[2];
	}     *dets;
	size_t num_dets;
};

void circ_hamil_init(struct circ_hamil *h);
void circ_hamil_destroy(struct circ_hamil *h);
int  circ_hamil_from_file(struct circ_hamil *h, data_id fid);

void circ_multidet_init(struct circ_multidet *md);
void circ_multidet_destroy(struct circ_multidet *md);
int  circ_multidet_from_file(struct circ_multidet *md, data_id fid);

#endif // CIRC_H