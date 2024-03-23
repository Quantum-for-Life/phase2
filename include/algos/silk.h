#ifndef PHASE2_SILK_H
#define PHASE2_SILK_H

#include "circ.h"
#include "data2.h"

#define SILK_NAME "silk"
#define SILK_NUM_MEA_QB (0)
#define SILK_NUM_ANC_QB (0)

struct silk_data_multidet {
	size_t num_dets;
	struct {
		long long	index;
		_Complex double coeff;
	} *dets;
};

struct silk_data {
	struct circ_hamil	  hamil;
	struct silk_data_multidet multidet;
	double			  time_factor;
	size_t			  num_steps;
	_Complex double		 *trotter_steps;
};

int silk_data_init(struct silk_data *rd, size_t num_steps);

void silk_data_destroy(struct silk_data *rd);

int silk_data_from_data(struct silk_data *rd, data2_id fid);

int silk_simulate(const struct silk_data *rd);

#endif // PHASE2_silk_H
