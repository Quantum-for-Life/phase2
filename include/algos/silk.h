#ifndef PHASE2_silk_H
#define PHASE2_silk_H

#include "circ.h"
#include "data2.h"

#define silk_NAME "silk"
#define silk_NUM_MEA_QB (1)
#define silk_NUM_ANC_QB (0)

struct silk_data_multidet {
	size_t num_dets;
	struct {
		long long	index;
		_Complex double coeff;
	} * dets;
};

struct silk_data_times {
	size_t num_steps;
	struct {
		double		t;
		_Complex double val;
	} * steps;
};

struct silk_data {
	struct circ_hamil	   hamil;
	struct silk_data_multidet multidet;
	struct silk_data_times	   times;
};

void silk_data_init(struct silk_data *rd);

void silk_data_destroy(struct silk_data *rd);

int silk_data_from_data(struct silk_data *rd, data2_id fid);

int silk_data_times_write(data2_id fid, struct silk_data_times *ts);

int silk_simulate(const struct silk_data *rd);

#endif // PHASE2_silk_H
