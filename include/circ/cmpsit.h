#ifndef CMPSIT_H
#define CMPSIT_H

#include <stddef.h>

#include "phase2.h"
#include "xoshiro256ss.h"

#define CMPSIT_TRUNC_DIST (16UL)

struct cmpsit_data {
	size_t depth;
	size_t length;
	size_t samples;
	double step_size;
	size_t steps;
};

/* Probability distribution to sample from */
struct cmpsit_pd {
	double *x;
	size_t len;
	double lambda_r;
};

/* Sampled circuit. */
struct cmpsit_rct {
	struct circ_hamil_term *trm;
	size_t len;
};

/* Results */
struct cmpsit_samples {
	_Complex double *a;
	size_t len;
};

struct cmpsit {
	struct circ circ;
	struct cmpsit_data data;
	struct cmpsit_pd pd;
	struct cmpsit_rct rct;
	struct cmpsit_samples samples;
	struct xoshiro256ss rng;
};

int cmpsit_init(struct cmpsit *ct, const struct cmpsit_data *data, data_id fid);

void cmpsit_destroy(struct cmpsit *ct);

#endif // CMPSIT_H
