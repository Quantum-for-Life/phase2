#ifndef CIRC_CMPSIT_H
#define CIRC_CMPSIT_H

#include <stddef.h>

#include "phase2.h"
#include "xoshiro256ss.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CIRC_CMPSIT_TRUNC_DIST (16UL)

struct circ_cmpsit_data {
	size_t depth;
	size_t length;
	size_t samples;
	double step_size;
	size_t steps;
};

/* Probability distribution to sample from */
struct circ_cmpsit_pd {
	double *x;
	size_t len;
	double lambda_r;
};

/* Sampled circuit. */
struct circ_cmpsit_rct {
	struct circ_hamil_term *trm;
	size_t len;
};

/* Results */
struct circ_cmpsit_samples {
	_Complex double *a;
	size_t len;
};

struct circ_cmpsit {
	struct circ circ;
	struct circ_cmpsit_data data;
	struct circ_cmpsit_pd pd;
	struct circ_cmpsit_rct rct;
	struct circ_cmpsit_samples samples;
	struct xoshiro256ss rng;
};

int circ_cmpsit_init(struct circ_cmpsit *ct,
	const struct circ_cmpsit_data *data, data_id fid);

void circ_cmpsit_destroy(struct circ_cmpsit *ct);

#ifdef __cplusplus
}
#endif

#endif // CIRC_CMPSIT_H
