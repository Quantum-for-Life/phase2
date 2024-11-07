#ifndef CIRC_QDRIFT_H
#define CIRC_QDRIFT_H

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

struct circ_cmpsit {
	struct circ circ;

	size_t depth;
	size_t length;
	double step_size;
	size_t steps;

	struct xoshiro256ss rng;

	/* Sampled circuit. */
	struct {
		double cf;
		struct paulis op;
	} *smpl_ct;
	size_t nsmpl_ct;

	struct {
		_Complex double *samples;
		size_t nsamples;
	} res;
};

int circ_cmpsit_init(
	struct circ_cmpsit *ct, struct circ_cmpsit_data *data, data_id fid);

void circ_cmpsit_destroy(struct circ_cmpsit *ct);

#ifdef __cplusplus
}
#endif

#endif // CIRC_QDRIFT_H
