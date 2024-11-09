#ifndef CIRC_QDRIFT_H
#define CIRC_QDRIFT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "phase2/circ.h"
#include "phase2/data.h"
#include "xoshiro256ss.h"

struct circ_qdrift_data {
	size_t depth;
	double step_size;
	size_t nsamples;
};

struct circ_qdrift {
	struct circ circ;

	size_t depth;
	double step_size;

	struct xoshiro256ss rng;
	size_t *smpl;

	struct {
		_Complex double *samples;
		size_t nsamples;
	} res;
};

int circ_qdrift_init(
	struct circ_qdrift *qd, struct circ_qdrift_data *data, data_id fid);

void circ_qdrift_destroy(struct circ_qdrift *qd);

#ifdef __cplusplus
}
#endif

#endif // CIRC_QDRIFT_H
