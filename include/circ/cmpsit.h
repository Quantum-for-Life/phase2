#ifndef CIRC_QDRIFT_H
#define CIRC_QDRIFT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "phase2/circ.h"
#include "phase2/data.h"
#include "xoshiro256ss.h"

struct qdrift_data {
	size_t depth;
	double step_size;
	size_t nsamples;
};

struct qdrift {
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

int qdrift_init(
	struct qdrift *qd, struct qdrift_data *data, data_id fid);

void qdrift_destroy(struct qdrift *qd);

#ifdef __cplusplus
}
#endif

#endif // CIRC_QDRIFT_H
