#ifndef CIRC_QDRIFT_H
#define CIRC_QDRIFT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct circ_qdrift_data {
	size_t depth;
	double step_size;
	size_t nsamples;
};

struct circ_qdrift_res {
	size_t nsamples;
	_Complex double *samples;
};

#ifdef __cplusplus
}
#endif

#endif // CIRC_QDRIFT_H
