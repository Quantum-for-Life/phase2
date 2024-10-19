#ifndef CIRC_QDRIFT_H
#define CIRC_QDRIFT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct circ_qdrift_data {
	size_t depth;
};

struct circ_drift_res {
	double *samples[2];
	size_t nsamples;
};

#ifdef __cplusplus
}
#endif

#endif // CIRC_QDRIFT_H
