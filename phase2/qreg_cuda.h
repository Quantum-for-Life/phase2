#ifndef QREG_CUDA_H
#define QREG_CUDA_H

#include "c23_compat.h"
#include <stdint.h>

#include <cuda_runtime_api.h>

#ifdef __cplusplus
extern "C" {
#endif

struct qreg_cuda {
	cuDoubleComplex *damp, *dbuf;
};

void qreg_paulirot_lo(struct qreg *reg, const struct paulis *codes_lo,
	const double *angles, const size_t ncodes, _Complex double buf_mul);

#ifdef __cplusplus
}
#endif

#endif /* QREG_CUDA_H */
