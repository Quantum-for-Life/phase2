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

/* qreg_backend_paulirot_lo lives in qreg_cuda_lo.cu;
 * declared in phase2/qreg.h for backend-neutral dispatch. */

#ifdef __cplusplus
}
#endif

#endif /* QREG_CUDA_H */
