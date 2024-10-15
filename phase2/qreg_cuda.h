#ifndef QREG_CUDA_H
#define QREG_CUDA_H

#include <stdint.h>

#include <cuda_runtime_api.h>

#ifdef __cplusplus
extern "C" {
#endif

struct qreg_cuQuantum {
	uint32_t num_qubits;
	cuDoubleComplex *d_sv, *d_buf;
	int32_t targs[QREG_MAX_WIDTH];
	uint32_t num_targs;
};

void qreg_paulirot_local(struct qreg *reg, const struct paulis *codes_lo,
	const double *angles, const size_t ncodes, _Complex double buf_mul);

#ifdef __cplusplus
}
#endif

#endif /* QREG_CUDA_H */
