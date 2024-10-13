#include <complex.h>
#include <stddef.h>

#include <cuda_runtime.h>
#include <cuComplex.h>

#include "custatevec.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "qreg_cuQuantum.h"
#include "world_cuQuantum.h"

const size_t threadPerBlock = 512;

__global__ void kernelAdd(cuDoubleComplex *a, cuDoubleComplex *b, size_t n)
{
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i < n)
		a[i] = cuCadd(a[i], b[i]);

}

__global__ void kernelMul(cuDoubleComplex *a, cuDoubleComplex b, size_t n)
{
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i < n)
		a[i] = cuCmul(a[i], b);
}

__global__ void kernelMix(cuDoubleComplex *a, cuDoubleComplex *b, size_t n)
{
	cuDoubleComplex z1, z2;
	const cuDoubleComplex half = { .x = 0.5, .y = 0.0 };

	size_t i = threadIdx.x + blockIdx.x * blockDim.x;

	if (i < n) {
		z1 = cuCadd(a[i], b[i]);
		z2 = cuCsub(a[i], b[i]);
		a[i] = cuCmul(z1, half);
		b[i] = cuCmul(z2, half);
	}
}



void qreg_paulirot_local(struct qreg *reg, custatevecHandle_t handle,
	       const struct paulis *codes_lo, const double *angles,
	       const size_t num_codes, double _Complex buf_mul)
{
 	const size_t blocks = (reg->num_amps + threadPerBlock - 1) / threadPerBlock;

	const struct qreg_cuQuantum *cu =
		(const struct qreg_cuQuantum *)reg->data;
	custatevecPauli_t paulis[QREG_MAX_WIDTH];

	cudaMemcpy(cu->d_sv, reg->amp, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyHostToDevice);
	cudaMemcpy(cu->d_buf, reg->buf, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyHostToDevice);

	/* Note that we're taking the conjugation of buf_mul. */
	const cuDoubleComplex b = { .x = creal(buf_mul), .y = -cimag(buf_mul) };
	kernelMul<<<blocks, threadPerBlock>>>(cu->d_buf, b, reg->num_amps);
	kernelMix<<<blocks, threadPerBlock>>>(cu->d_sv, cu->d_buf, reg->num_amps);

	for (size_t k = 0; k < num_codes; k++) {
		/* cuQuantum Paulis are the same as ours:
		 * CUSTATEVEC_PAULI_I = PAULI_I = 0, etc.
		 */
		for (size_t i = 0; i < cu->num_qubits; i++)
			paulis[i] = (custatevecPauli_t)paulis_get(codes_lo[k], i);

		// apply exponential
		custatevecApplyPauliRotation(
			handle, cu->d_sv, CUDA_C_64F,
			cu->num_qubits, angles[k], paulis,
			cu->targs, cu->num_targs,
	    		nullptr, nullptr, 0);
		custatevecApplyPauliRotation(
			handle, cu->d_buf, CUDA_C_64F,
			cu->num_qubits, -angles[k], paulis,
			cu->targs, cu->num_targs,
	    		nullptr, nullptr, 0);
	}

        kernelAdd<<<blocks, threadPerBlock>>>(cu->d_sv, cu->d_buf, reg->num_amps);

	cudaMemcpy(reg->amp, cu->d_sv, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyDeviceToHost);
	cudaMemcpy(reg->buf, cu->d_buf, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyDeviceToHost);

}
