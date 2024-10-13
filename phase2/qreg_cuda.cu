#include <complex.h>
#include <stddef.h>

#include <cuda_runtime.h>
#include <cuComplex.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "qreg_cuda.h"
#include "world_cuda.h"

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

__global__ void kernelPauliRot(cuDoubleComplex *a, size_t n, cuDoubleComplex eip,
		struct paulis code)
{
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= n)
		return;

	size_t j = i ^ code.pak[0];
	if (j < i)
		return;

	int minus = __popcll(j & code.pak[1]);
	int root4 = (__popcll(code.pak[0] & code.pak[1]) + 2*minus) & 0x3;
	cuDoubleComplex z;
	switch (root4) {
		case 0:
			z.x = 1.0;
			z.y = 0.0;
			break;
		case 1:
			z.x = 0.0;
			z.y = -1.0;
			break;
		case 2:
			z.x = -1.0;
			z.y = 0.0;
			break;
		case 3:
			z.x = 0.0;
			z.y = 1.0;
			break;
		default:
			__builtin_unreachable();
	}

	cuDoubleComplex zi, zj;
	zi = a[i];
	zj = a[j];

	cuDoubleComplex rc = { .x = cuCreal(eip), .y = 0.0 };
	cuDoubleComplex is = { .x = 0.0, .y = cuCimag(eip) };
	a[i] = cuCadd(cuCmul(rc, zi), cuCmul(is, cuCmul(cuConj(z), zj)));
	a[j] = cuCadd(cuCmul(rc, zj), cuCmul(is, cuCmul(z, zi)));
}

void qreg_paulirot_local(struct qreg *reg,
	       const struct paulis *codes_lo, const double *angles,
	       const size_t num_codes, double _Complex buf_mul)
{
 	const size_t blocks = (reg->num_amps + threadPerBlock - 1) / 
				threadPerBlock;

	const struct qreg_cuQuantum *cu =
		(const struct qreg_cuQuantum *)reg->data;

	/* Note that we're taking the conjugation of buf_mul. */
	const cuDoubleComplex b = { .x = creal(buf_mul), .y = -cimag(buf_mul) };
	kernelMul<<<blocks, threadPerBlock>>>(cu->d_buf, b, reg->num_amps);
	kernelMix<<<blocks, threadPerBlock>>>(cu->d_sv, cu->d_buf, reg->num_amps);

	cudaDeviceSynchronize();
	for (size_t k = 0; k < num_codes; k++) {
		cuDoubleComplex eip = {
			.x = cos(angles[k]),
			.y = sin(angles[k])
		};
		kernelPauliRot<<<blocks, threadPerBlock>>>
			(cu->d_sv, reg->num_amps, eip, codes_lo[k]);
		kernelPauliRot<<<blocks, threadPerBlock>>>
			(cu->d_buf, reg->num_amps, cuConj(eip), codes_lo[k]);

	}

	/* We mix again d_sv and d_buf. Sync them first. */
	cudaDeviceSynchronize();
        kernelAdd<<<blocks, threadPerBlock>>>(cu->d_sv, cu->d_buf,
			reg->num_amps);
}
