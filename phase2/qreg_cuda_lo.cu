#include <complex.h>
#include <stddef.h>

#include <cuComplex.h>
#include <cuda_runtime.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "qreg_cuda.h"

constexpr size_t threadPerBlock = 512;

constexpr cuDoubleComplex half = { .x = 0.5, .y = 0.0 };

__global__ void kernelMix(cuDoubleComplex *__restrict__ a,
	cuDoubleComplex *__restrict__ b, cuDoubleComplex bm, size_t n)
{
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= n)
		return;

	b[i] = cuCmul(b[i], bm);

	const cuDoubleComplex z1 = cuCadd(a[i], b[i]);
	const cuDoubleComplex z2 = cuCsub(a[i], b[i]);
	a[i] = cuCmul(z1, half);
	b[i] = cuCmul(z2, half);
}

// Assume z != nullptr
__device__ uint64_t paulisEffect(
	const struct paulis code, const uint64_t i, cuDoubleComplex *z)
{
	int mi = __popcll(i & code.pak[1]); // no. of minuses
	int is = __popcll(code.pak[0] & code.pak[1]); // no. of i's
	int r4 = (is + 2 * mi) & 0x3; // 4th root of unity
	switch (r4) {
	case 0:
		break;
	case 1:
		*z = cuCmul(*z, (cuDoubleComplex){ .x = 0.0, .y = 1.0 });
		break;
	case 2:
		*z = cuCmul(*z, (cuDoubleComplex){ .x = -1.0, .y = 0.0 });
		break;
	case 3:
		*z = cuCmul(*z, (cuDoubleComplex){ .x = 0.0, .y = -1.0 });
		break;
	default:
		unreachable();
	}

	return i ^ code.pak[0];
}

__global__ void kernelPauliRot(
	cuDoubleComplex *a, size_t n, struct paulis code, double c, double s)
{
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= n)
		return;

	cuDoubleComplex z = { .x = 1.0, .y = 0.0 };
	const uint64_t j = paulisEffect(code, i, &z);
	if (j < i)
		return;

	const cuDoubleComplex zi = a[i];
	const cuDoubleComplex zj = a[j];
	const cuDoubleComplex rc = { .x = c, .y = 0.0 };
	const cuDoubleComplex is = { .x = 0.0, .y = s };

	a[i] = cuCadd(cuCmul(rc, zi), cuCmul(is, cuCmul(cuConj(z), zj)));
	a[j] = cuCadd(cuCmul(rc, zj), cuCmul(is, cuCmul(z, zi)));
}

__global__ void kernelAdd(cuDoubleComplex *__restrict__ a,
	cuDoubleComplex *__restrict__ b, size_t n)
{
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= n)
		return;

	a[i] = cuCadd(a[i], b[i]);
}

void qreg_paulirot_lo(struct qreg *reg, const struct paulis *codes_lo,
	const double *angles, const size_t ncodes, double _Complex bm)
{
	const size_t blocks = (reg->namp + threadPerBlock - 1) / threadPerBlock;
	const struct qreg_cuda *cu = (const struct qreg_cuda *)reg->data;

	const cuDoubleComplex z = { .x = creal(bm), .y = cimag(bm) };
	kernelMix<<<blocks, threadPerBlock>>>(cu->damp, cu->dbuf, z, reg->namp);

	for (size_t k = 0; k < ncodes; k++) {
		double c = cos(angles[k]), s = sin(angles[k]);
		kernelPauliRot<<<blocks, threadPerBlock>>>(
			cu->damp, reg->namp, codes_lo[k], c, s);
		kernelPauliRot<<<blocks, threadPerBlock>>>(
			cu->dbuf, reg->namp, codes_lo[k], c, -s);
	}

	// We mix again damp and dbuf. Sync them first.
	cudaDeviceSynchronize();
	kernelAdd<<<blocks, threadPerBlock>>>(cu->damp, cu->dbuf, reg->namp);
}
