/*
 * CUDA kernels for the lo-part of qreg_paulirot.
 * One thread per amplitude; grid sized to ceil(namp /
 * 512).  kernelPauliRot mirrors the CPU kernel_rot
 * with the same j < i pairing guard, so each coupled
 * pair is touched once and no atomics are needed.
 * Launches go into the default stream; implicit
 * ordering keeps mix -> rotate -> add correct.
 */

#include <complex.h>
#include <stddef.h>

#include <cuComplex.h>
#include <cuda_runtime.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "qreg_cuda.h"

constexpr size_t threadPerBlock = 512;

constexpr cuDoubleComplex half = { .x = 0.5, .y = 0.0 };

/*
 * kernelMix - GPU equivalent of kernel_mix (see qreg_qreg.c).
 *
 * One thread per amplitude.  Grid geometry:
 *   blocks = ceil(namp / 512),  threads per block = 512.
 * Each thread applies the eigenspace projection for its
 * index independently (no inter-thread data dependency).
 */
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

/*
 * kernelPauliRot - GPU equivalent of kernel_rot (qreg_qreg.c).
 *
 * One thread per amplitude index.  Same j < i guard: each
 * coupled pair (i, j) is processed by the thread with the
 * smaller index only, avoiding write conflicts without
 * atomics or synchronisation.
 *
 * The popcount + phase-exponent computation comes from
 * paulis_effect_raw (paulis.h), which is shared with the
 * host paulis_effect implementation.  Only the
 * cuDoubleComplex multiply lives here.
 */
__global__ void kernelPauliRot(
	cuDoubleComplex *a, size_t n, struct paulis code, double c, double s)
{
	size_t i = threadIdx.x + blockIdx.x * blockDim.x;
	if (i >= n)
		return;

	int r4;
	const uint64_t j = paulis_effect_raw(code, i, &r4);
	if (j < i)
		return;

	const cuDoubleComplex phase[4] = {
		{ .x = 1.0, .y = 0.0 },
		{ .x = 0.0, .y = 1.0 },
		{ .x = -1.0, .y = 0.0 },
		{ .x = 0.0, .y = -1.0 },
	};
	const cuDoubleComplex z = phase[r4];

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

/*
 * qreg_paulirot_lo - host function launching the lo-qubit
 * rotation kernels on the GPU.
 *
 * All kernels are launched sequentially into the default
 * CUDA stream.  Implicit stream ordering guarantees that
 * kernelMix completes before kernelPauliRot, and all
 * rotations complete before kernelAdd.  No explicit
 * synchronisation is needed between launches.
 */
void qreg_paulirot_lo(struct qreg *reg, const struct paulis *codes_lo,
	const double *angles, const size_t ncodes, double _Complex bm)
{
	const size_t blocks = (reg->namp + threadPerBlock - 1) / threadPerBlock;
	const struct qreg_cuda *cu = (const struct qreg_cuda *)reg->backend;

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
	kernelAdd<<<blocks, threadPerBlock>>>(cu->damp, cu->dbuf, reg->namp);
}
