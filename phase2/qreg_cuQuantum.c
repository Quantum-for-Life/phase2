/* qreg implementation using cuQuantum library */
#include <complex.h>
#include <stdlib.h>

#include <cuda_runtime_api.h> // cudaMalloc, cudaMemcpy, etc.
#include <cuComplex.h>        // cuDoubleComplex
#include "custatevec.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "world_cuQuantum.h"

typedef _Complex double c64;

#define MAX_COUNT (1 << 29)

/* Local copy of the world info. Initialized by qreg_init() */
static struct world WD;

struct qreg_cuQuantum {
	uint32_t num_qubits;
	cuDoubleComplex *d_sv, *d_buf;
	int32_t targs[QREG_MAX_WIDTH];
	uint32_t num_targs;
};


int qreg_init(struct qreg *reg, const uint32_t num_qubits)
{
	struct qreg_cuQuantum *cu;
	cuDoubleComplex *d_sv, *d_buf;

	if (world_info(&WD) != WORLD_READY)
		return -1;

	uint32_t qb_hi = 0, nrk = WD.size;
	while (nrk >>= 1)
		qb_hi++;
	if (qb_hi >= num_qubits)
		return -1;
	const uint32_t qb_lo = num_qubits - qb_hi;
	const uint64_t num_amps = UINT64_C(1) << qb_lo;

	const int msg_count = num_amps < MAX_COUNT ? num_amps : MAX_COUNT;
	const size_t num_reqs = num_amps / msg_count;

	MPI_Request *reqs = malloc(sizeof *reqs * num_reqs * 2);
	if (!reqs)
		goto err_reqs_alloc;
	c64 *amp = malloc(sizeof *amp * num_amps * 2);
	if (!amp)
		goto err_amp_alloc;

	if (!(cu = malloc(sizeof *cu)))
		goto err_cu_alloc;
	if (cudaMalloc((void**)&d_sv, num_amps*sizeof *d_sv) != cudaSuccess)
		goto err_cuda_malloc_dsv;
	if (cudaMalloc((void**)&d_buf, num_amps*sizeof *d_sv) != cudaSuccess)
		goto err_cuda_malloc_dbuf;
	cu->d_sv = d_sv;
	cu->d_buf = d_buf;
	for (size_t i = 0; i < num_qubits; i++)
		cu->targs[i] = i;
	cu->num_targs = num_qubits;
	cu->num_qubits = num_qubits;

	reg->qb_lo = qb_lo;
	reg->qb_hi = qb_hi;
	reg->amp = amp;
	reg->buf = amp + num_amps;
	reg->num_amps = num_amps;
	reg->msg_count = msg_count;
	reg->reqs_snd = reqs;
	reg->reqs_rcv = reqs + num_reqs;
	reg->num_reqs = num_reqs;
	reg->data = cu;

	return 0;

	cudaFree(d_buf);
err_cuda_malloc_dbuf:
	cudaFree(d_sv);
err_cuda_malloc_dsv:
	free(cu);
err_cu_alloc:
	free(amp);
err_amp_alloc:
	free(reqs);
err_reqs_alloc:
	return -1;
}

void qreg_destroy(struct qreg *reg)
{
	struct qreg_cuQuantum *cu = reg->data;
	cudaFree(cu->d_buf);
	cudaFree(cu->d_sv);
	free(cu);

	if (reg->amp)
		free(reg->amp);
	reg->amp = NULL;
	reg->buf = NULL;

	if (reg->reqs_snd)
		free(reg->reqs_snd);
	reg->reqs_snd = NULL;
	reg->reqs_rcv = NULL;
}

static void qb_split(uint64_t n, const uint32_t qb_lo, const uint32_t qb_hi,
	uint64_t *lo, uint64_t *hi)
{
	const uint64_t mask_lo = (UINT64_C(1) << qb_lo) - 1;
	const uint64_t mask_hi = (UINT64_C(1) << qb_hi) - 1;

	*lo = n & mask_lo;
	n >>= qb_lo;
	*hi = n & mask_hi;
}

void qreg_getamp(const struct qreg *reg, const uint64_t i, c64 *z)
{
	uint64_t rank, loci;
	qb_split(i, reg->qb_lo, reg->qb_hi, &loci, &rank);

	if (WD.rank == (int)rank)
		*z = reg->amp[loci];
	MPI_Bcast(z, 2, MPI_DOUBLE, rank, MPI_COMM_WORLD);
}

void qreg_setamp(struct qreg *reg, const uint64_t i, c64 z)
{
	uint64_t rank, loci;
	qb_split(i, reg->qb_lo, reg->qb_hi, &loci, &rank);

	if (WD.rank == (int)rank)
		reg->amp[loci] = z;
	MPI_Barrier(MPI_COMM_WORLD);
}

void qreg_zero(struct qreg *reg)
{
	c64 *z = reg->amp;
	while (z < reg->amp + reg->num_amps)
		*z++ = 0.0;
}


static void qreg_exchbuf_init(struct qreg *reg, const int rnk_rem)
{
	const int nr = reg->num_reqs;

	for (int i = 0; i < nr; i++) {
		const size_t offset = i * reg->msg_count;

		MPI_Isend(reg->amp + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_snd + i);
		MPI_Irecv(reg->buf + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_rcv + i);
	}
}

static void qreg_exchbuf_waitall(struct qreg *reg)
{
	const int nr = reg->num_reqs;

	MPI_Waitall(nr, reg->reqs_snd, MPI_STATUSES_IGNORE);
	MPI_Waitall(nr, reg->reqs_rcv, MPI_STATUSES_IGNORE);
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t num_codes)
{
	/* Compute permutation from outer qubits */
	struct paulis hi = code_hi;
	paulis_shr(&hi, reg->qb_lo);

	const uint64_t rnk_loc = WD.rank;
	const uint64_t rnk_rem = paulis_effect(hi, rnk_loc, NULL);
	qreg_exchbuf_init(reg, rnk_rem);

	/* Compute multiplication factor for the buffer
	   code_hi acts on the value of rank_remote now
	   (as if from receiving end). We discard the result. */
	c64 buf_mul = 1.0;
	paulis_effect(hi, rnk_loc, &buf_mul);

	qreg_exchbuf_waitall(reg);

	/* Compute permutation from inner qubits */
	for (uint64_t i = 0; i < reg->num_amps; i++) {
		reg->buf[i] *= conj(buf_mul);

		_Complex a = reg->amp[i], b = reg->buf[i];
		reg->amp[i] = (a + b) / 2.0;
		reg->buf[i] = (a - b) / 2.0;
	}


	const struct world_cuQuantum *w = WD.data;
	const struct qreg_cuQuantum *cu = reg->data;
	custatevecPauli_t paulis[QREG_MAX_WIDTH];

	cudaMemcpy(cu->d_sv, reg->amp, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyHostToDevice);
	cudaMemcpy(cu->d_buf, reg->buf, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyHostToDevice);

	for (size_t k = 0; k < num_codes; k++) {
		/* cuQuantum Paulis are the same as ours:
		 * CUSTATEVEC_PAULI_I = PAULI_I = 0, etc.
		 */
		for (size_t i = 0; i < cu->num_qubits; i++)
			paulis[i] = (int)paulis_get(codes_lo[k], i);

		// apply exponential
		custatevecApplyPauliRotation(
			w->handle, cu->d_sv, CUDA_C_64F,
			cu->num_qubits, angles[k], paulis,
			cu->targs, cu->num_targs,
	    		(void *)0, (void *)0, 0);
		custatevecApplyPauliRotation(
			w->handle, cu->d_buf, CUDA_C_64F,
			cu->num_qubits, -angles[k], paulis,
			cu->targs, cu->num_targs,
	    		(void *)0, (void *)0, 0);
	}

	cudaMemcpy(reg->amp, cu->d_sv, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyDeviceToHost);
	cudaMemcpy(reg->buf, cu->d_buf, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyDeviceToHost);

	for (uint64_t i = 0; i < reg->num_amps; i++)
		reg->amp[i] += reg->buf[i];
}
