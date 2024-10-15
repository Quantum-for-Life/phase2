#include <complex.h>
#include <stdlib.h>

#include <cuComplex.h>
#include <cuda_runtime_api.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "qreg_cuda.h"
#include "world_cuda.h"

typedef _Complex double c64;

#define MAX_COUNT (1 << 29)

/* Local copy of the world info. Initialized by qreg_init() */
static struct world WD;

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
	if (cudaMalloc((void **)&d_sv, num_amps * sizeof *d_sv) != cudaSuccess)
		goto err_cuda_malloc_dsv;
	if (cudaMalloc((void **)&d_buf, num_amps * sizeof *d_sv) != cudaSuccess)
		goto err_cuda_malloc_dbuf;
	cu->d_sv = d_sv;
	cu->d_buf = d_buf;
	for (size_t i = 0; i < qb_lo; i++)
		cu->targs[i] = i;
	cu->num_targs = qb_lo;
	cu->num_qubits = qb_lo;

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
	struct qreg_cuQuantum *cu = reg->data;

	qb_split(i, reg->qb_lo, reg->qb_hi, &loci, &rank);

	cudaDeviceSynchronize();
	if (WD.rank == (int)rank)
		cudaMemcpy(z, cu->d_sv + loci, sizeof(cuDoubleComplex),
			cudaMemcpyDeviceToHost);

	MPI_Bcast(z, 2, MPI_DOUBLE, rank, MPI_COMM_WORLD);
}

void qreg_setamp(struct qreg *reg, const uint64_t i, c64 z)
{
	uint64_t rank, loci;
	struct qreg_cuQuantum *cu = reg->data;

	qb_split(i, reg->qb_lo, reg->qb_hi, &loci, &rank);

	cudaDeviceSynchronize();
	if (WD.rank == (int)rank)
		cudaMemcpy(cu->d_sv + loci, &z, sizeof(cuDoubleComplex),
			cudaMemcpyHostToDevice);

	MPI_Barrier(MPI_COMM_WORLD);
}

void qreg_zero(struct qreg *reg)
{
	struct qreg_cuQuantum *cu = reg->data;

	/* cuDoubleComplex zero representation is all bits set to zero */
	cudaMemset(cu->d_sv, 0, reg->num_amps * sizeof(cuDoubleComplex));
	cudaDeviceSynchronize();
}

static void qreg_exchbuf_init(struct qreg *reg, const int rnk_rem)
{
	struct qreg_cuQuantum *cu = reg->data;

	const int nr = reg->num_reqs;

	for (int i = 0; i < nr; i++) {
		const size_t offset = i * reg->msg_count;

		MPI_Isend(cu->d_sv + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_snd + i);
		MPI_Irecv(cu->d_buf + offset, reg->msg_count * 2, MPI_DOUBLE,
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

	qreg_paulirot_local(reg, codes_lo, angles, num_codes, buf_mul);
}
