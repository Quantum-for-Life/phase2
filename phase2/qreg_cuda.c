#include "c23_compat.h"
#include <complex.h>
#include <stdlib.h>

#include <cuComplex.h>
#include <cuda_runtime_api.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#include "qreg_cuda.h"
#include "qreg_impl.h"
#include "world_cuda.h"

typedef _Complex double c64;

int qreg_backend_init(struct qreg *reg)
{
	cuDoubleComplex *damp, *dbuf;

	struct qreg_cuda *cu = malloc(sizeof *cu);
	if (!cu)
		goto err_cu_alloc;
	if (cudaMalloc((void **)&damp, reg->namp * sizeof *damp) != cudaSuccess)
		goto err_cuda_malloc_damp;
	if (cudaMalloc((void **)&dbuf, reg->namp * sizeof *dbuf) != cudaSuccess)
		goto err_cuda_malloc_dbuf;
	cu->damp = damp;
	cu->dbuf = dbuf;
	reg->data = cu;

	return 0;

	cudaFree(dbuf);
err_cuda_malloc_dbuf:
	cudaFree(damp);
err_cuda_malloc_damp:
	free(cu);
err_cu_alloc:

	return -1;
}

void qreg_backend_destroy(struct qreg *reg)
{
	struct qreg_cuda *cu = reg->data;
	cudaFree(cu->dbuf);
	cudaFree(cu->damp);
	free(cu);
}

void qreg_getamp(const struct qreg *reg, const uint64_t i, c64 *z)
{
	struct qreg_cuda *cu = reg->data;

	const uint64_t i_lo = qreg_getilo(reg, i);
	const uint64_t rank = qreg_getihi(reg, i);

	cudaDeviceSynchronize();
	if (reg->wd.rank == (int)rank)
		cudaMemcpy(z, cu->damp + i_lo, sizeof(cuDoubleComplex),
			cudaMemcpyDeviceToHost);

	MPI_Bcast(z, 2, MPI_DOUBLE, rank, MPI_COMM_WORLD);
}

void qreg_setamp(struct qreg *reg, const uint64_t i, c64 z)
{
	struct qreg_cuda *cu = reg->data;

	const uint64_t i_lo = qreg_getilo(reg, i);
	const uint64_t rank = qreg_getihi(reg, i);

	cudaDeviceSynchronize();
	if (reg->wd.rank == (int)rank)
		cudaMemcpy(cu->damp + i_lo, &z, sizeof(cuDoubleComplex),
			cudaMemcpyHostToDevice);

	MPI_Barrier(MPI_COMM_WORLD);
}

void qreg_zero(struct qreg *reg)
{
	struct qreg_cuda *cu = reg->data;

	/* cuDoubleComplex zero representation is all bits set to zero */
	cudaMemset(cu->damp, 0, reg->namp * sizeof(cuDoubleComplex));
	cudaDeviceSynchronize();
}

static void exch_init(struct qreg *reg, const int rnk_rem)
{
	struct qreg_cuda *cu = reg->data;

	const int nr = reg->nreqs;
	for (int i = 0; i < nr; i++) {
		const size_t offset = i * reg->msg_count;

		MPI_Isend(cu->damp + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_snd + i);
		MPI_Irecv(cu->dbuf + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_rcv + i);
	}
}

static void exch_waitall(struct qreg *reg)
{
	const int nr = reg->nreqs;

	MPI_Waitall(nr, reg->reqs_snd, MPI_STATUSES_IGNORE);
	MPI_Waitall(nr, reg->reqs_rcv, MPI_STATUSES_IGNORE);
}

static void qreg_paulirot_hi(struct qreg *reg, struct paulis code_hi, c64 *bm)
{
	paulis_shr(&code_hi, reg->nqb_lo);
	const uint64_t rnk_rem = paulis_effect(code_hi, reg->wd.rank, nullptr);

	exch_init(reg, rnk_rem);
	paulis_effect(code_hi, rnk_rem, bm);
	exch_waitall(reg);
}

/* Defined in qreg_cuda_lo.cu */
void qreg_paulirot_lo(struct qreg *reg, const struct paulis *codes_lo,
	const double *angles, const size_t ncodes, c64 bm);

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *phis,
	const size_t ncodes)
{
	c64 bm = 1.0;
	qreg_paulirot_hi(reg, code_hi, &bm);
	qreg_paulirot_lo(reg, codes_lo, phis, ncodes, bm);
}
