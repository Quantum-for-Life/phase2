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

static struct world WD;

struct qreg_cuQuantum {
	uint32_t num_qubits;
	cuDoubleComplex *d_sv;
	int32_t targs[QREG_MAX_WIDTH];
	uint32_t num_targs;
};

int qreg_init(struct qreg *reg, const uint32_t num_qubits)
{
	const size_t num_amps = 1UL << num_qubits;
	struct qreg_cuQuantum *cu;
	cuDoubleComplex *d_sv;
	
	if (world_info(&WD) != WORLD_READY)
		return -1;
	uint32_t qb_hi = 0, nrk = WD.size;
	while (nrk >>= 1)
		qb_hi++;
	if (qb_hi >= num_qubits)
		return -1;
	const uint32_t qb_lo = num_qubits - qb_hi;

	if (!(cu = malloc(sizeof *cu)))
		return -1;

	if (cudaMalloc((void**)&d_sv, num_amps * sizeof(cuDoubleComplex)) !=
		cudaSuccess)
		goto err_cuda_malloc;

	cu->num_qubits = num_qubits;	
	cu->d_sv = d_sv;
	for (size_t i = 0; i < num_qubits; i++) 
		cu->targs[i] = i;
	cu->num_targs = num_qubits;

	reg->qb_lo = qb_lo;
	reg->qb_hi = qb_hi;
	reg->num_amps = num_amps;
	reg->data = cu;

	return 0;

err_cuda_malloc:
	free(cu);
	return -1;
}

void qreg_destroy(struct qreg *reg)
{
	const struct qreg_cuQuantum *cu = reg->data;
	cudaFree(cu->d_sv);
	free(cu);
}

void qreg_getamp(const struct qreg *reg, const uint64_t i, _Complex double *z)
{
	const struct world_cuQuantum *w = WD.data;
	const struct qreg_cuQuantum *cu = reg->data;
	
	cudaDeviceSynchronize();
	cudaMemcpy(z, &(cu->d_sv[i]), sizeof *z, cudaMemcpyDeviceToHost);
}

void qreg_setamp(struct qreg *reg, const uint64_t i, _Complex double z)
{
	const struct world_cuQuantum *w = WD.data;
	const struct qreg_cuQuantum *cu = reg->data;
	
	cudaDeviceSynchronize();
	cudaMemcpy(&(cu->d_sv[i]), &z, sizeof z, cudaMemcpyHostToDevice);
}

void qreg_zero(struct qreg *reg)
{
	const struct qreg_cuQuantum *cu = reg->data;

	/* Double representation of +0.0 is all bits set to 0 */
	/* This function is synchronous. */
	cudaMemset(cu->d_sv, 0, reg->num_amps * sizeof(cuDoubleComplex));
	cudaDeviceSynchronize();
}


void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t num_codes)
{
	const struct world_cuQuantum *w = WD.data;
	const struct qreg_cuQuantum *cu = reg->data;
	custatevecPauli_t paulis[QREG_MAX_WIDTH];

	for (size_t k = 0; k < num_codes; k++) {
		struct paulis code;
		paulis_merge(&code, reg->qb_lo, reg->qb_hi,
			codes_lo[k], code_hi);
		/* cuQuantum Paulis are the same as ours:
		 * CUSTATEVEC_PAULI_I = PAULI_I = 0, etc.
		 */
		for (size_t i = 0; i < cu->num_qubits; i++)
			paulis[i] = paulis_get(code, i);
	
		// apply exponential
		custatevecApplyPauliRotation(
			w->handle, cu->d_sv, CUDA_C_64F,
			cu->num_qubits, angles[k], paulis,
			cu->targs, cu->num_targs,
	    		(void *)0, (void *)0, 0);
	}
}
