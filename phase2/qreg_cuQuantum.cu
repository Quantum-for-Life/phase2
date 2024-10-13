#include <complex.h>

#include <cuda_runtime.h>
#include <cuComplex.h>

#include "custatevec.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "qreg_cuQuantum.h"
#include "world_cuQuantum.h"

void qreg_paulirot_local(struct qreg *reg, custatevecHandle_t handle,
	       const struct paulis *codes_lo, const double *angles,
	       const size_t num_codes, double _Complex buf_mul)
{
	/* Compute permutation from inner qubits */
	for (uint64_t i = 0; i < reg->num_amps; i++) {
		reg->buf[i] *= conj(buf_mul);

		_Complex double a = reg->amp[i], b = reg->buf[i];
		reg->amp[i] = (a + b) / 2.0;
		reg->buf[i] = (a - b) / 2.0;
	}

	const struct qreg_cuQuantum *cu =
		(const struct qreg_cuQuantum *)reg->data;
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

	cudaMemcpy(reg->amp, cu->d_sv, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyDeviceToHost);
	cudaMemcpy(reg->buf, cu->d_buf, sizeof(double) * 2 * reg->num_amps,
			cudaMemcpyDeviceToHost);

	for (uint64_t i = 0; i < reg->num_amps; i++)
		reg->amp[i] += reg->buf[i];

}
