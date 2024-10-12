/* qreg implementation using cuQuantum library */
#include <complex.h>
#include <stdlib.h>

#include "custatevec.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "world_cuQuantum.h"

static struct world WD;

struct qreg_cuQuantum {};

int qreg_init(struct qreg *reg, const uint32_t num_qubits)
{
	return 0;
}

void qreg_destroy(struct qreg *reg)
{

}

void qreg_getamp(const struct qreg *reg, const uint64_t i, _Complex double *z)
{
}

void qreg_setamp(struct qreg *reg, const uint64_t i, _Complex double z)
{
}

void qreg_zero(struct qreg *reg)
{
}


void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t num_codes)
{
	struct paulis code;

	for (size_t k = 0; k < num_codes; k++) {
		paulis_merge(&code, reg->qb_lo, reg->qb_hi,
			codes_lo[k], code_hi);
		
		// (...)
	}
}
