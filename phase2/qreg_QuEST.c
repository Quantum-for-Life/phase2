/* qreg implementation using QuEST library */
#include <complex.h>
#include <stdlib.h>

#include "QuEST.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "world_QuEST.h"

static struct world WD;

struct qreg_QuEST {
	size_t	num_qubits;
	Qureg	reg;
	int tar_qb[QREG_MAX_WIDTH];
	pauli_op_t tar_op[QREG_MAX_WIDTH];
};

int qreg_init(struct qreg *reg, const uint32_t num_qubits)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;
	uint32_t qb_hi = 0, nrk = WD.size;
	while (nrk >>= 1)
		qb_hi++;
	if (qb_hi >= num_qubits)
		return -1;
	const uint32_t qb_lo = num_qubits - qb_hi;
	const uint64_t num_amps = UINT64_C(1) << qb_lo;
	
	struct qreg_QuEST *qr = malloc(sizeof *qr);
	if (!reg)
		return -1; 

	struct world_QuEST *qe = WD.data;
	qr->reg = createQureg(num_qubits, qe->env);
	qr->num_qubits = num_qubits;
	for (size_t i = 0; i < num_qubits; i++) {
		qr->tar_qb[i] = i;
		qr->tar_op[i] = PAULI_I;
	}

	reg->num_amps = num_amps; 
	reg->qb_lo = qb_lo;
	reg->qb_hi = qb_hi;
	reg->data = qr;

	return 0;
}

void qreg_destroy(struct qreg *reg)
{
	struct world_QuEST *qe = WD.data;
	struct qreg_QuEST *qr = reg->data;
	
	destroyQureg(qr->reg, qe->env);
}

void qreg_getamp(const struct qreg *reg, const uint64_t i, _Complex double *z)
{
	struct qreg_QuEST *qr = reg->data;
	Complex cz = getAmp(qr->reg, i);
	*z = cz.real + I * cz.imag;
}

void qreg_setamp(struct qreg *reg, const uint64_t i, _Complex double z)
{
	struct qreg_QuEST *qr = reg->data;
	
	double x = creal(z), y = cimag(z);
	setAmps(qr->reg, i, &x, &y, 1);
}

void qreg_zero(struct qreg *reg)
{
	struct qreg_QuEST *qr = reg->data;

	initBlankState(qr->reg);
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t num_codes)
{
	struct qreg_QuEST *qr = reg->data;

	for (size_t k = 0; k < num_codes; k++) {
		struct paulis code;
		paulis_merge(&code, reg->qb_lo, reg->qb_hi,
			codes_lo[k], code_hi);

		for (size_t i = 0; i < qr->num_qubits; i++)
			qr->tar_op[i] = paulis_get(code, i);
			
		multiRotatePauli(qr->reg, qr->tar_qb, qr->tar_op,
			qr->num_qubits, angles[k] * -2.0);
		
	}
}
