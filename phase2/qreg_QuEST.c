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
	size_t num_qubits;
	Qureg reg;
	int tar_qb[QREG_MAX_WIDTH];
	int tar_op[QREG_MAX_WIDTH];
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
	if (!qr)
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

/* This is based on QuEST's implementation of statevec_multiRotatePauli() */
static void quest_paulirot(
	Qureg qureg, int *tar_qb, int *tar_op, size_t num_tar, double angle)
{
	size_t mask = 0;
	const qreal f = 1 / sqrt(2);
	Complex uRxAlpha = { .real = f, .imag = 0 };
	Complex uRxBeta = { .real = 0, .imag = -f };
	Complex uRyAlpha = { .real = f, .imag = 0 };
	Complex uRyBeta = { .real = -f, .imag = 0 };

	/* rotate basis so that exp(Z) will effect exp(Y) and exp(X) */
	for (int i = 0; i < num_tar; i++) {
		switch (tar_op[i]) {
		case PAULI_I:
			continue;
		case PAULI_X:
			compactUnitary(qureg, tar_qb[i], uRyAlpha, uRyBeta);
			break;
		case PAULI_Y:
			compactUnitary(qureg, tar_qb[i], uRxAlpha, uRxBeta);
			break;
		case PAULI_Z:
			break;
		default:
			__builtin_unreachable();
		}
		mask |= UINT64_C(1) << i;
	}

	/* Apply diagonal operator */
	_Complex double z, z_ph = cexp(I * angle);
	for (size_t i = 0; i < qureg.numAmpsPerChunk; i++) {
		z = qureg.stateVec.real[i] + I * qureg.stateVec.imag[i];

		size_t idx = qureg.numAmpsPerChunk * qureg.chunkId + i;
		if (__builtin_popcountll(mask & idx) % 2 == 0)
			z *= z_ph;
		else
			z *= conj(z_ph);

		qureg.stateVec.real[i] = creal(z);
		qureg.stateVec.imag[i] = cimag(z);
	}

	/* undo X and Y basis rotations */
	uRxBeta.imag *= -1;
	uRyBeta.real *= -1;
	for (int i = 0; i < num_tar; i++) {
		switch (tar_op[i]) {
		case PAULI_I:
			break;
		case PAULI_X:
			compactUnitary(qureg, tar_qb[i], uRyAlpha, uRyBeta);
			break;
		case PAULI_Y:
			compactUnitary(qureg, tar_qb[i], uRxAlpha, uRxBeta);
			break;
		case PAULI_Z:
			break;
		default:
			__builtin_unreachable();
		}
	}
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t num_codes)
{
	struct paulis code;
	struct qreg_QuEST *qr = reg->data;

	for (size_t k = 0; k < num_codes; k++) {
		paulis_merge(
			&code, reg->qb_lo, reg->qb_hi, codes_lo[k], code_hi);
		for (size_t i = 0; i < qr->num_qubits; i++)
			qr->tar_op[i] = paulis_get(code, i);

		quest_paulirot(qr->reg, qr->tar_qb, qr->tar_op, qr->num_qubits,
			angles[k]);
	}
}
