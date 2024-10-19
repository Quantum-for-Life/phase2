#include "c23_compat.h"
#include <complex.h>
#include <stdlib.h>

#include "QuEST.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

#include "qreg_impl.h"
#include "world_quest.h"

struct qreg_quest {
	Qureg qureg;
	int tg_qb[QREG_MAX_WIDTH];
	int tg_op[QREG_MAX_WIDTH];
};

int qreg_backend_init(struct qreg *reg)
{
	struct qreg_quest *const q = malloc(sizeof *q);
	if (!q)
		return -1;

	/* QuEST will allocate memory for the entire state.
	   We need to free the memory allocated by qreg during qreg_init
	   to make room for QuEST's Qureg. */
	free(reg->amp);
	reg->amp = nullptr;

	struct world_quest *const w = reg->wd.data;
	uint32_t nqb = reg->nqb_lo + reg->nqb_hi;
	q->qureg = createQureg(nqb, w->env);
	for (size_t i = 0; i < nqb; i++) {
		q->tg_qb[i] = i;
		q->tg_op[i] = PAULI_I;
	}

	reg->data = q;

	return 0;
}

void qreg_backend_destroy(struct qreg *reg)
{
	struct world_quest *w = reg->wd.data;
	struct qreg_quest *q = reg->data;

	destroyQureg(q->qureg, w->env);
}

void qreg_getamp(const struct qreg *reg, const uint64_t i, _Complex double *z)
{
	struct qreg_quest *q = reg->data;

	const Complex cz = getAmp(q->qureg, i);
	*z = cz.real + I * cz.imag;
}

void qreg_setamp(struct qreg *reg, const uint64_t i, _Complex double z)
{
	struct qreg_quest *q = reg->data;

	double x = creal(z), y = cimag(z);
	setAmps(q->qureg, i, &x, &y, 1);
}

void qreg_zero(struct qreg *reg)
{
	struct qreg_quest *q = reg->data;

	initBlankState(q->qureg);
}

/* This is based on QuEST's implementation of statevec_multiRotatePauli() */
static void paulirot_onecode(
	Qureg qureg, int *tg_qb, int *tg_op, size_t ntg, double angle)
{
	size_t mask = 0;
	const qreal f = 1.0 / sqrt(2);
	Complex ux_alpha = { .real = f, .imag = 0 };
	Complex ux_beta = { .real = 0, .imag = -f };
	Complex uy_alpha = { .real = f, .imag = 0 };
	Complex uy_beta = { .real = -f, .imag = 0 };

	/* rotate basis so that exp(Z) will effect exp(Y) and exp(X) */
	for (size_t i = 0; i < ntg; i++) {
		switch (tg_op[i]) {
		case PAULI_I:
			continue;
		case PAULI_X:
			compactUnitary(qureg, tg_qb[i], uy_alpha, uy_beta);
			break;
		case PAULI_Y:
			compactUnitary(qureg, tg_qb[i], ux_alpha, ux_beta);
			break;
		case PAULI_Z:
			break;
		default:
			unreachable();
		}
		mask |= UINT64_C(1) << i;
	}

	/* Apply diagonal operator */
	_Complex double z, z_ph = cexp(I * angle);
	for (int64_t i = 0; i < qureg.numAmpsPerChunk; i++) {
		z = qureg.stateVec.real[i] + I * qureg.stateVec.imag[i];

		size_t idx = qureg.numAmpsPerChunk * qureg.chunkId + i;
		if (stdc_count_ones_ul(mask & idx) % 2 == 0)
			z *= z_ph;
		else
			z *= conj(z_ph);

		qureg.stateVec.real[i] = creal(z);
		qureg.stateVec.imag[i] = cimag(z);
	}

	/* undo X and Y basis rotations */
	ux_beta.imag *= -1.0;
	uy_beta.real *= -1.0;
	for (size_t i = 0; i < ntg; i++) {
		switch (tg_op[i]) {
		case PAULI_I:
			break;
		case PAULI_X:
			compactUnitary(qureg, tg_qb[i], uy_alpha, uy_beta);
			break;
		case PAULI_Y:
			compactUnitary(qureg, tg_qb[i], ux_alpha, ux_beta);
			break;
		case PAULI_Z:
			break;
		default:
			unreachable();
		}
	}
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t ncodes)
{
	const size_t nqb = reg->nqb_lo + reg->nqb_hi;
	struct qreg_quest *q = reg->data;

	struct paulis code;
	for (size_t k = 0; k < ncodes; k++) {
		paulis_merge(
			&code, reg->nqb_lo, reg->nqb_hi, codes_lo[k], code_hi);
		for (size_t i = 0; i < nqb; i++)
			q->tg_op[i] = paulis_get(code, i);

		paulirot_onecode(q->qureg, q->tg_qb, q->tg_op, nqb, angles[k]);
	}
}
