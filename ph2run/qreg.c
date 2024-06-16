#include <complex.h>
#include <math.h>
#include <stdlib.h>

#include "mpi.h"

#include "qreg.h"

#define MAX_COUNT (1 << 29)

int ev_init(struct qreg_ev *ev)
{
	int init, nrk, rk;

	MPI_Initialized(&init);
	if (!init && MPI_Init(NULL, NULL) != MPI_SUCCESS)
		return -1;

	if (MPI_Comm_size(MPI_COMM_WORLD, &nrk) != MPI_SUCCESS)
		return -1;
	if (nrk == 0)
		return -1;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rk) != MPI_SUCCESS)
		return -1;

	ev->num_ranks = nrk;
	ev->rank      = rk;

	return 0;
}

struct paulis paulis_new(void)
{
	struct paulis code = {
		.pak = {0, 0}
	};

	return code;
}

static int paulis_countis(struct paulis code)
{
	return __builtin_popcountll(code.pak[0] & code.pak[1]);
}

void paulis_set(
	struct paulis *code, const enum pauli_op pauli, const uint32_t n)
{
	const uint64_t n_mask = UINT64_C(1) << n;

	switch (pauli) {
	case PAULI_I:
		code->pak[0] &= ~n_mask;
		code->pak[1] &= ~n_mask;
		break;
	case PAULI_X:
		code->pak[0] |= n_mask;
		code->pak[1] &= ~n_mask;
		break;
	case PAULI_Y:
		code->pak[0] |= n_mask;
		code->pak[1] |= n_mask;
		break;
	case PAULI_Z:
		code->pak[0] &= ~n_mask;
		code->pak[1] |= n_mask;
		break;
	}
}

enum pauli_op paulis_get(const struct paulis code, const uint32_t n)
{
	int pa = 0;
	pa |= code.pak[0] >> n & 1;
	pa |= (code.pak[1] >> n & 1) << 1;

	switch (pa) {
	case 0:
		return PAULI_I;
	case 1:
		return PAULI_X;
	case 2:
		return PAULI_Z;
	case 3:
		return PAULI_Y;
	default:
		__builtin_unreachable();
	}
}

int paulis_eq(const struct paulis code1, const struct paulis code2)
{
	return code1.pak[0] == code2.pak[0] && code1.pak[1] == code2.pak[1];
}

void paulis_mask(struct paulis *code, const uint64_t mask)
{
	code->pak[0] &= mask;
	code->pak[1] &= mask;
}

void paulis_shr(struct paulis *code, const uint32_t n)
{
	code->pak[0] >>= n;
	code->pak[1] >>= n;
}

uint64_t paulis_effect(
	const struct paulis code, const uint64_t i, _Complex double *z)
{
	uint64_t j = i ^ code.pak[0];
	if (z != NULL) {
		const int minuses = __builtin_popcountll(j & code.pak[1]);
		int	  root4	  = (paulis_countis(code) + 2 * minuses) & 0x3;
		switch (root4) {
		case 0:
			break;
		case 1:
			*z *= _Complex_I;
			break;
		case 2:
			*z *= -1.0;
			break;
		case 3:
			*z *= -_Complex_I;
			break;
		default:
			__builtin_unreachable();
		}
	}

	return j;
}

void paulis_split(const struct paulis code, const uint32_t qb_lo,
	const uint32_t qb_hi, struct paulis *lo, struct paulis *hi)
{
	const uint64_t mask_lo = ((uint64_t)1 << qb_lo) - 1;
	const uint64_t mask_hi = ((uint64_t)1 << (qb_hi + qb_lo)) - 1;

	lo->pak[0] = code.pak[0];
	lo->pak[1] = code.pak[1];
	paulis_mask(lo, mask_lo);

	hi->pak[0] = code.pak[0];
	hi->pak[1] = code.pak[1];
	paulis_mask(hi, mask_hi - mask_lo);
}

int qreg_init(struct qreg *reg, const uint32_t num_qubits)
{
	if (ev_init(&reg->ev) < 0)
		goto err_ev;

	uint32_t qb_hi = 0, nrk = reg->ev.num_ranks;
	while (nrk >>= 1)
		qb_hi++;
	if (qb_hi >= num_qubits)
		return -1;
	const uint32_t qb_lo	= num_qubits - qb_hi;
	const uint64_t num_amps = UINT64_C(1) << qb_lo;

	const int    msg_count = num_amps < MAX_COUNT ? num_amps : MAX_COUNT;
	const size_t num_reqs  = num_amps / msg_count;

	MPI_Request *reqs = malloc(sizeof *reqs * num_reqs * 2);
	if (!reqs)
		goto err_reqs_alloc;
	_Complex double *amp = malloc(sizeof *amp * num_amps * 2);
	if (!amp)
		goto err_amp_alloc;

	reg->qb_lo     = qb_lo;
	reg->qb_hi     = qb_hi;
	reg->amp       = amp;
	reg->buf       = amp + num_amps;
	reg->num_amps  = num_amps;
	reg->msg_count = msg_count;
	reg->reqs_snd  = reqs;
	reg->reqs_rcv  = reqs + num_reqs;
	reg->num_reqs  = num_reqs;

	return 0;

	// free(amp);
err_amp_alloc:
	free(reqs);
err_reqs_alloc:
err_ev:
	return -1;
}

void qreg_destroy(struct qreg *reg)
{
	if (reg->amp)
		free(reg->amp);
	reg->amp = NULL;
	reg->buf = NULL;

	if (reg->reqs_snd)
		free(reg->reqs_snd);
	reg->reqs_snd = NULL;
	reg->reqs_rcv = NULL;
}

struct remote_idx {
	int	 rank;
	uint64_t i;
};

static struct remote_idx calc_remote_idx(
	const uint32_t qb_lo, const uint32_t qb_hi, const uint64_t i)
{
	const uint64_t mask_lo = (UINT64_C(1) << qb_lo) - 1;
	const uint64_t mask_hi = (UINT64_C(1) << (qb_hi + qb_lo)) - 1;

	struct remote_idx di;
	di.rank = (i & mask_hi) >> qb_lo;
	di.i	= i & mask_lo;

	return di;
}

void qreg_getamp(const struct qreg *reg, const uint64_t i, _Complex double *z)
{
	const struct remote_idx di = calc_remote_idx(reg->qb_lo, reg->qb_hi, i);

	if (reg->ev.rank == di.rank)
		*z = reg->amp[di.i];
	MPI_Bcast(z, 2, MPI_DOUBLE, di.rank, MPI_COMM_WORLD);
}

void qreg_setamp(struct qreg *reg, const uint64_t i, _Complex double z)
{
	const struct remote_idx di = calc_remote_idx(reg->qb_lo, reg->qb_hi, i);

	if (reg->ev.rank == di.rank)
		reg->amp[di.i] = z;
	MPI_Barrier(MPI_COMM_WORLD);
}

void qreg_zero(struct qreg *reg)
{
	_Complex double *z = reg->amp;
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

static void kernel_rot(_Complex double *amp, const uint64_t i,
	const _Complex double zi, const uint64_t j, const _Complex double zj,
	const _Complex double eip)
{
	_Complex double xi, xj;

	xi     = amp[i];
	xj     = amp[j];
	amp[i] = creal(eip) * xi + _Complex_I * cimag(eip) * zi * xj;
	amp[j] = creal(eip) * xj + _Complex_I * cimag(eip) * zj * xi;
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t num_codes)
{
	/* Compute permutation from outer qubits */
	const uint64_t rnk_loc = reg->ev.rank;
	const uint64_t rnk_rem = paulis_effect(code_hi, rnk_loc, NULL);
	qreg_exchbuf_init(reg, rnk_rem);

	/* Compute multiplication factor for the buffer
	   code_hi acts on the value of rank_remote now
	   (as if from receiving end). We discard the result. */
	_Complex double buf_mul = 1.0;
	paulis_effect(code_hi, rnk_loc, &buf_mul);

	qreg_exchbuf_waitall(reg);

	/* Compute permutation from inner qubits */
	for (uint64_t i = 0; i < reg->num_amps; i++) {
		reg->buf[i] *= buf_mul;

		_Complex a = reg->amp[i], b = reg->buf[i];
		reg->amp[i] = (a + b) / 2.0;
		reg->buf[i] = (a - b) / 2.0;
	}
	for (size_t k = 0; k < num_codes; k++) {
		const _Complex double eip = cexp(_Complex_I * angles[k]);
		for (uint64_t i = 0; i < reg->num_amps; i++) {
			_Complex double z = 1.0;
			const uint64_t	j = paulis_effect(codes_lo[k], i, &z);
			if (j < i)
				continue;

			kernel_rot(reg->amp, i, z, j, conj(z), eip);
			kernel_rot(reg->buf, i, z, j, conj(z), conj(eip));
		}
	}
	for (uint64_t i = 0; i < reg->num_amps; i++)
		reg->amp[i] += reg->buf[i];
}
