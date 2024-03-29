#include <math.h>
#include <stdlib.h>

#include "mpi.h"

#include "qreg.h"

#define MAX_COUNT (1 << 30)

int
ev_init(struct qreg_ev *ev)
{
	int initialized, num_ranks, rank;

	MPI_Initialized(&initialized);
	if (!initialized && MPI_Init(NULL, NULL) != MPI_SUCCESS)
		return -QREG_EMPI;

	if (MPI_Comm_size(MPI_COMM_WORLD, &num_ranks) != MPI_SUCCESS)
		return -QREG_EMPI;
	if (num_ranks == 0)
		return -QREG_ESIZE;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
		return -QREG_EMPI;

	ev->num_ranks = num_ranks;
	ev->rank      = rank;

	return QREG_OK;
}

struct paulis
paulis_new(void)
{
	struct paulis code;
	code.pak[0] = (u64)0;
	code.pak[1] = (u64)0;
	code.is	    = R0;

	return code;
}

static void
paulis_update_iscount(struct paulis *code)
{
	code->is = __builtin_popcountll(code->pak[0] & code->pak[1]) & 0x3;
}

void
paulis_set(struct paulis *code, const enum pauli_op pauli, const u32 n)
{
	const u64 n_mask = (u64)1 << n;

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
	paulis_update_iscount(code);
}

enum pauli_op
paulis_get(const struct paulis code, const u32 n)
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

int
paulis_eq(const struct paulis code1, const struct paulis code2)
{
	return code1.pak[0] == code2.pak[0] && code1.pak[1] == code2.pak[1];
}

void
paulis_mask(struct paulis *code, const u64 mask)
{
	code->pak[0] &= mask;
	code->pak[1] &= mask;
	paulis_update_iscount(code);
}

void
paulis_shr(struct paulis *code, const u32 n)
{
	code->pak[0] >>= n;
	code->pak[1] >>= n;
	paulis_update_iscount(code);
}

u64
paulis_effect(const struct paulis code, const u64 i, root4 *z)
{
	u64 j = i ^ code.pak[0];
	if (z != NULL) {
		const int minuses = __builtin_popcountll(j & code.pak[1]);

		*z = (code.is + 2 * minuses) & 0x3;
	}

	return j;
}

void
paulis_split(const struct paulis code, const u32 qb_lo, const u32 qb_hi,
	struct paulis *lo, struct paulis *hi)
{
	const u64 mask_lo = ((u64)1 << qb_lo) - 1;

	*lo = code;
	paulis_mask(lo, mask_lo);
	paulis_update_iscount(lo);

	*hi = code;
	paulis_mask(hi, ((u64)1 << (qb_lo + qb_hi)) - mask_lo);
	paulis_update_iscount(hi);
}

// Assume n > 0
static u32
ulog2(const u64 n)
{
	u32 l = 0;
	for (u64 i = n >> 1; i > 0; i >>= 1)
		l++;

	return l;
}

int
qreg_init(struct qreg *reg, const u32 num_qubits)
{
	int ret;

	if (ev_init(&reg->ev) < 0) {
		ret = -QREG_EMPI;
		goto err_ev;
	}

	const u32 qb_hi = ulog2(reg->ev.num_ranks);
	if (qb_hi >= num_qubits)
		return -QREG_ESIZE;
	const u32 qb_lo	   = num_qubits - qb_hi;
	const u64 num_amps = (u64)1 << qb_lo;

	const int    msg_count = num_amps < MAX_COUNT ? num_amps : MAX_COUNT;
	const size_t num_reqs  = num_amps / msg_count;

	MPI_Request *const reqs = malloc(sizeof(MPI_Request) * num_reqs * 4);
	if (reqs == NULL) {
		ret = -QREG_ENOMEM;
		goto err_reqs_alloc;
	}
	fl *const amp = malloc(sizeof(fl) * num_amps * 4);
	if (amp == NULL) {
		ret = -QREG_ENOMEM;
		goto err_amp_alloc;
	}

	reg->qb_lo     = qb_lo;
	reg->qb_hi     = qb_hi;
	reg->amp[0]    = amp;
	reg->amp[1]    = amp + num_amps;
	reg->buf[0]    = amp + 2 * num_amps;
	reg->buf[1]    = amp + 3 * num_amps;
	reg->num_amps  = num_amps;
	reg->msg_count = msg_count;
	reg->reqs_snd  = reqs;
	reg->reqs_rcv  = reqs + 2 * num_reqs;
	reg->num_reqs  = num_reqs;

	return QREG_OK;

	// free(amp);
err_amp_alloc:
	free(reqs);
err_reqs_alloc:
err_ev:
	return ret;
}

void
qreg_destroy(struct qreg *reg)
{
	if (reg->amp[0])
		free(reg->amp[0]);
	reg->amp[0] = NULL;
	reg->amp[1] = NULL;
	reg->buf[0] = NULL;
	reg->buf[1] = NULL;

	if (reg->reqs_snd)
		free(reg->reqs_snd);
	reg->reqs_snd = NULL;
	reg->reqs_rcv = NULL;
}

struct remote_idx {
	int rank;
	u64 i;
};

static struct remote_idx
calc_remote_idx(const u32 qb_lo, const u32 qb_hi, const u64 i)
{
	const u64 mask_lo = ((u64)1 << qb_lo) - 1;
	const u64 mask_hi = ((u64)1 << (qb_hi + qb_lo)) - mask_lo;

	struct remote_idx di;
	di.rank = (i & mask_hi) >> qb_lo;
	di.i	= i & mask_lo;

	return di;
}

void
qreg_getamp(const struct qreg *reg, const u64 i, fl (*z)[2])
{
	const struct remote_idx di = calc_remote_idx(reg->qb_lo, reg->qb_hi, i);

	if (reg->ev.rank == di.rank) {
		(*z)[0] = reg->amp[0][di.i];
		(*z)[1] = reg->amp[1][di.i];
	};
	MPI_Bcast(z, 2, MPI_DOUBLE, di.rank, MPI_COMM_WORLD);
}

void
qreg_setamp(struct qreg *reg, const u64 i, const fl z[2])
{
	const struct remote_idx di = calc_remote_idx(reg->qb_lo, reg->qb_hi, i);

	if (reg->ev.rank == di.rank) {
		reg->amp[0][di.i] = z[0];
		reg->amp[1][di.i] = z[1];
	};
	MPI_Barrier(MPI_COMM_WORLD);
}

void
qreg_zero(struct qreg *reg)
{
	for (u64 i = 0; i < reg->num_amps; i++) {
		reg->amp[0][i] = 0.0;
		reg->amp[1][i] = 0.0;
	}
}

static void
qreg_exchbuf_init(struct qreg *reg, const int rnk_rem)
{
	const int nr = reg->num_reqs;

	for (int i = 0; i < nr; i++) {
		const size_t offset = i * reg->msg_count;

		MPI_Isend(reg->amp[0] + offset, reg->msg_count, MPI_DOUBLE,
			rnk_rem, 2 * i, MPI_COMM_WORLD, reg->reqs_snd + 2 * i);
		MPI_Isend(reg->amp[1] + offset, reg->msg_count, MPI_DOUBLE,
			rnk_rem, 2 * i + 1, MPI_COMM_WORLD,
			reg->reqs_snd + 2 * i + 1);

		MPI_Irecv(reg->buf[0] + offset, reg->msg_count, MPI_DOUBLE,
			rnk_rem, 2 * i, MPI_COMM_WORLD, reg->reqs_rcv + 2 * i);
		MPI_Irecv(reg->buf[1] + offset, reg->msg_count, MPI_DOUBLE,
			rnk_rem, 2 * i + 1, MPI_COMM_WORLD,
			reg->reqs_rcv + 2 * i + 1);
	}
}

static void
qreg_exchbuf_waitall(struct qreg *reg)
{
	const int nr = reg->num_reqs;

	MPI_Waitall(2 * nr, reg->reqs_snd, MPI_STATUSES_IGNORE);
	MPI_Waitall(2 * nr, reg->reqs_rcv, MPI_STATUSES_IGNORE);
}

static void
kernel_mul(fl *amp[2], const u64 i, const root4 z)
{
	const fl x = amp[0][i];
	const fl y = amp[1][i];

	switch (z) {
	case R0:
		break;
	case R1:
		amp[0][i] = -y;
		amp[1][i] = x;
		break;
	case R2:
		amp[0][i] = -x;
		amp[1][i] = -y;
		break;
	case R3:
		amp[0][i] = y;
		amp[1][i] = -x;
		break;
	}
}

static void
kernel_sep(fl *amp, fl *buf, const u64 i)
{
	const fl a = amp[i];
	const fl b = buf[i];

	amp[i] = (a + b) / 2.0;
	buf[i] = (a - b) / 2.0;
}

static void
kernel_rot(fl *amp[2], const u64 i, const root4 zi, const u64 j, const root4 zj,
	const fl eip[2])
{
	const fl xi_re = amp[0][i], xi_im = amp[1][i];
	const fl xj_re = amp[0][j], xj_im = amp[1][j];

	switch (zi) {
	case R0:
		amp[0][i] = eip[0] * xi_re - eip[1] * xj_im;
		amp[1][i] = eip[0] * xi_im + eip[1] * xj_re;
		break;
	case R1:
		amp[0][i] = eip[0] * xi_re - eip[1] * xj_re;
		amp[1][i] = eip[0] * xi_im - eip[1] * xj_im;
		break;
	case R2:
		amp[0][i] = eip[0] * xi_re + eip[1] * xj_im;
		amp[1][i] = eip[0] * xi_im - eip[1] * xj_re;
		break;
	case R3:
		amp[0][i] = eip[0] * xi_re + eip[1] * xj_re;
		amp[1][i] = eip[0] * xi_im + eip[1] * xj_im;
		break;
	}

	switch (zj) {
	case R0:
		amp[0][j] = eip[0] * xj_re - eip[1] * xi_im;
		amp[1][j] = eip[0] * xj_im + eip[1] * xi_re;
		break;
	case R1:
		amp[0][j] = eip[0] * xj_re - eip[1] * xi_re;
		amp[1][j] = eip[0] * xj_im - eip[1] * xi_im;
		break;
	case R2:
		amp[0][j] = eip[0] * xj_re + eip[1] * xi_im;
		amp[1][j] = eip[0] * xj_im - eip[1] * xi_re;
		break;
	case R3:
		amp[0][j] = eip[0] * xj_re + eip[1] * xi_re;
		amp[1][j] = eip[0] * xj_im + eip[1] * xi_im;
		break;
	}
}

/* All codes are assumed to share the same hi code */
void
qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, const size_t num_codes)
{
	/* Compute permutation from outer qubits */
	const u64 rnk_loc = reg->ev.rank;
	const u64 rnk_rem = paulis_effect(code_hi, rnk_loc, NULL);
	qreg_exchbuf_init(reg, rnk_rem);

	/* Compute multiplication factor for the buffer
	   code_hi acts on the value of rank_remote now
	   (as if from receiving end). We discard the result. */
	root4 buf_mul;
	paulis_effect(code_hi, rnk_rem, &buf_mul);

	qreg_exchbuf_waitall(reg);

	/* Compute permutation from inner qubits */
	for (u64 i = 0; i < reg->num_amps; i++) {
		kernel_mul(reg->buf, i, buf_mul);
		kernel_sep(reg->amp[0], reg->buf[0], i);
		kernel_sep(reg->amp[1], reg->buf[1], i);
	}
	for (size_t k = 0; k < num_codes; k++) {
		const fl eip_amp[2] = { cos(angles[k]), sin(angles[k]) };
		const fl eip_buf[2] = { eip_amp[0], -eip_amp[1] };

		for (u64 i = 0; i < reg->num_amps; i++) {
			root4 zi, zj;

			const u64 j = paulis_effect(codes_lo[k], i, &zi);
			if (j < i)
				continue;
			paulis_effect(codes_lo[k], j, &zj);

			kernel_rot(reg->amp, i, zi, j, zj, eip_amp);
			kernel_rot(reg->buf, i, zi, j, zj, eip_buf);
		}
	}
	for (u64 i = 0; i < reg->num_amps; i++) {
		reg->amp[0][i] += reg->buf[0][i];
		reg->amp[1][i] += reg->buf[1][i];
	}
}
