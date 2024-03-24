#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "error.h"
#include "types.h"
#include "qreg.h"


#define MAX_COUNT (1 << 30)


int ev_init(struct ev *ev)
{
	int initialized, num_ranks, rank;

	MPI_Initialized(&initialized);
	if (!initialized && MPI_Init(NULL, NULL) != MPI_SUCCESS)
		return -EMPI;

	if (MPI_Comm_size(MPI_COMM_WORLD, &num_ranks) != MPI_SUCCESS)
		return -EMPI;
	if (num_ranks == 0)
		return -ESIZE;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
		return -EMPI;

	ev->num_ranks = num_ranks;
	ev->rank      = rank;

	return OK;
}

int ev_destroy(struct ev *ev)
{
	(void)ev;

	int finalized;
	MPI_Finalized(&finalized);
	if (!finalized && MPI_Finalize() != MPI_SUCCESS)
		return -EMPI;

	ev->num_ranks = 0;

	return OK;
}


struct paulis paulis_new(void)
{
	struct paulis code;
	code.pak[0] = (u64)0;
	code.pak[1] = (u64)0;

	return code;
}

struct paulis paulis_from_ops(const enum pauli_op *paulis, const u32 num_paulis)
{
	struct paulis code = paulis_new();
	for (u32 i = 0; i < num_paulis; i++)
		paulis_set(&code, paulis[i], i);

	return code;
}

void paulis_print(const struct paulis code)
{
	for (int i = 0; i < PAULI_MAX_WIDTH; i++) {
		const enum pauli_op pauli = paulis_get(code, i);
		printf("%c", PAULI_LABEL[pauli]);
	}
	printf("\n");
}

void paulis_set(struct paulis *code, const enum pauli_op pauli, const u32 n)
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
}

enum pauli_op paulis_get(const struct paulis code, const u32 n)
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

void paulis_mask(struct paulis *code, const u64 mask)
{
	code->pak[0] &= mask;
	code->pak[1] &= mask;
}

void paulis_shr(struct paulis *code, const u32 n)
{
	code->pak[0] >>= n;
	code->pak[1] >>= n;
}

u64 paulis_effect(struct paulis code, u64 i, root4 *z)
{
	const u64 *p = code.pak;
	const u64  j = i ^ p[0];

	if (z == NULL)
		return j;

	int root4 = 0;
	for (u64 yp = ~i & p[0] & p[1]; yp > 0; yp >>= 1)
		if (yp & 1)
			root4 += 1;
	for (u64 zm = i & ~p[0] & p[1]; zm > 0; zm >>= 1)
		if (zm & 1)
			root4 += 2;
	for (u64 ym = i & p[0] & p[1]; ym > 0; ym >>= 1)
		if (ym & 1)
			root4 += 3;
	*z = root4 & 3;

	return j;
}

void paulis_split(const struct paulis code, const u32 qb_lo, const u32 qb_hi,
	struct paulis *lo, struct paulis *hi)
{
	const u64 mask_lo = ((u64)1 << qb_lo) - 1;

	*lo = code;
	paulis_mask(lo, mask_lo);

	*hi = code;
	paulis_mask(hi, ((u64)1 << (qb_lo + qb_hi)) - mask_lo);
}


// Assume n > 0
static u32 ulog2(const u64 n)
{
	u32 l = 0;
	for (u64 i = n >> 1; i > 0; i >>= 1)
		l++;

	return l;
}

int qreg_init(struct qreg *reg, const u32 num_qubits, const struct ev *ev)
{
	int ret;

	const u32 qb_hi = ulog2(ev->num_ranks);
	if (qb_hi >= num_qubits)
		return -ESIZE;
	const u32 qb_lo	   = num_qubits - qb_hi;
	const u64 num_amps = (u64)1 << qb_lo;

	const int    msg_count = num_amps < MAX_COUNT ? num_amps : MAX_COUNT;
	const size_t num_reqs  = num_amps / msg_count;

	MPI_Request *const reqs = malloc(sizeof(MPI_Request) * num_reqs * 4);
	if (reqs == NULL) {
		ret = -EMEM;
		goto err_reqs_alloc;
	}
	fl *const amp = malloc(sizeof(fl) * num_amps * 4);
	if (amp == NULL) {
		ret = -EMEM;
		goto err_amp_alloc;
	}

	reg->ev	       = ev;
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

	return OK;

	// free(amp);
err_amp_alloc:
	free(reqs);
err_reqs_alloc:
	return ret;
}

void qreg_destroy(struct qreg *reg)
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

struct distidx {
	int    rank;
	size_t i;
};

static struct distidx distidx_fromn(
	const u32 qb_lo, const u32 qb_hi, const u64 n)
{
	const u64 mask_lo = ((u64)1 << qb_lo) - 1;
	const u64 mask_hi = ((u64)1 << (qb_hi + qb_lo)) - mask_lo;

	struct distidx di;
	di.rank = (n & mask_hi) >> qb_lo;
	di.i	= n & mask_lo;

	return di;
}

void qreg_getamp(const struct qreg *reg, const u64 n, fl (*z)[2])
{
	const struct distidx di = distidx_fromn(reg->qb_lo, reg->qb_hi, n);

	if (reg->ev->rank == di.rank) {
		(*z)[0] = reg->amp[0][di.i];
		(*z)[1] = reg->amp[1][di.i];
	};
	MPI_Bcast(z, 2, MPI_DOUBLE, di.rank, MPI_COMM_WORLD);
}

void qreg_setamp(struct qreg *reg, const u64 n, const fl z[2])
{
	const struct distidx di = distidx_fromn(reg->qb_lo, reg->qb_hi, n);

	if (reg->ev->rank == di.rank) {
		reg->amp[0][di.i] = z[0];
		reg->amp[1][di.i] = z[1];
	};

	MPI_Barrier(MPI_COMM_WORLD);
}

void qreg_blank(struct qreg *reg)
{
	for (u64 i = 0; i < reg->num_amps; i++) {
		reg->amp[0][i] = 0.0;
		reg->amp[1][i] = 0.0;
	}
}

void qreg_zerostate(struct qreg *reg)
{
	qreg_blank(reg);
	qreg_setamp(reg, 0, (fl[]){ 1.0, 0.0 });
}

static void qreg_exchbuf_init(struct qreg *reg, const int rnk_rem)
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

static void qreg_exchbuf_waitall(struct qreg *reg)
{
	const int nr = reg->num_reqs;

	MPI_Waitall(2 * nr, reg->reqs_snd, MPI_STATUSES_IGNORE);
	MPI_Waitall(2 * nr, reg->reqs_rcv, MPI_STATUSES_IGNORE);
}

static void kernel_sep(fl *amp, fl *buf, const u64 i)
{
	const fl a = amp[i];
	const fl b = buf[i];

	amp[i] = (a + b) / 2.0;
	buf[i] = (a - b) / 2.0;
}

static void kernel_rot(fl *amp[2], const struct paulis code, const u64 i,
	const u64 j, const fl eip[2])
{
	root4 zi, zj;
	paulis_effect(code, i, &zj);
	paulis_effect(code, j, &zi);

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

static void kernel_add(fl *amp, const fl *buf, const u64 i)
{
	amp[i] += buf[i];
}

static void kernel_mul_cpx_scalar(fl *amp[2], const u64 i, const root4 z)
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

/* All codes are assumed to share the same hi code */
void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const fl *angles, const size_t num_codes)
{
	/* Compute permutation from outer qubits */
	const u64 rnk_loc = reg->ev->rank;
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
		kernel_mul_cpx_scalar(reg->buf, i, buf_mul);
		kernel_sep(reg->amp[0], reg->buf[0], i);
		kernel_sep(reg->amp[1], reg->buf[1], i);
	}
	for (size_t k = 0; k < num_codes; k++) {
		const fl eip_amp[2] = { cos(angles[k]), sin(angles[k]) };
		const fl eip_buf[2] = { eip_amp[0], -eip_amp[1] };

		for (u64 i = 0; i < reg->num_amps; i++) {
			const u64 j = paulis_effect(codes_lo[k], i, NULL);
			if (j < i)
				continue;

			kernel_rot(reg->amp, codes_lo[k], i, j, eip_amp);
			kernel_rot(reg->buf, codes_lo[k], i, j, eip_buf);
		}
	}
	for (u64 i = 0; i < reg->num_amps; i++) {
		kernel_add(reg->amp[0], reg->buf[0], i);
		kernel_add(reg->amp[1], reg->buf[1], i);
	}

	// MPI_Barrier(MPI_COMM_WORLD);
}
