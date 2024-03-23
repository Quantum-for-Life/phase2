#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "error.h"
#include "ev.h"
#include "paulis.h"
#include "qreg.h"

#define MAX_COUNT (1 << 30)

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
		ret = -ENOMEM;
		goto err_reqs_alloc;
	}
	fl *const amp = malloc(sizeof(fl) * num_amps * 4);
	if (amp == NULL) {
		ret = -ENOMEM;
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

void qreg_blankstate(struct qreg *reg)
{
	for (u64 i = 0; i < reg->num_amps; i++) {
		reg->amp[0][i] = 0.0;
		reg->amp[1][i] = 0.0;
	}
}

void qreg_zerostate(struct qreg *reg)
{
	qreg_blankstate(reg);
	const fl one[2] = { 1.0, 0.0 };
	qreg_setamp(reg, 0, one);
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

static void kernel_rot(
	fl *amp[2], struct paulis code, const u64 i, const fl eip[2])
{
	// TODO: simplify this. We know that z is a 4th root of unity
	fl  zj[2] = { 1, 0 };
	fl  zi[2] = { 1, 0 };
	u64 j	  = paulis_effect(code, i, &zj);
	paulis_effect(code, j, &zi);

	const fl xi_re = amp[0][i];
	const fl xi_im = amp[1][i];
	const fl xj_re = amp[0][j];
	const fl xj_im = amp[1][j];

	fl c1	  = eip[1] * zi[0];
	fl c2	  = eip[1] * zi[1];
	amp[0][i] = eip[0] * xi_re - c1 * xj_im - c2 * xj_re;
	amp[1][i] = eip[0] * xi_im + c1 * xj_re - c2 * xj_im;

	c1	  = eip[1] * zj[0];
	c2	  = eip[1] * zj[1];
	amp[0][j] = eip[0] * xj_re - c1 * xi_im - c2 * xi_re;
	amp[1][j] = eip[0] * xj_im + c1 * xi_re - c2 * xi_im;
}

static void kernel_add(fl *amp, const fl *buf, const u64 i)
{
	amp[i] += buf[i];
}

static void kernel_mul_cpx_scalar(fl *amp[2], const u64 i, fl z[2])
{
	const fl x = amp[0][i];
	const fl y = amp[1][i];

	amp[0][i] = x * z[0] - y * z[1];
	amp[1][i] = x * z[1] + y * z[0];
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
	fl buf_mul[2] = { 1.0, 0.0 };
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
			if (paulis_effect(codes_lo[k], i, NULL) < i)
				continue;

			kernel_rot(reg->amp, codes_lo[k], i, eip_amp);
			kernel_rot(reg->buf, codes_lo[k], i, eip_buf);
		}
	}
	for (u64 i = 0; i < reg->num_amps; i++) {
		kernel_add(reg->amp[0], reg->buf[0], i);
		kernel_add(reg->amp[1], reg->buf[1], i);
	}

	// MPI_Barrier(MPI_COMM_WORLD);
}
