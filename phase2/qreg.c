#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"

/* --------------------------------------------------------------------------*/
/* Remove this section when C23 arrives.                                     */
#include <stdbool.h>
#define nullptr (void *)0
#define unreachable() (__builtin_unreachable())
/* --------------------------------------------------------------------------*/

#define MAX_COUNT (1 << 29)

typedef _Complex double c64;

/* Local copy of the world info. Initialized by qreg_init() */
static struct world WD;

int qreg_init(struct qreg *reg, const uint32_t nqb)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	uint32_t nqb_hi = 0, nrk = WD.size;
	/* nqb_hi = (int) log2(nrk) */
	while (nrk >>= 1)
		nqb_hi++;
	if (nqb_hi >= nqb)
		return -1;
	const uint32_t nqb_lo = nqb - nqb_hi;
	const uint64_t namp = UINT64_C(1) << nqb_lo;

	const int msg_count = namp < MAX_COUNT ? namp : MAX_COUNT;
	const size_t nreqs = namp / msg_count;

	MPI_Request *const reqs = malloc(sizeof(MPI_Request) * nreqs * 2);
	if (reqs == nullptr)
		goto err_reqs_alloc;
	c64 *const amp = malloc(sizeof(c64) * namp * 2);
	if (amp == nullptr)
		goto err_amp_alloc;

	reg->nqb_lo = nqb_lo;
	reg->nqb_hi = nqb_hi;
	reg->amp = amp;
	reg->buf = amp + namp;
	reg->namp = namp;
	reg->msg_count = msg_count;
	reg->reqs_snd = reqs;
	reg->reqs_rcv = reqs + nreqs;
	reg->nreqs = nreqs;

	return 0;

	// free(amp);
err_amp_alloc:
	free(reqs);
err_reqs_alloc:
	return -1;
}

void qreg_destroy(struct qreg *reg)
{
	if (reg->amp != nullptr)
		free(reg->amp);
	if (reg->reqs_snd != nullptr)
		free(reg->reqs_snd);
}

static void qb_split(uint64_t n, const uint32_t nqb_lo, const uint32_t nqb_hi,
	uint64_t *lo, uint64_t *hi)
{
	const uint64_t mask_lo = (UINT64_C(1) << nqb_lo) - 1;
	const uint64_t mask_hi = (UINT64_C(1) << nqb_hi) - 1;

	*lo = n & mask_lo;
	n >>= nqb_lo;
	*hi = n & mask_hi;
}

void qreg_getamp(const struct qreg *reg, const uint64_t i, c64 *z)
{
	uint64_t rank, loci;
	qb_split(i, reg->nqb_lo, reg->nqb_hi, &loci, &rank);

	if (WD.rank == (int)rank)
		*z = reg->amp[loci];
	MPI_Bcast(z, 2, MPI_DOUBLE, rank, MPI_COMM_WORLD);
}

void qreg_setamp(struct qreg *reg, const uint64_t i, c64 z)
{
	uint64_t rank, loci;
	qb_split(i, reg->nqb_lo, reg->nqb_hi, &loci, &rank);

	if (WD.rank == (int)rank)
		reg->amp[loci] = z;
	MPI_Barrier(MPI_COMM_WORLD);
}

void qreg_zero(struct qreg *reg)
{
	memset(reg->amp, 0, reg->namp * sizeof(c64));
}

static void exch_init(struct qreg *reg, const int rnk_rem)
{
	const int nr = reg->nreqs;

	for (int i = 0; i < nr; i++) {
		const size_t offset = i * reg->msg_count;

		MPI_Isend(reg->amp + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_snd + i);
		MPI_Irecv(reg->buf + offset, reg->msg_count * 2, MPI_DOUBLE,
			rnk_rem, i, MPI_COMM_WORLD, reg->reqs_rcv + i);
	}
}

static void exch_wait_amp(struct qreg *reg)
{
	MPI_Waitall(reg->nreqs, reg->reqs_snd, MPI_STATUSES_IGNORE);
}

static void exch_wait_buf(struct qreg *reg)
{
	MPI_Waitall(reg->nreqs, reg->reqs_rcv, MPI_STATUSES_IGNORE);
}

/* These kernels can be easily ported to CUDA. */
static __inline__ void kernel_mul(size_t i, c64 *a, c64 b)
{
	a[i] *= b;
}

static __inline__ void kernel_mix(size_t i, c64 *restrict a, c64 *restrict b)
{
	const c64 x = a[i];
	const c64 y = b[i];
	a[i] = (x + y) / 2.0;
	b[i] = (x - y) / 2.0;
}

static __inline__ void kernel_rot(
	size_t i, c64 *a, struct paulis code, double c, double s)
{
	c64 sz = s;
	const uint64_t j = paulis_effect(code, i, &sz);
	if (j < i)
		return;

	const c64 xi = a[i], xj = a[j];
	a[i] = c * xi + I * conj(sz) * xj;
	a[j] = c * xj + I * sz * xi;
}

static __inline__ void kernel_add(size_t i, c64 *restrict a, c64 *restrict b)
{
	a[i] += b[i];
}

static void qreg_paulirot_hi(struct qreg *reg, struct paulis code_hi)
{
	c64 bm = 1.0;

	paulis_shr(&code_hi, reg->nqb_lo);
	const uint64_t rnk_rem = paulis_effect(code_hi, WD.rank, &bm);

	exch_init(reg, rnk_rem);

	exch_wait_buf(reg);
	for (size_t i = 0; i < reg->namp; i++)
		kernel_mul(i, reg->buf, conj(bm));

	exch_wait_amp(reg);
}

static void qreg_paulirot_lo(struct qreg *reg, const struct paulis *codes_lo,
	const double *angles, const size_t ncodes)
{
	for (uint64_t i = 0; i < reg->namp; i++)
		kernel_mix(i, reg->amp, reg->buf);

	for (size_t k = 0; k < ncodes; k++) {
		const double c = cos(angles[k]), s = sin(angles[k]);
		for (size_t i = 0; i < reg->namp; i++)
			kernel_rot(i, reg->amp, codes_lo[k], c, s);
		for (size_t i = 0; i < reg->namp; i++)
			kernel_rot(i, reg->buf, codes_lo[k], c, -s);
	}

	for (uint64_t i = 0; i < reg->namp; i++)
		kernel_add(i, reg->amp, reg->buf);
}

void qreg_paulirot(struct qreg *reg, const struct paulis code_hi,
	const struct paulis *codes_lo, const double *angles,
	const size_t ncodes)
{
	/* Compute the action on hi qubits. */
	qreg_paulirot_hi(reg, code_hi);

	/* Compute the action on lo qubits. */
	qreg_paulirot_lo(reg, codes_lo, angles, ncodes);
}
