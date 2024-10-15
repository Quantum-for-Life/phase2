#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "bench.h"

#define REPS_MAX (9UL)

#define NQB_MAX (29)
static uint32_t NQB_MIN = 1;

#define WD_SEED UINT64_C(0x18c9ee04abeee30c)
static struct world WD;

#define SEED UINT64_C(0x2d1da81dc94cf64f)
static struct xoshiro256ss RNG;

static double rand_double(void)
{
	return (double)(xoshiro256ss_next(&RNG) >> 11) * 0x1.0p-53;
}

static pauli_op_t rand_pauli(void)
{
	return (pauli_op_t)(xoshiro256ss_next(&RNG) % 4);
}

struct b_qreg_init {
	uint32_t n;
};

static int b_qreg_init(void *data)
{
	struct b_qreg_init *q = data;
	struct qreg reg;

	qreg_init(&reg, q->n);
	qreg_destroy(&reg);

	return 0;
}

static void measure_b_qreg_init(void)
{
	for (uint32_t n = NQB_MIN; n <= NQB_MAX; n++) {
		struct bench b;
		struct b_qreg_init q = { .n = n };

		bench_mark(&b, REPS_MAX, b_qreg_init, &q);
		log_info("b_qreg_init [nqb,t_ms]: %2u,%11.6f",
				n, bench_msrep(b));
	}
}

struct b_qreg_get {
	struct qreg *reg;
	_Complex double amp;
	uint64_t i;
};

static int b_qreg_get(void *data)
{
	struct b_qreg_get *d = data;

	qreg_getamp(d->reg, d->i, &d->amp);

	return 0;
}

static void measure_b_qreg_get(void)
{

	for (uint32_t n = NQB_MIN; n <= NQB_MAX; n++) {
		struct bench b;
		struct qreg reg;
		struct b_qreg_get q;

		qreg_init(&reg, n);
		qreg_setamp(&reg, q.i, q.amp);

		q.reg = &reg;
		q.amp = 0.03 + I * 0.72;
		q.i = (1 << n) - 1;

		bench_mark(&b, REPS_MAX, b_qreg_get, &q);
		log_info("b_qreg_get [nqb,t_ms]: %2u,%11.6f",
				n, bench_msrep(b));

		qreg_destroy(&reg);
	}
}

struct b_qreg_zero {
	struct qreg *reg;
};

static int b_qreg_zero(void *data)
{
	struct b_qreg_zero *d = data;

	qreg_zero(d->reg);

	return 0;
}

static void measure_b_qreg_zero(void)
{
	struct bench b;

	for (uint64_t n = NQB_MIN; n <= NQB_MAX; n++) {
		struct qreg reg;
		struct b_qreg_zero q;

		qreg_init(&reg, n);
		q.reg = &reg;

		bench_mark(&b, REPS_MAX, b_qreg_zero, &q);
		log_info("b_qreg_zero [nqb,t_ms]: %2u,%11.6f",
				n, bench_msrep(b));

		qreg_destroy(&reg);
	}
}

struct b_qreg_paulirot {
	struct qreg *reg;
	struct paulis code_hi, *codes_lo;
	double *angles;
	size_t ncodes;
};

static void b_qreg_paulirot_init(struct b_qreg_paulirot *q,
	struct qreg *reg, size_t ncodes)
{
	q->reg = reg;
	q->code_hi = paulis_new();
	if (!(q->codes_lo = malloc(sizeof(struct paulis) * ncodes)))
		exit(-1);
	if (!(q->angles = malloc(sizeof(double) * ncodes)))
		exit(-1);
	for (size_t i = 0; i < ncodes; i++) {
		q->codes_lo[i] = paulis_new();
		q->angles[i] = 0.0;
	}
	q->ncodes = ncodes;
}

static void b_qreg_paulirot_destroy(struct b_qreg_paulirot *q)
{
	free(q->codes_lo);
	free(q->angles);
}

static void b_qreg_paulirot_rand(struct b_qreg_paulirot *q)
{
	for (size_t k = q->reg->qb_lo; k < q->reg->qb_lo + q->reg->qb_hi; k++) {
		paulis_set(&q->code_hi, rand_pauli(), k);
	}
	for (size_t i = 0; i < q->ncodes; i++) {
		for (size_t k = 0; k < q->reg->qb_lo; k++)
			paulis_set(q->codes_lo + i, rand_pauli(), k);
		q->angles[i] = rand_double() * 2.0 - 1.0;
	}
}

static int b_qreg_paulirot(void *data)
{
	struct b_qreg_paulirot *q = data;

	qreg_paulirot(q->reg, q->code_hi, q->codes_lo, q->angles, q->ncodes);

	return 0;
}

static void measure_b_qreg_paulirot(void)
{
	for (uint64_t n = NQB_MIN; n <= NQB_MAX; n++) {
		for (size_t ncodes = 1; ncodes <= 100; ncodes *= 10) {

			struct bench b;
			struct qreg reg;
			struct b_qreg_paulirot q;

			qreg_init(&reg, n);
			b_qreg_paulirot_init(&q, &reg, ncodes);
			b_qreg_paulirot_rand(&q);

			bench_mark(&b, REPS_MAX, b_qreg_paulirot, &q);
			log_info("b_qreg_paulirot [nqb,ncodes,t_ms]: "
					"%2u,%3u,%11.6f",
				n, ncodes, bench_msrep(b));

			b_qreg_paulirot_destroy(&q);
			qreg_destroy(&reg);
		}
	}
}

int main(int argc, char **argv)
{
	(void)argc;
	(void)argv;

	world_init((void *)0, (void *)0, WD_SEED);
	world_info(&WD);
	log_info("Initialize world");

	log_info("backend: %s", WORLD_BACKEND);
	log_info("MPI World size: %d", WD.size);


	uint64_t siz = WD.size;
	while (siz >>= 1) NQB_MIN++;
	log_info("nqb_min= %u, nqb_max= %u", NQB_MIN, NQB_MAX);

	xoshiro256ss_init(&RNG, SEED);
	log_info("BENCHES >>>");

	measure_b_qreg_init();
	measure_b_qreg_get();
	measure_b_qreg_zero();
	measure_b_qreg_paulirot();

	log_info("<<< END BENCHES");

	log_info("Destroy world");
	world_destroy();

	return 0;
}
