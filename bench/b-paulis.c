#include <stdint.h>

#include "phase2/paulis.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "bench.h"

#define REPS_MAX (999UL)

#define WD_SEED UINT64_C(0x0094be8e2d4ba8eb)
static struct world WD;

#define SEED UINT64_C(0xa208312e4065b1df)
static struct xoshiro256ss RNG;

#define WIDTH_MAX (64)

static int rand_pauli(void)
{
	return (int)(xoshiro256ss_next(&RNG) % 4);
}

struct b_paulis_set {
	struct paulis *p;
	int op;
	uint32_t n;
};

static int b_paulis_set(void *data)
{
	struct b_paulis_set *d = data;

	paulis_set(d->p, d->op, d->n);

	return 0;
}

static void measure_b_paulis_set(void)
{
	struct bench b;
	struct b_paulis_set d;
	struct paulis p = paulis_new();

	d.p = &p;
	d.op = rand_pauli();
	d.n = xoshiro256ss_next(&RNG) % WIDTH_MAX;

	bench_mark(&b, REPS_MAX, b_paulis_set, &d);
	log_info("b_paulis_set [t_ms]: %.6f", bench_msrep(b));
}

struct b_paulis_get {
	struct paulis p;
	uint32_t n;
};

static int b_paulis_get(void *data)
{
	struct b_paulis_get *d = data;

	paulis_get(d->p, d->n);

	return 0;
}

static void measure_b_paulis_get(void)
{
	struct bench b;
	struct b_paulis_get d;
	struct paulis p = paulis_new();

	d.p = p;
	d.n = xoshiro256ss_next(&RNG) % WIDTH_MAX;

	bench_mark(&b, REPS_MAX, b_paulis_get, &d);
	log_info("b_paulis_get [t_ms]: %.6f", bench_msrep(b));
}

struct b_paulis_effect {
	struct paulis p;
	uint64_t i;
	_Complex double *z;
};

static int b_paulis_effect(void *data)
{
	struct b_paulis_effect *d = data;

	paulis_effect(d->p, d->i, d->z);

	return 0;
}

static void measure_b_paulis_effect(void)
{
	struct bench b;
	struct b_paulis_effect d;
	struct paulis p = paulis_new();
	_Complex double z = rand_dbl01(&RNG);

	for (uint32_t k = 0; k < WIDTH_MAX; k++)
		paulis_set(&p, rand_pauli(), k);
	d.p = p;
	d.i = xoshiro256ss_next(&RNG);
	d.z = &z;

	bench_mark(&b, REPS_MAX, b_paulis_effect, &d);
	log_info("b_paulis_effect [t_ms]: %.6f", bench_msrep(b));
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

	xoshiro256ss_init(&RNG, SEED);
	log_info("BENCHES >>>");

	measure_b_paulis_set();
	measure_b_paulis_get();
	measure_b_paulis_effect();

	log_info("<<< END BENCHES");
	log_info("Destroy world");
	world_destroy();

	return 0;
}
