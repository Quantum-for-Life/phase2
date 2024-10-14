#include <complex.h>
#include <stdint.h>
#include <stdio.h>

#include "phase2/qreg.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "bench.h"

#define REPS_MAX (999UL)

#define NQB_MAX (29)
static uint32_t NQB_MIN = 1;

#define WD_SEED UINT64_C(0x18c9ee04abeee30c)
static struct world WD;

#define SEED UINT64_C(0x2d1da81dc94cf64f)
static struct xoshiro256ss RNG;

static int b_qreg_init(void *data)
{
	uint32_t *nqb = data;
	struct qreg reg;

	qreg_init(&reg, *nqb);
	qreg_destroy(&reg);

	return 0;
}

static void b_qreg_init_measure(void)
{
	struct bench b;

	log_info("n_qb,t_ms");
	for (uint32_t n = NQB_MIN; n <= NQB_MAX; n++) {
		bench_mark(&b, REPS_MAX, b_qreg_init, &n);
		log_info("%3u,%.6f", n, bench_msrep(b));
	}
}

struct b_qreg_get {
	struct qreg reg;
	_Complex double amp;
	uint64_t i;
};

static int b_qreg_get(void *data)
{
	struct b_qreg_get *d = data;

	qreg_getamp(&d->reg, d->i, &d->amp);

	return 0;
}

static void b_qreg_get_measure(void)
{
	struct bench b;
	
	log_info("n_qb,t_ms");
	for (uint32_t n = NQB_MIN; n <= NQB_MAX; n++) {
		struct b_qreg_get qg;
		qg.amp = 0.03 + I * 0.72;
		qg.i = (1 << n) - 1;

		qreg_init(&qg.reg, n);
		qreg_setamp(&qg.reg, qg.i, qg.amp);

		bench_mark(&b, REPS_MAX, b_qreg_get, &qg);
		log_info("%3u,%.6f", n, bench_msrep(b));

		qreg_destroy(&qg.reg);
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

	log_info("BENCHES >>>");
	xoshiro256ss_init(&RNG, SEED);

	log_info("b_qreg_init:");
	b_qreg_init_measure();

	log_info("b_qreg_get:");
	b_qreg_get_measure();

	log_info("<<< END BENCHES");

	log_info("Destroy world");
	world_destroy();

	return 0;
}
