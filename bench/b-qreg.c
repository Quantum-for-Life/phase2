#include <stdint.h>
#include <stdio.h>

#include "phase2/qreg.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#include "bench.h"

#define REPS_MAX (999UL)
#define NQUBITS_MAX (29)


#define WD_SEED UINT64_C(0x18c9ee04abeee30c)
static struct world WD;

#define SEED UINT64_C(0x2d1da81dc94cf64f)
static struct xoshiro256ss RNG;

int b_qreg_init(void *data)
{
	uint32_t *nqb = data;
	struct qreg reg;

	qreg_init(&reg, *nqb);
	qreg_destroy(&reg);

	return 0;
}

int main(int argc, char **argv)
{
	(void)argc;
	(void)argv;

	struct bench b;

	world_init((void *)0, (void *)0, WD_SEED);
	world_info(&WD);

	log_info("MPI World size: %d", WD.size);
	uint32_t nqb_min = 1;
	uint64_t siz = WD.size;
	while (siz >>= 1)
		nqb_min++;
	log_info("nqb_min= %u, nqb_max= %u", nqb_min, NQUBITS_MAX);

	xoshiro256ss_init(&RNG, SEED);

	log_info("BENCHES >>>");
	log_info("b_qreg_init:");
	log_info("n_qb,t_ms");
	for (uint32_t n = nqb_min; n <= NQUBITS_MAX; n++) {
		bench_mark(&b, REPS_MAX, b_qreg_init, &n);
		log_info("%3u,%.6f", n, bench_msrep(b));
	}

	log_info("<<< END BENCHES");
	world_destroy();

	return 0;
}
