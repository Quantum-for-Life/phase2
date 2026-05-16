/*
 * b-qreg -- end-to-end micro-benchmarks for the
 * MPI register and the Pauli-rotation kernel
 * (phase2/qreg.h).
 *
 * Two fixed qubit sizes cross the CPU cache
 * hierarchy:
 *   - nqb=14 (16 K amps, ~256 KB, in L2 on most
 *     cores),
 *   - nqb=18 (262 K amps, ~4 MB, RAM-bound).
 *
 * Scenarios:
 *   - qreg_init  (one per size)
 *   - paulirot   (each size x ncodes in {1, 10, 100})
 *
 * Total bench time on a quiet host: ~3-4 s.
 */

#include "bench.h"

#include <complex.h>
#include <stdlib.h>

#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#define WD_SEED  UINT64_C(0x18c9ee04abeee30c)
#define SEED     UINT64_C(0x2d1da81dc94cf64f)

#define NUM_RUNS 11

static struct xoshiro256ss RNG;

static enum pauli_op rand_pauli(void)
{
	return (enum pauli_op)(xoshiro256ss_next(&RNG) % 4);
}

static void record(FILE *out, const struct bench_prov *prov, int mpi_ranks,
	const char *name, const char *params_json,
	double *samples, int num_runs)
{
	const struct bench_stats st = bench_compute_stats(samples, num_runs);

	char path[256];
	bench_runs_path(prov, path, sizeof path);

	struct bench_baseline bl;
	const bool has_bl = bench_find_baseline(path, prov->hostname,
		name, params_json, &bl);

	char label[80];
	snprintf(label, sizeof label, "%s %s", name, params_json);
	bench_print_row(label, &st, &bl, has_bl);

	bench_append_jsonl(out, prov, BENCH_BACKEND, mpi_ranks,
		name, params_json, num_runs, &st);
}


/* -- qreg_init -- alloc + zero of the amp vector --------------------------- */

static void measure_qreg_init(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, uint32_t nqb)
{
	/* Warm-up: one init/free pair to fault the heap pages.  Without
	 * this the first timed call carries the first-touch cost and
	 * skews the max. */
	struct qreg warm;
	qreg_init(&warm, nqb);
	qreg_free(&warm);

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++) {
		const double t0 = bench_now_ns();
		struct qreg reg;
		qreg_init(&reg, nqb);
		const double t1 = bench_now_ns();
		samples[r] = t1 - t0;
		qreg_free(&reg);
	}

	char params[32];
	snprintf(params, sizeof params, "{\"nqb\":%u}", nqb);
	record(out, prov, mpi_ranks, "qreg_init", params, samples, NUM_RUNS);
}


/* -- paulirot -- the central rotation kernel ------------------------------- */

struct rot_setup {
	struct paulis code_hi;
	struct paulis *codes_lo;
	double *angles;
	size_t ncodes;
};

static int rot_setup_init(struct rot_setup *s, const struct qreg *reg,
	size_t ncodes)
{
	s->code_hi = paulis_new();
	s->codes_lo = malloc(sizeof *s->codes_lo * ncodes);
	s->angles = malloc(sizeof *s->angles * ncodes);
	if (!s->codes_lo || !s->angles) {
		free(s->codes_lo);
		free(s->angles);
		return -1;
	}
	s->ncodes = ncodes;

	for (uint32_t k = reg->qb_lo; k < reg->qb_lo + reg->qb_hi; k++)
		paulis_set(&s->code_hi, rand_pauli(), k);
	for (size_t i = 0; i < ncodes; i++) {
		s->codes_lo[i] = paulis_new();
		for (uint32_t k = 0; k < reg->qb_lo; k++)
			paulis_set(&s->codes_lo[i], rand_pauli(), k);
		s->angles[i] = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
	}
	return 0;
}

static void rot_setup_free(struct rot_setup *s)
{
	free(s->codes_lo);
	free(s->angles);
}

static void measure_paulirot(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, uint32_t nqb, size_t ncodes)
{
	struct qreg reg;
	if (qreg_init(&reg, nqb) < 0) {
		fprintf(stderr, "b-qreg: qreg_init(%u) failed\n", nqb);
		return;
	}
	qreg_zero(&reg);

	struct rot_setup s;
	if (rot_setup_init(&s, &reg, ncodes) < 0) {
		fprintf(stderr, "b-qreg: rot_setup_init failed\n");
		qreg_free(&reg);
		return;
	}

	/* Warm-up: one untimed call brings page caches and the cache batch
	 * into a representative state. */
	qreg_paulirot(&reg, s.code_hi, s.codes_lo, s.angles, s.ncodes);

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++) {
		const double t0 = bench_now_ns();
		qreg_paulirot(&reg, s.code_hi, s.codes_lo, s.angles,
			s.ncodes);
		const double t1 = bench_now_ns();
		samples[r] = t1 - t0;
	}

	char params[64];
	snprintf(params, sizeof params, "{\"nqb\":%u,\"ncodes\":%zu}",
		nqb, ncodes);
	record(out, prov, mpi_ranks, "paulirot", params, samples, NUM_RUNS);

	rot_setup_free(&s);
	qreg_free(&reg);
}


/* -- main ------------------------------------------------------------------ */

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);
	struct world_info wd;
	world_info(&wd);
	xoshiro256ss_init(&RNG, SEED);

	struct bench_prov prov;
	bench_prov_init(&prov);

	FILE *out = NULL;
	if (wd.rank == 0) {
		char path[256];
		if (bench_runs_path(&prov, path, sizeof path) < 0) {
			fprintf(stderr,
				"b-qreg: cannot create bench/runs/\n");
			world_free();
			return 1;
		}
		out = fopen(path, "a");
		if (!out) {
			fprintf(stderr,
				"b-qreg: cannot open %s\n", path);
			world_free();
			return 1;
		}
		bench_print_header();
	}

	/* qreg requires nqb >= log2(world_size); skip the small size on
	 * larger MPI runs. */
	uint32_t nqb_min = 1;
	int sz = wd.size;
	while (sz >>= 1)
		nqb_min++;

	/* Two qubit sizes crossing the cache hierarchy.  Per-size ncodes
	 * arrays differ: at nqb=18 the ncodes=100 case is ~600 ms / call
	 * and 11 samples blows the fast-bench budget; ncodes={1, 10} is
	 * enough to observe RAM-bound regressions. */
	const uint32_t sizes[]            = {  14,  18 };
	const size_t   ncodes_per_size[2][3] = {
		{ 1, 10, 100 },   /* nqb=14: in L2, all three cheap */
		{ 1, 10,   0 },   /* nqb=18: skip ncodes=100 */
	};

	for (size_t si = 0; si < sizeof sizes / sizeof sizes[0]; si++) {
		const uint32_t nqb = sizes[si];
		if (nqb < nqb_min)
			continue;
		measure_qreg_init(out, &prov, wd.size, nqb);
		for (size_t ci = 0; ci < 3; ci++) {
			const size_t nc = ncodes_per_size[si][ci];
			if (nc == 0)
				continue;
			measure_paulirot(out, &prov, wd.size, nqb, nc);
		}
	}

	if (out)
		fclose(out);
	world_free();
	return 0;
}
