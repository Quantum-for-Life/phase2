/*
 * b-prob -- micro-benchmarks for the CDF kernel
 * (include/prob.h).
 *
 * Two operations:
 *
 *   cdf_from_array   build a CDF from n weighted
 *                    samples; O(n) two-pass.
 *
 *   cdf_inverse      inverse-CDF lookup for a
 *                    uniform draw y in [0, 1];
 *                    O(log n) binary search.
 *
 * Sized at n=100 and n=1000 to bracket typical
 * Hamiltonian term counts in qdrift / cmpsit.
 *
 * The inverse cell pre-draws a fixed array of
 * uniform y values at warm-up so the timed loop
 * measures only the binary search, not the
 * xoshiro RNG.
 *
 * Total bench time on a quiet host: ~8 s
 * (~7 s MPI init + ~500 ms measurements).
 */

#define LOG_SUBSYS "b-prob"

#include "bench.h"

#include <stdlib.h>

#include "phase2/world.h"
#include "prob.h"
#include "xoshiro256ss.h"

#define WD_SEED  UINT64_C(0x53e1f4b6c8a209d7)
#define SEED     UINT64_C(0x9a4c712f8e5036bd)

#define NUM_RUNS 11
#define N_SMALL  100
#define N_LARGE  1000

#define MOM_K_BUILD    10
#define INNER_BUILD    1000

#define MOM_K_INVERSE  100
#define INNER_INVERSE  10000

#define Y_DRAWS        4096


static struct xoshiro256ss RNG;

/* Sink so the compiler can't fold prob_cdf_inverse's
 * return value when it's the only use. */
static volatile size_t g_sink;


/* -- table row + JSONL append ----------------------------------------------*/

static void record(FILE *out, const struct bench_prov *prov, int mpi_ranks,
	const char *name, const char *params_json, const char *display,
	double *samples, int sub_samples)
{
	const struct bench_stats st = bench_compute_stats(samples, NUM_RUNS);

	char path[256];
	bench_runs_path(prov, path, sizeof path);

	struct bench_baseline bl;
	const bool has_bl = bench_find_baseline(path, prov->hostname,
		name, params_json, sub_samples, &bl);

	bench_print_row(display, sub_samples, &st, &bl, has_bl);

	bench_append_jsonl(out, prov, BENCH_BACKEND, mpi_ranks,
		name, params_json, NUM_RUNS, sub_samples, &st);
}


/* -- scenarios -------------------------------------------------------------*/

static void measure_cdf_from_array(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, size_t n)
{
	double *weights = malloc(sizeof *weights * n);
	if (!weights) {
		fprintf(stderr, "b-prob: oom (weights n=%zu)\n", n);
		return;
	}
	for (size_t i = 0; i < n; i++)
		weights[i] = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;

	struct prob_cdf cdf;
	if (prob_cdf_init(&cdf, n) < 0) {
		fprintf(stderr, "b-prob: prob_cdf_init failed (n=%zu)\n", n);
		free(weights);
		return;
	}

	/* Warm-up. */
	prob_cdf_from_array_strided(&cdf, weights, sizeof weights[0], NULL);

	double samples[NUM_RUNS];
	BENCH_SAMPLE_LOOP(samples, NUM_RUNS, MOM_K_BUILD, INNER_BUILD, ({
		prob_cdf_from_array_strided(&cdf, weights,
			sizeof weights[0], NULL);
	}));

	char params[40];
	snprintf(params, sizeof params, "{\"n\":%zu}", n);
	char display[40];
	snprintf(display, sizeof display, "cdf_from_array n=%zu", n);
	record(out, prov, mpi_ranks, "cdf_from_array", params, display,
		samples, MOM_K_BUILD);

	prob_cdf_free(&cdf);
	free(weights);
}

static void measure_cdf_inverse(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, size_t n)
{
	double *weights = malloc(sizeof *weights * n);
	if (!weights) {
		fprintf(stderr, "b-prob: oom (weights n=%zu)\n", n);
		return;
	}
	for (size_t i = 0; i < n; i++)
		weights[i] = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;

	struct prob_cdf cdf;
	if (prob_cdf_init(&cdf, n) < 0) {
		fprintf(stderr, "b-prob: prob_cdf_init failed (n=%zu)\n", n);
		free(weights);
		return;
	}
	prob_cdf_from_array_strided(&cdf, weights, sizeof weights[0], NULL);

	/* Pre-draw uniform y values so the timed loop measures only
	 * prob_cdf_inverse, not the RNG. */
	double ys[Y_DRAWS];
	for (int k = 0; k < Y_DRAWS; k++)
		ys[k] = xoshiro256ss_dbl01(&RNG);

	/* Warm-up. */
	g_sink = prob_cdf_inverse(&cdf, ys[0]);

	double samples[NUM_RUNS];
	BENCH_SAMPLE_LOOP(samples, NUM_RUNS, MOM_K_INVERSE, INNER_INVERSE, ({
		g_sink = prob_cdf_inverse(&cdf, ys[_bi & (Y_DRAWS - 1)]);
	}));

	char params[40];
	snprintf(params, sizeof params, "{\"n\":%zu}", n);
	char display[40];
	snprintf(display, sizeof display, "cdf_inverse n=%zu", n);
	record(out, prov, mpi_ranks, "cdf_inverse", params, display,
		samples, MOM_K_INVERSE);

	prob_cdf_free(&cdf);
	free(weights);
}


/* -- main ------------------------------------------------------------------*/

int main(void)
{
	setenv("PHASE2_LOG", "warn", 0);

	bench_pin_cpu(0);

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
			fprintf(stderr, "b-prob: cannot create bench/runs/\n");
			world_free();
			return 1;
		}
		out = fopen(path, "a");
		if (!out) {
			fprintf(stderr, "b-prob: cannot open %s\n", path);
			world_free();
			return 1;
		}
		bench_print_banner("b-prob");
		bench_print_header();
	}

	measure_cdf_from_array(out, &prov, wd.size, N_SMALL);
	measure_cdf_from_array(out, &prov, wd.size, N_LARGE);
	measure_cdf_inverse  (out, &prov, wd.size, N_SMALL);
	measure_cdf_inverse  (out, &prov, wd.size, N_LARGE);

	if (out)
		fclose(out);
	world_free();
	return 0;
}
