/*
 * state_prep_coeff: _expand + _inner timing for two
 * sizes covering N4-class and N8-class fixtures.
 * Programmatic fixture: random C_alpha / C_beta from
 * xoshiro256ss; no HDF5 dependency.  Total bench wall
 * ~1 s on a quiet host.
 */

#define LOG_SUBSYS "b-state-prep-coeff"

#include "bench.h"

#include <complex.h>
#include <stdlib.h>

#include "phase2.h"
#include "phase2/state_prep_coeff.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#define WD_SEED  UINT64_C(0xb47c0a82d51e9637)
#define SEED     UINT64_C(0x9e3d5a17c0b46f82)

#define NUM_RUNS 11

static struct xoshiro256ss RNG;

static volatile _Complex double g_sink_z;


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

/*
 * For each (n_sites, n_alpha, n_beta) cell: build a
 * random C_alpha (and reuse it as C_beta), allocate
 * qreg + scratch, time expand and inner.
 */
static void measure_size(FILE *out, const struct bench_prov *prov, int mpi_ranks,
	uint32_t n_sites, uint32_t n_alpha, uint32_t n_beta,
	int expand_K, int expand_inner, int inner_K, int inner_inner)
{
	const size_t Ca_len = (size_t)n_sites * (n_alpha ? n_alpha : 1);
	double *Ca = malloc(sizeof *Ca * Ca_len);
	if (!Ca) {
		fprintf(stderr, "b-state-prep-coeff: oom Ca\n");
		return;
	}
	for (size_t i = 0; i < Ca_len; i++)
		Ca[i] = xoshiro256ss_dbl01(&RNG) - 0.5;

	struct qreg reg;
	if (qreg_init(&reg, 2 * n_sites) < 0) {
		fprintf(stderr, "b-state-prep-coeff: qreg_init nqb=%u failed\n",
			2 * n_sites);
		free(Ca);
		return;
	}

	struct state_prep_coeff_scratch sc;
	if (state_prep_coeff_scratch_init(&sc, n_sites, n_alpha, n_beta) < 0) {
		fprintf(stderr, "b-state-prep-coeff: scratch_init failed\n");
		qreg_free(&reg);
		free(Ca);
		return;
	}

	/* Warm-up: prime the qreg pages + cache state. */
	state_prep_coeff_expand(&reg, &sc, Ca, NULL, 1.0, 0, 0);

	double samples[NUM_RUNS];

	for (int r = 0; r < NUM_RUNS; r++) {
		BENCH_SAMPLE_MIN(samples[r], expand_K, ({
			for (int b = 0; b < expand_inner; b++)
				state_prep_coeff_expand(&reg, &sc,
					Ca, NULL, 1.0, 0, 0);
		}));
		samples[r] /= (double)expand_inner;
	}

	char params[80];
	snprintf(params, sizeof params,
		"{\"n_sites\":%u,\"n_alpha\":%u,\"n_beta\":%u}",
		n_sites, n_alpha, n_beta);
	char display[40];
	snprintf(display, sizeof display, "expand ns=%u na=%u nb=%u",
		n_sites, n_alpha, n_beta);
	record(out, prov, mpi_ranks, "expand", params, display,
		samples, expand_K);

	/* state is already populated by warm-up + measurement
	 * loop above; ready for inner. */
	for (int r = 0; r < NUM_RUNS; r++) {
		BENCH_SAMPLE_MIN(samples[r], inner_K, ({
			for (int b = 0; b < inner_inner; b++) {
				_Complex double z = 0.0;
				state_prep_coeff_inner(&reg, &sc, Ca, NULL,
					1.0, 0, &z);
				g_sink_z = z;
			}
		}));
		samples[r] /= (double)inner_inner;
	}

	snprintf(display, sizeof display, "inner ns=%u na=%u nb=%u",
		n_sites, n_alpha, n_beta);
	record(out, prov, mpi_ranks, "inner", params, display,
		samples, inner_K);

	state_prep_coeff_scratch_free(&sc);
	qreg_free(&reg);
	free(Ca);
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
			fprintf(stderr,
				"b-state-prep-coeff: cannot create bench/runs/\n");
			world_free();
			return 1;
		}
		out = fopen(path, "a");
		if (!out) {
			fprintf(stderr,
				"b-state-prep-coeff: cannot open %s\n", path);
			world_free();
			return 1;
		}
		bench_print_banner("b-state-prep-coeff");
		bench_print_header();
	}

	/* N4 class (n_sites=4, na=nb=2): 36-term outer
	 * product, expand ~1 us per call.  K=100 INNER=10.
	 * N8 class (n_sites=8, na=nb=4): 4900-term outer,
	 * expand ~hundreds of us per call.  K=10 INNER=1. */
	measure_size(out, &prov, wd.size, 4, 2, 2,
		100, 10, 100, 10);
	measure_size(out, &prov, wd.size, 8, 4, 4,
		10, 1, 10, 1);

	if (out)
		fclose(out);
	world_free();
	return 0;
}
