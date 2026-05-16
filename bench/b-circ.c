/*
 * b-circ -- micro-benchmarks for the per-step
 * Trotter kernel (circ_step) and overlap
 * measurement (circ_measure) from
 * phase2/circ.h.
 *
 * Fixtures are built programmatically (no HDF5):
 * a Hamiltonian of `nterms` random Pauli terms
 * and a multidet state of `ndets` random basis
 * states with random amplitudes.  Sizes stay
 * cache-resident (nqb = 12, 4 K amplitudes).
 *
 * Scenarios:
 *
 *   circ_step nq=12 t=10  d=10  -- ~ms / call
 *   circ_step nq=12 t=100 d=10  -- ~10 ms / call (scaling)
 *   circ_measure nq=12 d=10     -- ~us / call
 *
 * circ_init takes ownership of the Hamiltonian
 * and multidet buffers; circ_free releases them
 * along with the cache and register.  Each
 * scenario builds + tears down its own circ
 * context.
 *
 * Total bench time on a quiet host: ~10 s.
 */

#define LOG_SUBSYS "b-circ"

#include "bench.h"

#include <complex.h>
#include <stdlib.h>

#include "phase2/circ.h"
#include "phase2/paulis.h"
#include "phase2/qreg.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#define WD_SEED  UINT64_C(0x6e2da8cb71f04539)
#define SEED     UINT64_C(0x8c75d3b0e4a2917b)

#define NUM_RUNS 11

static struct xoshiro256ss RNG;

/* Volatile sink so the compiler can't drop circ_measure's
 * return value when its only use is the sample timing. */
static volatile _Complex double g_sink_z;


/* -- fixture builders ------------------------------------------------------*/

static int build_hamil(struct circ_hamil *hm, uint32_t nqb, size_t nterms)
{
	if (circ_hamil_init(hm, nqb, nterms) < 0)
		return -1;
	for (size_t i = 0; i < nterms; i++) {
		hm->terms[i].op = paulis_new();
		for (uint32_t k = 0; k < nqb; k++)
			paulis_set(&hm->terms[i].op,
				(enum pauli_op)(xoshiro256ss_next(&RNG) % 4),
				k);
		hm->terms[i].cf = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
	}
	circ_hamil_sort_lex(hm);
	return 0;
}

static int build_muldet(struct circ_muldet *md, uint32_t nqb, size_t ndets)
{
	if (circ_muldet_init(md, ndets) < 0)
		return -1;
	const uint64_t mask = (nqb >= 64)
		? ~UINT64_C(0)
		: ((UINT64_C(1) << nqb) - 1);
	for (size_t i = 0; i < ndets; i++) {
		md->dets[i].idx = xoshiro256ss_next(&RNG) & mask;
		const double re = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
		const double im = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
		md->dets[i].cf = re + im * I;
	}
	return 0;
}


/* -- table row + JSONL append ----------------------------------------------*/

static void record(FILE *out, const struct bench_prov *prov, int mpi_ranks,
	const char *name, const char *params_json, const char *display,
	double *samples, int num_runs, int sub_samples)
{
	const struct bench_stats st = bench_compute_stats(samples, num_runs);

	char path[256];
	bench_runs_path(prov, path, sizeof path);

	struct bench_baseline bl;
	const bool has_bl = bench_find_baseline(path, prov->hostname,
		name, params_json, sub_samples, &bl);

	bench_print_row(display, sub_samples, &st, &bl, has_bl);

	bench_append_jsonl(out, prov, BENCH_BACKEND, mpi_ranks,
		name, params_json, num_runs, sub_samples, &st);
}


/* -- scenarios -------------------------------------------------------------*/

static void measure_circ_step(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, uint32_t nqb, size_t nterms, size_t ndets, int K)
{
	struct circ_hamil hm;
	if (build_hamil(&hm, nqb, nterms) < 0) {
		fprintf(stderr, "b-circ: build_hamil failed\n");
		return;
	}
	struct circ_muldet md;
	if (build_muldet(&md, nqb, ndets) < 0) {
		fprintf(stderr, "b-circ: build_muldet failed\n");
		circ_hamil_free(&hm);
		return;
	}

	struct circ ct;
	if (circ_init(&ct, hm, STPREP_MULTIDET, &md, 1) < 0) {
		fprintf(stderr, "b-circ: circ_init failed\n");
		return; /* hm / md adopted; circ_init frees on error */
	}
	if (circ_prepst(&ct) < 0) {
		fprintf(stderr, "b-circ: circ_prepst failed\n");
		circ_free(&ct);
		return;
	}

	/* Warm-up: untimed step brings the cache + page state
	 * into a steady configuration. */
	circ_step(&ct, &ct.hm, 0.01);

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++) {
		BENCH_SAMPLE_MIN(samples[r], K,
			circ_step(&ct, &ct.hm, 0.01));
	}

	char params[80];
	snprintf(params, sizeof params,
		"{\"nqb\":%u,\"nterms\":%zu,\"ndets\":%zu}",
		nqb, nterms, ndets);
	/* nqb is constant for this binary; the banner shows it.
	 * Per-row label carries the varying params only. */
	(void)nqb;
	char display[40];
	snprintf(display, sizeof display,
		"circ_step t=%zu d=%zu", nterms, ndets);
	record(out, prov, mpi_ranks, "circ_step", params, display,
		samples, NUM_RUNS, K);

	circ_free(&ct);
}

static void measure_circ_measure(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, uint32_t nqb, size_t nterms, size_t ndets, int K)
{
	struct circ_hamil hm;
	if (build_hamil(&hm, nqb, nterms) < 0) {
		fprintf(stderr, "b-circ: build_hamil failed\n");
		return;
	}
	struct circ_muldet md;
	if (build_muldet(&md, nqb, ndets) < 0) {
		fprintf(stderr, "b-circ: build_muldet failed\n");
		circ_hamil_free(&hm);
		return;
	}

	struct circ ct;
	if (circ_init(&ct, hm, STPREP_MULTIDET, &md, 1) < 0) {
		fprintf(stderr, "b-circ: circ_init failed\n");
		return;
	}
	if (circ_prepst(&ct) < 0) {
		fprintf(stderr, "b-circ: circ_prepst failed\n");
		circ_free(&ct);
		return;
	}

	/* One step so the register isn't a trivial basis state. */
	circ_step(&ct, &ct.hm, 0.01);

	/* Warm-up. */
	g_sink_z = circ_measure(&ct);

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++) {
		BENCH_SAMPLE_MIN(samples[r], K,
			g_sink_z = circ_measure(&ct));
	}

	char params[64];
	snprintf(params, sizeof params,
		"{\"nqb\":%u,\"ndets\":%zu}", nqb, ndets);
	(void)nqb;
	char display[40];
	snprintf(display, sizeof display,
		"circ_measure d=%zu", ndets);
	record(out, prov, mpi_ranks, "circ_measure", params, display,
		samples, NUM_RUNS, K);

	circ_free(&ct);
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
			fprintf(stderr, "b-circ: cannot create bench/runs/\n");
			world_free();
			return 1;
		}
		out = fopen(path, "a");
		if (!out) {
			fprintf(stderr, "b-circ: cannot open %s\n", path);
			world_free();
			return 1;
		}
		bench_print_banner("b-circ (nqb=12)");
		bench_print_header();
	}

	/* qreg requires nqb >= log2(world_size). */
	uint32_t nqb_min = 1;
	int sz = wd.size;
	while (sz >>= 1)
		nqb_min++;

	const uint32_t nqb = 12;
	if (nqb >= nqb_min) {
		measure_circ_step   (out, &prov, wd.size, nqb,  10, 10, 10);
		measure_circ_step   (out, &prov, wd.size, nqb, 100, 10, 10);
		measure_circ_measure(out, &prov, wd.size, nqb,  10, 10, 100);
	}

	if (out)
		fclose(out);
	world_free();
	return 0;
}
