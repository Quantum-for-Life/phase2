/*
 * b-paulis -- microbenchmarks for the bit-packed
 * Pauli-string encoding (phase2/paulis.h).
 *
 * Three scenarios: set / get / effect.  Each runs
 * NUM_RUNS samples of INNER_REPS tight-loop calls
 * and emits a JSONL record to bench/runs/<host>.jsonl.
 * Total bench time on a quiet host: ~0.5 s.
 */

#include "bench.h"

#include <complex.h>

#include "phase2/paulis.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#define WD_SEED  UINT64_C(0x0094be8e2d4ba8eb)
#define SEED     UINT64_C(0xa208312e4065b1df)

#define NUM_RUNS    11
#define INNER_REPS  100000
#define MOM_K       1000
#define NQB         64

static struct xoshiro256ss RNG;

/* Sink to defeat dead-code elimination on operations whose return value the
 * compiler would otherwise prove unused.  volatile so writes aren't sunk. */
static volatile uint64_t g_sink_u64;
static volatile int g_sink_int;
static volatile _Complex double g_sink_z;

/* Print a banner and the canonical params object so the operator can see
 * which scenario set is being run.  The params string is the EXACT byte
 * sequence the writer emits and the reader matches against. */
#define PARAMS_NQB64 "{\"nqb\":64}"

static void run_scenario(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, const char *name, const char *params_json,
	int sub_samples, double *samples)
{
	const struct bench_stats st = bench_compute_stats(samples, NUM_RUNS);

	char path[256];
	bench_runs_path(prov, path, sizeof path);

	struct bench_baseline bl;
	const bool has_bl = bench_find_baseline(path, prov->hostname,
		name, params_json, sub_samples, &bl);

	/* All b-paulis scenarios share nqb=64; the banner shows it once,
	 * so the row label is just the op name. */
	bench_print_row(name, sub_samples, &st, &bl, has_bl);

	bench_append_jsonl(out, prov, BENCH_BACKEND, mpi_ranks,
		name, params_json, NUM_RUNS, sub_samples, &st);
}

int main(void)
{
	/* Quiet phase2's INFO banners by default so the bench table is
	 * the only thing on stdout.  Don't overwrite an existing setting
	 * -- `PHASE2_LOG=debug make bench` still works. */
	setenv("PHASE2_LOG", "warn", 0);

	/* Pin to CPU 0 so the OS scheduler doesn't migrate us mid-sample.
	 * Best-effort: failures are ignored (containers may forbid it). */
	bench_pin_cpu(0);

	world_init(nullptr, nullptr, WD_SEED);
	struct world_info wd;
	world_info(&wd);
	xoshiro256ss_init(&RNG, SEED);

	struct bench_prov prov;
	bench_prov_init(&prov);

	char path[256];
	FILE *out = NULL;
	if (wd.rank == 0) {
		if (bench_runs_path(&prov, path, sizeof path) < 0) {
			fprintf(stderr,
				"b-paulis: cannot create bench/runs/\n");
			world_free();
			return 1;
		}
		out = fopen(path, "a");
		if (!out) {
			fprintf(stderr,
				"b-paulis: cannot open %s\n", path);
			world_free();
			return 1;
		}
		bench_print_banner("b-paulis (nqb=64)");
		bench_print_header();
	}

	double samples[NUM_RUNS];

	/* Set up a random Pauli string for the get / effect benches.  Use the
	 * same string for both so cache state is consistent between
	 * scenarios. */
	struct paulis ps = paulis_new();
	for (uint32_t k = 0; k < NQB; k++)
		paulis_set(&ps, (enum pauli_op)(xoshiro256ss_next(&RNG) % 4),
			k);

	/* --- paulis_set -- toggle every qubit, cycling through ops --- */
	{
		struct paulis p = paulis_new();
		BENCH_SAMPLE_LOOP(samples, NUM_RUNS, MOM_K, INNER_REPS, ({
			paulis_set(&p, (enum pauli_op)(_bi & 0x3),
				(uint32_t)(_bi % NQB));
		}));
		if (out)
			run_scenario(out, &prov, wd.size, "paulis_set",
				PARAMS_NQB64, MOM_K, samples);
	}

	/* --- paulis_get -- read every qubit --- */
	{
		BENCH_SAMPLE_LOOP(samples, NUM_RUNS, MOM_K, INNER_REPS, ({
			g_sink_int = paulis_get(ps, (uint32_t)(_bi % NQB));
		}));
		if (out)
			run_scenario(out, &prov, wd.size, "paulis_get",
				PARAMS_NQB64, MOM_K, samples);
	}

	/* --- paulis_effect -- per-amplitude rotation kernel --- */
	{
		_Complex double z = 1.0;
		uint64_t i = 0x123456789ABCDEFULL;
		BENCH_SAMPLE_LOOP(samples, NUM_RUNS, MOM_K, INNER_REPS, ({
			i = paulis_effect(ps, i, &z);
		}));
		g_sink_u64 = i;
		g_sink_z = z;
		if (out)
			run_scenario(out, &prov, wd.size, "paulis_effect",
				PARAMS_NQB64, MOM_K, samples);
	}

	if (out)
		fclose(out);
	world_free();
	return 0;
}
