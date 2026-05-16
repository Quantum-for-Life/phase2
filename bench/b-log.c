/*
 * b-log -- micro-benchmark for the log_at gating
 * macro (include/log.h).
 *
 * log_at expands to a branch on log_threshold
 * followed by a call into log_emit when taken.
 * Two scenarios:
 *
 *   filtered  level < log_threshold; the branch
 *             is not taken, log_emit is not
 *             called.  Headline cost ~1 ns
 *             (predicted branch + one load).
 *             This is the common production case
 *             for log_trace / log_debug.
 *
 *   emitted   level >= log_threshold; log_emit
 *             runs (snprintf + clock_gettime +
 *             localtime_r + fwrite).  Upper
 *             bound on the cost of an enabled
 *             log line.
 *
 * The emitted scenario redirects stderr to
 * /dev/null around the timed loop so I/O cost
 * is bounded and the bench table on stdout
 * stays clean.
 *
 * Total bench time on a quiet host: ~2 s.
 */

#define LOG_SUBSYS "b-log"

#include "bench.h"

#include <fcntl.h>

#include "log.h"
#include "phase2/world.h"

#define WD_SEED          UINT64_C(0x4f6b9c3d8e1a720d)

#define NUM_RUNS         11
#define MOM_K_FILTERED   1000
#define MOM_K_EMITTED    10
#define INNER_FILTERED   100000
#define INNER_EMITTED    1000

/* Sink so the compiler doesn't optimise away log_at when the
 * arguments include only constants and a loop variable. */
static volatile int g_sink;

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

	bench_print_row(name, sub_samples, &st, &bl, has_bl);

	bench_append_jsonl(out, prov, BENCH_BACKEND, mpi_ranks,
		name, params_json, NUM_RUNS, sub_samples, &st);
}

int main(void)
{
	/* world_init() invokes log_init() which honours PHASE2_LOG;
	 * keep that quiet so the bench table is the only thing on
	 * stdout. */
	setenv("PHASE2_LOG", "warn", 0);

	bench_pin_cpu(0);

	world_init(nullptr, nullptr, WD_SEED);
	struct world_info wd;
	world_info(&wd);

	struct bench_prov prov;
	bench_prov_init(&prov);

	FILE *out = NULL;
	if (wd.rank == 0) {
		char path[256];
		if (bench_runs_path(&prov, path, sizeof path) < 0) {
			fprintf(stderr, "b-log: cannot create bench/runs/\n");
			world_free();
			return 1;
		}
		out = fopen(path, "a");
		if (!out) {
			fprintf(stderr, "b-log: cannot open %s\n", path);
			world_free();
			return 1;
		}
		bench_print_banner("b-log");
		bench_print_header();
	}

	double samples[NUM_RUNS];

	/* --- filtered: log_at(LOG_TRACE) with threshold=LOG_INFO --- */
	{
		log_threshold = LOG_INFO;
		log_at(LOG_TRACE, "warmup %d", 0); /* warm-up */

		BENCH_SAMPLE_LOOP(samples, NUM_RUNS,
			MOM_K_FILTERED, INNER_FILTERED, ({
			log_at(LOG_TRACE, "msg %d", _bi);
			g_sink = _bi;
		}));
		if (out)
			run_scenario(out, &prov, wd.size, "log_at_filtered",
				"{}", MOM_K_FILTERED, samples);
	}

	/* --- emitted: log_at(LOG_WARN) with stderr -> /dev/null --- */
	{
		const int saved_stderr = dup(STDERR_FILENO);
		const int devnull = open("/dev/null", O_WRONLY);
		if (devnull >= 0) {
			dup2(devnull, STDERR_FILENO);
			close(devnull);
		}

		log_threshold = LOG_WARN;
		log_at(LOG_WARN, "warmup %d", 0); /* warm-up */

		BENCH_SAMPLE_LOOP(samples, NUM_RUNS,
			MOM_K_EMITTED, INNER_EMITTED, ({
			log_at(LOG_WARN, "msg %d", _bi);
			g_sink = _bi;
		}));

		/* Restore stderr before printing the row. */
		if (saved_stderr >= 0) {
			fflush(stderr);
			dup2(saved_stderr, STDERR_FILENO);
			close(saved_stderr);
		}

		if (out)
			run_scenario(out, &prov, wd.size, "log_at_emitted",
				"{}", MOM_K_EMITTED, samples);
	}

	if (out)
		fclose(out);
	world_free();
	return 0;
}
