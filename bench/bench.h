/*
 * bench/bench.h -- header-only micro-benchmark
 * framework for phase2.
 *
 * Each bench binary (bench/b-NAME.c) is a
 * standalone program that loops over scenarios.
 * For each scenario it:
 *
 *   1. gathers NUM_RUNS samples by timing a tight
 *      inner loop of INNER_REPS op invocations;
 *   2. computes median / min / max plus a "noisy"
 *      flag ((max - min) / median > 15%);
 *   3. looks up a baseline from the host's JSONL
 *      file in bench/runs/<hostname>.jsonl (if
 *      any) and prints a table row with the
 *      percent delta;
 *   4. appends the new record to that file.
 *
 * The JSONL schema (one record per line):
 *
 *     {
 *       "timestamp": "2026-05-16T13:42:30Z",
 *       "hostname":  "...",
 *       "commit":    "...",
 *       "compiler":  "...",
 *       "backend":   "qreg" | "cuda",
 *       "mpi_ranks": <int>,
 *       "scenario":  "...",
 *       "params":    { ... },
 *       "num_runs":  <int>,
 *       "median_ns": <double>,
 *       "min_ns":    <double>,
 *       "max_ns":    <double>,
 *       "noisy":     <bool>
 *     }
 *
 * Baseline lookup matches both (hostname,
 * scenario) and the exact params JSON object so
 * records taken at different sizes / batch counts
 * stay distinct.  jsmn (include/jsmn.h) parses the
 * candidate lines.
 *
 * Provenance hostname / commit / timestamp /
 * compiler are auto-collected via gethostname,
 * popen("git rev-parse"), strftime, __VERSION__.
 *
 * Header-only: every helper is static inline.
 * Each bench binary inlines what it uses; the
 * framework leaves no .o file behind.
 */
#ifndef PHASE2_BENCH_H
#define PHASE2_BENCH_H

#define _POSIX_C_SOURCE 200809L
#define _DEFAULT_SOURCE         /* for timegm() */
#define _GNU_SOURCE             /* for sched_setaffinity() */

#include "c23_compat.h"
#include <errno.h>
#include <math.h>
#include <sched.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#include "jsmn.h"


/* -- timing ----------------------------------------------------------------*/

static inline double bench_now_ns(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}


/* -- CPU pinning -----------------------------------------------------------*/

/*
 * Pin the current thread to one fixed CPU so the OS scheduler doesn't
 * migrate it mid-sample (a migration changes which L2 cache the thread
 * sees and is a common source of variance).  Returns 0 on success, -1
 * otherwise.  Best-effort: a sandbox or container that strips
 * CAP_SYS_NICE may forbid the call; we don't make it fatal.
 */
static inline int bench_pin_cpu(int cpu)
{
	cpu_set_t set;
	CPU_ZERO(&set);
	CPU_SET(cpu, &set);
	return sched_setaffinity(0, sizeof set, &set);
}


/* -- stats -----------------------------------------------------------------*/

struct bench_stats {
	double median, min, max;
	bool noisy; /* (max - min) / median > 0.15 */
};

static inline int bench_cmp_double(const void *a, const void *b)
{
	const double da = *(const double *)a;
	const double db = *(const double *)b;
	return (da > db) - (da < db);
}

static inline struct bench_stats bench_compute_stats(double *data, int n)
{
	struct bench_stats s;
	qsort(data, (size_t)n, sizeof(double), bench_cmp_double);
	s.min = data[0];
	s.max = data[n - 1];
	s.median = (n % 2) ? data[n / 2] : 0.5 * (data[n / 2 - 1] + data[n / 2]);
	s.noisy = s.median > 0.0 && (s.max - s.min) / s.median > 0.15;
	return s;
}


/* -- provenance ------------------------------------------------------------*/

struct bench_prov {
	char hostname[128];
	char commit[64];
	char timestamp[32];
	char compiler[32];
};

static inline void bench_prov_init(struct bench_prov *p)
{
	if (gethostname(p->hostname, sizeof p->hostname) != 0)
		snprintf(p->hostname, sizeof p->hostname, "unknown");
	p->hostname[sizeof p->hostname - 1] = '\0';

	FILE *fp = popen("git rev-parse --short HEAD 2>/dev/null", "r");
	if (fp) {
		if (!fgets(p->commit, (int)sizeof p->commit, fp))
			snprintf(p->commit, sizeof p->commit, "unknown");
		else
			p->commit[strcspn(p->commit, "\n")] = '\0';
		pclose(fp);
	} else {
		snprintf(p->commit, sizeof p->commit, "unknown");
	}

	const time_t t = time(NULL);
	struct tm tm;
	gmtime_r(&t, &tm);
	strftime(p->timestamp, sizeof p->timestamp, "%Y-%m-%dT%H:%M:%SZ", &tm);

#ifdef __VERSION__
	snprintf(p->compiler, sizeof p->compiler, "%s", __VERSION__);
#else
	snprintf(p->compiler, sizeof p->compiler, "unknown");
#endif
}


/* -- output path -----------------------------------------------------------*/

/* Place the host's JSONL file at bench/runs/<hostname>.jsonl, relative to
 * the project root (the caller's cwd when invoked via `make bench`).
 * Returns 0 on success, -1 if mkdir or path-formatting fails. */
static inline int bench_runs_path(const struct bench_prov *p,
	char *out, size_t out_sz)
{
	if (mkdir("bench/runs", 0755) != 0 && errno != EEXIST)
		return -1;
	const int n = snprintf(out, out_sz, "bench/runs/%s.jsonl",
		p->hostname);
	return (n > 0 && (size_t)n < out_sz) ? 0 : -1;
}


/* -- baseline lookup -------------------------------------------------------*/

struct bench_baseline {
	double min_ns;
	char timestamp[32];
};

/* Match a jsmn STRING token against a literal C string. */
static inline bool bench_tok_eq(const char *js, const jsmntok_t *t,
	const char *s)
{
	if (t->type != JSMN_STRING)
		return false;
	const size_t n = (size_t)(t->end - t->start);
	return strlen(s) == n && strncmp(js + t->start, s, n) == 0;
}

/* Compare a JSON string-token's contents against a literal C string. */
static inline bool bench_slice_eq(const char *js, const jsmntok_t *t,
	const char *s)
{
	const size_t n = (size_t)(t->end - t->start);
	return strlen(s) == n && strncmp(js + t->start, s, n) == 0;
}

/*
 * Scan a JSONL file for the LAST record whose (hostname, scenario, params)
 * triple matches the arguments.  Returns true and populates *out if found.
 *
 * params_json must be the byte-for-byte canonical form the writer emitted
 * (no whitespace, fixed key order).
 */
static inline bool bench_find_baseline(const char *path,
	const char *hostname, const char *scenario,
	const char *params_json, struct bench_baseline *out)
{
	FILE *fp = fopen(path, "r");
	if (!fp)
		return false;

	bool found = false;
	char line[2048];

	while (fgets(line, sizeof line, fp)) {
		jsmn_parser p;
		jsmntok_t toks[64];
		jsmn_init(&p);
		const int n = jsmn_parse(&p, line, strlen(line),
			toks, sizeof toks / sizeof toks[0]);
		if (n < 2 || toks[0].type != JSMN_OBJECT)
			continue;

		bool host_ok = false, scen_ok = false, params_ok = false;
		double min_v = 0.0;
		const char *ts_ptr = NULL;
		int ts_len = 0;

		/* Walk top-level key/value pairs; tokens after toks[0] are
		 * key, value, key, value, ... at the top level. */
		for (int i = 1; i < n - 1; i++) {
			const jsmntok_t *k = &toks[i];
			const jsmntok_t *v = &toks[i + 1];
			if (k->type != JSMN_STRING)
				continue;
			if (bench_tok_eq(line, k, "hostname")) {
				host_ok = bench_slice_eq(line, v, hostname);
			} else if (bench_tok_eq(line, k, "scenario")) {
				scen_ok = bench_slice_eq(line, v, scenario);
			} else if (bench_tok_eq(line, k, "params")) {
				/* params is an OBJECT; compare its byte
				 * slice (from `{` through matching `}`)
				 * against the canonical params_json the
				 * writer emits. */
				if (v->type == JSMN_OBJECT) {
					const size_t span =
						(size_t)(v->end - v->start);
					params_ok =
						strlen(params_json) == span
						&& strncmp(line + v->start,
							params_json, span)
							== 0;
				}
			} else if (bench_tok_eq(line, k, "min_ns")) {
				char buf[64];
				const size_t len =
					(size_t)(v->end - v->start);
				if (len < sizeof buf) {
					memcpy(buf, line + v->start, len);
					buf[len] = '\0';
					min_v = strtod(buf, NULL);
				}
			} else if (bench_tok_eq(line, k, "timestamp")) {
				ts_ptr = line + v->start;
				ts_len = v->end - v->start;
			}
		}

		if (host_ok && scen_ok && params_ok) {
			out->min_ns = min_v;
			if (ts_ptr && ts_len > 0 &&
				(size_t)ts_len < sizeof out->timestamp) {
				memcpy(out->timestamp, ts_ptr,
					(size_t)ts_len);
				out->timestamp[ts_len] = '\0';
			} else {
				out->timestamp[0] = '\0';
			}
			found = true;
		}
	}
	fclose(fp);
	return found;
}

/* Baseline is "stale" if older than `days` days.  Used to flag the row in
 * the console table when a recorded baseline is unlikely to reflect the
 * current toolchain / hardware. */
static inline bool bench_timestamp_stale(const char *ts, int days)
{
	struct tm tm;
	memset(&tm, 0, sizeof tm);
	if (sscanf(ts, "%d-%d-%dT%d:%d:%d",
		    &tm.tm_year, &tm.tm_mon, &tm.tm_mday,
		    &tm.tm_hour, &tm.tm_min, &tm.tm_sec) != 6)
		return false;
	tm.tm_year -= 1900;
	tm.tm_mon -= 1;
	const time_t then = timegm(&tm);
	return difftime(time(NULL), then) > (double)days * 86400.0;
}


/* -- JSONL output ----------------------------------------------------------*/

static inline void bench_append_jsonl(FILE *fp, const struct bench_prov *prov,
	const char *backend, int mpi_ranks,
	const char *scenario, const char *params_json,
	int num_runs, const struct bench_stats *st)
{
	fprintf(fp,
		"{\"timestamp\":\"%s\","
		"\"hostname\":\"%s\","
		"\"commit\":\"%s\","
		"\"compiler\":\"%s\","
		"\"backend\":\"%s\","
		"\"mpi_ranks\":%d,"
		"\"scenario\":\"%s\","
		"\"params\":%s,"
		"\"num_runs\":%d,"
		"\"median_ns\":%.2f,"
		"\"min_ns\":%.2f,"
		"\"max_ns\":%.2f,"
		"\"noisy\":%s}\n",
		prov->timestamp, prov->hostname, prov->commit, prov->compiler,
		backend, mpi_ranks, scenario, params_json, num_runs,
		st->median, st->min, st->max,
		st->noisy ? "true" : "false");
	fflush(fp);
}


/* -- console table ---------------------------------------------------------*/

static inline void bench_print_header(void)
{
	printf("%-32s %12s %12s %12s %12s %8s\n",
		"scenario / params",
		"min (ns)", "median (ns)", "max (ns)", "prev (ns)", "delta");
	printf("%-32s %12s %12s %12s %12s %8s\n",
		"-----------------",
		"--------", "-----------", "--------", "---------", "-----");
}

/*
 * The headline is `min` (delta = (st->min - bl->min_ns) / bl->min_ns).
 * Noise at the nanosecond scale is additive -- interrupts, cache
 * evictions, frequency dips only ever make a call slower -- so the
 * minimum is the best estimator of the unperturbed cost.  Median and
 * max stay in the table for context.
 */
static inline void bench_print_row(const char *label,
	const struct bench_stats *st, const struct bench_baseline *bl,
	bool has_baseline)
{
	char delta[16] = "      --";
	char warn[24] = "";

	if (has_baseline && bl->min_ns > 0.0) {
		const double pct = (st->min - bl->min_ns)
			/ bl->min_ns * 100.0;
		snprintf(delta, sizeof delta, "%+6.1f%%", pct);
		if (bench_timestamp_stale(bl->timestamp, 30))
			snprintf(warn, sizeof warn, " [stale]");
	}

	printf("%-32s %12.2f %12.2f %12.2f %12.2f %8s%s%s\n",
		label, st->min, st->median, st->max,
		has_baseline ? bl->min_ns : 0.0,
		delta,
		st->noisy ? " [noisy]" : "", warn);
}


/* -- backend tag (compile-time) --------------------------------------------*/

#if PHASE2_BACKEND == 0
#define BENCH_BACKEND "qreg"
#elif PHASE2_BACKEND == 2
#define BENCH_BACKEND "cuda"
#else
#define BENCH_BACKEND "unknown"
#endif


/* -- sample-gathering loop -------------------------------------------------*/

/*
 * Gather NUM_RUNS samples by timing a tight loop of INNER_REPS
 * invocations of the op, returning per-op nanoseconds in samples[].
 * The op is called once with `warmup_arg` before the timed loop fires;
 * pass the same pointer for both warm-up and timed paths to keep cache
 * state realistic.
 *
 * Defined as a macro so the op call inlines at the bench site.
 */
#define BENCH_SAMPLE_LOOP(samples, num_runs, inner_reps, op_call)	\
	do {								\
		for (int _br = 0; _br < (num_runs); _br++) {		\
			const double _t0 = bench_now_ns();		\
			for (int _bi = 0; _bi < (inner_reps); _bi++) {	\
				op_call;				\
			}						\
			const double _t1 = bench_now_ns();		\
			(samples)[_br] = (_t1 - _t0) / (double)(inner_reps); \
		}							\
	} while (0)


#endif /* PHASE2_BENCH_H */
