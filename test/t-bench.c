/*
 * t-bench -- unit tests for the bench framework
 * (bench/bench.h).
 *
 * Covers:
 *   bench_compute_stats: median / min / max, 15%
 *     noisy threshold in both directions, even-count
 *     median.
 *   bench_now_ns: monotone non-decreasing.
 *   bench_pin_cpu: best-effort; returns 0 or -1.
 *   bench_runs_path: path layout + mkdir idempotence
 *     in a tmp cwd.
 *   bench_timestamp_stale: 30-day boundary using a
 *     fixed-past anchor and an isoformatted "now".
 *   bench_find_baseline: hostname / scenario / params
 *     / sub_samples match keys, last-record wins,
 *     pre-MOM records never match, missing file ->
 *     false.
 *   bench_color_on + print helpers: no ANSI escape
 *     codes when stdout isn't a TTY; row carries
 *     label, K, delta, [noisy] flag.
 *
 * Stdout is captured via a self-pipe (dup2) so the
 * binary works whether the runner pipes its stdout or
 * not.  The first capture also locks
 * bench_color_on's cached state to "off" by ensuring
 * the first isatty() check sees a pipe, not a TTY.
 */

/* bench.h sets _GNU_SOURCE etc.; include it before any
 * system header so feature macros take effect. */
#include "bench.h"

#include "c23_compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#include "phase2/world.h"

#include "test.h"

#define WD_SEED   UINT64_C(0x9c3a5b8e4f7d2106)
#define CAP_BYTES (16 * 1024)


/* -- stdout capture --------------------------------------------------------*/

struct capture {
	int target_fd;
	int saved_fd;
	int pipe_r;
	int pipe_w;
};

static int capture_start(struct capture *c, int target_fd)
{
	int pfd[2];
	if (pipe(pfd) < 0)
		return -1;
	c->pipe_r = pfd[0];
	c->pipe_w = pfd[1];
	c->target_fd = target_fd;
	c->saved_fd = dup(target_fd);
	if (c->saved_fd < 0)
		goto err;
	if (dup2(c->pipe_w, target_fd) < 0)
		goto err;
	return 0;
err:
	close(c->pipe_r);
	close(c->pipe_w);
	return -1;
}

static int capture_end(struct capture *c, char *buf, size_t cap)
{
	if (c->target_fd == STDOUT_FILENO)
		fflush(stdout);
	else if (c->target_fd == STDERR_FILENO)
		fflush(stderr);

	if (dup2(c->saved_fd, c->target_fd) < 0)
		return -1;
	close(c->saved_fd);
	close(c->pipe_w);

	size_t total = 0;
	while (total + 1 < cap) {
		ssize_t n = read(c->pipe_r, buf + total, cap - 1 - total);
		if (n <= 0)
			break;
		total += (size_t)n;
	}
	buf[total] = '\0';
	close(c->pipe_r);
	return (int)total;
}


/* -- temp file helpers -----------------------------------------------------*/

static void write_file(const char *path, const char *content)
{
	FILE *fp = fopen(path, "w");
	TEST_ASSERT(fp != NULL, "fopen %s", path);
	fputs(content, fp);
	fclose(fp);
}


/* -- tests: stats ----------------------------------------------------------*/

static void t_stats_basic(void)
{
	double a[] = { 3.0, 1.0, 2.0 };
	struct bench_stats s = bench_compute_stats(a, 3);
	TEST_NEAR(s.min, 1.0, 1e-12);
	TEST_NEAR(s.max, 3.0, 1e-12);
	TEST_NEAR(s.median, 2.0, 1e-12);
	/* spread (3-1)/median 2 = 1.0 > 0.15 -> noisy */
	TEST_EQ(s.noisy, true);
}

static void t_stats_even_count(void)
{
	double a[] = { 4.0, 2.0, 3.0, 1.0 };
	struct bench_stats s = bench_compute_stats(a, 4);
	TEST_NEAR(s.median, 2.5, 1e-12);
}

static void t_stats_quiet(void)
{
	double a[] = { 100.0, 100.0, 100.0 };
	struct bench_stats s = bench_compute_stats(a, 3);
	TEST_EQ(s.noisy, false);
	TEST_NEAR(s.median, 100.0, 1e-12);
}

static void t_stats_noisy_threshold(void)
{
	/* Median = 100; spread <= 15 not noisy, > 15 noisy. */
	double below[] = { 92.5, 100.0, 107.4 }; /* spread 14.9 */
	TEST_EQ(bench_compute_stats(below, 3).noisy, false);

	double above[] = { 92.0, 100.0, 108.0 }; /* spread 16.0 */
	TEST_EQ(bench_compute_stats(above, 3).noisy, true);
}


/* -- tests: timing and pinning --------------------------------------------*/

static void t_now_ns_monotone(void)
{
	const double t0 = bench_now_ns();
	const double t1 = bench_now_ns();
	TEST_ASSERT(t1 >= t0,
		"bench_now_ns went backwards: t0=%g t1=%g", t0, t1);
}

static void t_pin_cpu(void)
{
	const int r = bench_pin_cpu(0);
	TEST_ASSERT(r == 0 || r == -1,
		"bench_pin_cpu returned unexpected %d", r);
}


/* -- tests: runs_path ------------------------------------------------------*/

static void t_runs_path_basic(void)
{
	char tmpl[] = "/tmp/t-bench.XXXXXX";
	TEST_ASSERT(mkdtemp(tmpl) != NULL, "mkdtemp");
	char cwd[256];
	TEST_ASSERT(getcwd(cwd, sizeof cwd) != NULL, "getcwd");
	TEST_ASSERT(chdir(tmpl) == 0, "chdir %s", tmpl);

	/* bench_runs_path mkdirs "bench/runs"; "bench" itself
	 * must already exist (the bench binaries are run from
	 * the project root where it does).  Mirror that here. */
	TEST_ASSERT(mkdir("bench", 0755) == 0, "mkdir bench");

	struct bench_prov prov;
	memset(&prov, 0, sizeof prov);
	snprintf(prov.hostname, sizeof prov.hostname, "host42");

	char path[256];
	TEST_EQ(bench_runs_path(&prov, path, sizeof path), 0);
	TEST_ASSERT(strcmp(path, "bench/runs/host42.jsonl") == 0,
		"path mismatch: %s", path);

	/* Idempotent: second call ok even though bench/runs exists. */
	TEST_EQ(bench_runs_path(&prov, path, sizeof path), 0);

	rmdir("bench/runs");
	rmdir("bench");
	TEST_ASSERT(chdir(cwd) == 0, "chdir back");
	rmdir(tmpl);
}


/* -- tests: timestamp_stale ------------------------------------------------*/

static void t_timestamp_stale_far_past(void)
{
	TEST_EQ(bench_timestamp_stale("2000-01-01T00:00:00", 1), true);
}

static void t_timestamp_stale_now(void)
{
	char now_iso[32];
	const time_t t = time(NULL);
	struct tm tm;
	gmtime_r(&t, &tm);
	strftime(now_iso, sizeof now_iso, "%Y-%m-%dT%H:%M:%S", &tm);
	TEST_EQ(bench_timestamp_stale(now_iso, 30), false);
}


/* -- tests: find_baseline --------------------------------------------------*/

static void t_find_baseline_match(void)
{
	char tmpl[] = "/tmp/t-bench.XXXXXX";
	TEST_ASSERT(mkdtemp(tmpl) != NULL, "mkdtemp");
	char path[256];
	snprintf(path, sizeof path, "%s/test.jsonl", tmpl);

	write_file(path,
		"{\"timestamp\":\"2026-05-15T00:00:00Z\","
		"\"hostname\":\"host42\",\"scenario\":\"paulis_set\","
		"\"params\":{\"nqb\":64},"
		"\"num_runs\":11,\"sub_samples\":1000,"
		"\"median_ns\":1.50,\"min_ns\":1.45,\"max_ns\":1.60,"
		"\"noisy\":false}\n"
		"{\"timestamp\":\"2026-05-16T00:00:00Z\","
		"\"hostname\":\"host42\",\"scenario\":\"paulis_get\","
		"\"params\":{\"nqb\":64},"
		"\"num_runs\":11,\"sub_samples\":1000,"
		"\"median_ns\":2.00,\"min_ns\":1.90,\"max_ns\":2.10,"
		"\"noisy\":false}\n");

	struct bench_baseline bl;
	TEST_EQ(bench_find_baseline(path, "host42", "paulis_set",
		"{\"nqb\":64}", 1000, &bl), true);
	TEST_NEAR(bl.min_ns, 1.45, 1e-12);
	TEST_ASSERT(strcmp(bl.timestamp, "2026-05-15T00:00:00Z") == 0,
		"ts %s", bl.timestamp);

	/* Wrong K -> miss. */
	TEST_EQ(bench_find_baseline(path, "host42", "paulis_set",
		"{\"nqb\":64}", 100, &bl), false);

	/* Wrong host -> miss. */
	TEST_EQ(bench_find_baseline(path, "other", "paulis_set",
		"{\"nqb\":64}", 1000, &bl), false);

	/* Wrong params -> miss. */
	TEST_EQ(bench_find_baseline(path, "host42", "paulis_set",
		"{\"nqb\":32}", 1000, &bl), false);

	unlink(path);
	rmdir(tmpl);
}

static void t_find_baseline_last_wins(void)
{
	char tmpl[] = "/tmp/t-bench.XXXXXX";
	TEST_ASSERT(mkdtemp(tmpl) != NULL, "mkdtemp");
	char path[256];
	snprintf(path, sizeof path, "%s/test.jsonl", tmpl);

	write_file(path,
		"{\"timestamp\":\"2026-01-01T00:00:00Z\","
		"\"hostname\":\"h\",\"scenario\":\"s\","
		"\"params\":{},\"num_runs\":1,\"sub_samples\":1,"
		"\"median_ns\":10,\"min_ns\":10,\"max_ns\":10,"
		"\"noisy\":false}\n"
		"{\"timestamp\":\"2026-05-01T00:00:00Z\","
		"\"hostname\":\"h\",\"scenario\":\"s\","
		"\"params\":{},\"num_runs\":1,\"sub_samples\":1,"
		"\"median_ns\":20,\"min_ns\":20,\"max_ns\":20,"
		"\"noisy\":false}\n");

	struct bench_baseline bl;
	TEST_EQ(bench_find_baseline(path, "h", "s", "{}", 1, &bl), true);
	TEST_NEAR(bl.min_ns, 20.0, 1e-12);
	TEST_ASSERT(strcmp(bl.timestamp, "2026-05-01T00:00:00Z") == 0,
		"ts %s", bl.timestamp);

	unlink(path);
	rmdir(tmpl);
}

static void t_find_baseline_pre_mom(void)
{
	/* Old record without sub_samples must not match -- the
	 * MOM rework intentionally invalidates pre-K records. */
	char tmpl[] = "/tmp/t-bench.XXXXXX";
	TEST_ASSERT(mkdtemp(tmpl) != NULL, "mkdtemp");
	char path[256];
	snprintf(path, sizeof path, "%s/test.jsonl", tmpl);

	write_file(path,
		"{\"timestamp\":\"2026-05-15T00:00:00Z\","
		"\"hostname\":\"h\",\"scenario\":\"s\","
		"\"params\":{},\"num_runs\":1,"
		"\"median_ns\":1,\"min_ns\":1,\"max_ns\":1,"
		"\"noisy\":false}\n");

	struct bench_baseline bl;
	TEST_EQ(bench_find_baseline(path, "h", "s", "{}", 1, &bl), false);

	unlink(path);
	rmdir(tmpl);
}

static void t_find_baseline_missing_file(void)
{
	struct bench_baseline bl;
	TEST_EQ(bench_find_baseline("/tmp/t-bench-nonexistent.jsonl",
		"h", "s", "{}", 1, &bl), false);
}


/* -- tests: colour gating and print helpers --------------------------------*/

/*
 * Lock bench_color_on's cache to "off" by triggering its first
 * isatty(STDOUT_FILENO) check while stdout is dup2'd to a pipe.
 * After this, every later test sees the same cached value
 * regardless of how the binary's stdout is connected.
 */
static void prime_color_off(void)
{
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	(void)bench_color_on();
	char drain[64];
	capture_end(&cap, drain, sizeof drain);
}

static void t_color_off(void)
{
	TEST_EQ(bench_color_on(), false);
}

static void t_print_header_no_escapes(void)
{
	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	bench_print_header();
	const int n = capture_end(&cap, buf, sizeof buf);

	TEST_ASSERT(n > 0, "captured nothing");
	TEST_ASSERT(strchr(buf, '\x1b') == NULL,
		"ANSI escape leaked into non-TTY output: %s", buf);
	TEST_ASSERT(strstr(buf, "scenario") != NULL,
		"missing 'scenario' header: %s", buf);
	TEST_ASSERT(strstr(buf, "min(ns)") != NULL,
		"missing 'min(ns)' header: %s", buf);
}

static void t_print_banner_no_escapes(void)
{
	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	bench_print_banner("b-test");
	capture_end(&cap, buf, sizeof buf);

	TEST_ASSERT(strchr(buf, '\x1b') == NULL,
		"banner leaked ANSI: %s", buf);
	TEST_ASSERT(strstr(buf, "b-test") != NULL,
		"banner missing name: %s", buf);
	TEST_ASSERT(strstr(buf, "==") != NULL,
		"banner missing rule: %s", buf);
}

static void t_print_row_baseline(void)
{
	struct bench_stats st = {
		.median = 1.50, .min = 1.45, .max = 1.55,
		.noisy = false,
	};
	struct bench_baseline bl = {
		.min_ns = 1.50,
	};
	snprintf(bl.timestamp, sizeof bl.timestamp,
		"2026-05-15T00:00:00Z");

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	bench_print_row("paulis_set", 1000, &st, &bl, true);
	capture_end(&cap, buf, sizeof buf);

	TEST_ASSERT(strstr(buf, "paulis_set") != NULL,
		"row missing label: %s", buf);
	TEST_ASSERT(strstr(buf, "1000") != NULL,
		"row missing K: %s", buf);
	/* delta = (1.45 - 1.50) / 1.50 * 100 = -3.3% */
	TEST_ASSERT(strstr(buf, "-3.3%") != NULL,
		"row missing delta: %s", buf);
	TEST_ASSERT(strstr(buf, "[noisy]") == NULL,
		"row spuriously noisy: %s", buf);
}

static void t_print_row_noisy_no_baseline(void)
{
	struct bench_stats st = {
		.median = 1.50, .min = 1.20, .max = 2.00,
		.noisy = true,
	};
	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	bench_print_row("paulis_set", 1000, &st, NULL, false);
	capture_end(&cap, buf, sizeof buf);

	TEST_ASSERT(strstr(buf, "[noisy]") != NULL,
		"row missing [noisy]: %s", buf);
	TEST_ASSERT(strstr(buf, "--") != NULL,
		"row missing '--' placeholder: %s", buf);
}


/* -- main ------------------------------------------------------------------*/

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	prime_color_off();
	t_color_off();

	t_stats_basic();
	t_stats_even_count();
	t_stats_quiet();
	t_stats_noisy_threshold();
	t_now_ns_monotone();
	t_pin_cpu();
	t_runs_path_basic();
	t_timestamp_stale_far_past();
	t_timestamp_stale_now();
	t_find_baseline_match();
	t_find_baseline_last_wins();
	t_find_baseline_pre_mom();
	t_find_baseline_missing_file();
	t_print_banner_no_escapes();
	t_print_header_no_escapes();
	t_print_row_baseline();
	t_print_row_noisy_no_baseline();

	world_free();
	return 0;
}
