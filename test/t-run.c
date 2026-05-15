/*
 * t-run -- self-tests for the parallel runner at
 * test/run.
 *
 * The runner has no compile-time hooks; this test
 * exercises it as a black box by invoking the built
 * `./test/run` binary via popen(3) against the synthetic
 * fixtures under `test/run-fixtures/`.  Each scenario
 * checks the runner's two observable contracts:
 *
 *   - its exit code (0 iff every dispatched test
 *     passed; 1 otherwise),
 *   - the surface of its stdout (status line, summary
 *     counts, failure dumps, ANSI suppression).
 *
 * Fixtures:
 *   pass    exits 0
 *   fail    exits 1
 *   sleep   usleep(250 ms), exits 0
 *   abort   abort() -- signal-killed (SIGABRT)
 *   banner  writes a distinctive line to stdout AND
 *           stderr, then exits 1
 *
 * The runner is invoked under `--no-color` for every
 * case except the dedicated colour test, so the
 * assertions can grep plain ASCII without escaping ANSI
 * sequences.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include "test.h"

/*
 * Runner binary, fixture binaries, and PH2_TESTDIR
 * (used by sibling tests) all live under build/.
 * t-run is invoked by the runner with cwd = repo
 * root (see test/Makefile's `check` target), so
 * relative paths under ./build/test work.
 */
#define TBIN       "./build/test"
#define FIXDIR     TBIN "/run-fixtures"
#define FIX_PASS   FIXDIR "/pass"
#define FIX_FAIL   FIXDIR "/fail"
#define FIX_SLEEP  FIXDIR "/sleep"
#define FIX_ABORT  FIXDIR "/abort"
#define FIX_BANNER FIXDIR "/banner"

#define RUN_BIN    TBIN "/run"


/* -- popen wrapper --------------------------------------------------- */

struct run_out {
	char buf[16384];	/* enough for any fixture run */
	size_t len;
	int exit_code;		/* normalised: WEXITSTATUS or
				 * 128+sig, -1 on harness error */
};

/*
 * Invoke `cmd` via popen, drain its stdout (with stderr
 * folded in via shell `2>&1`), and report the normalised
 * child exit status.  Used so every scenario in this
 * file can assert against a single captured stream.
 */
static struct run_out run_cmd(const char *cmd)
{
	struct run_out r = { 0 };
	char full[1024];
	const int n = snprintf(full, sizeof full, "%s 2>&1", cmd);
	TEST_ASSERT(n > 0 && (size_t)n < sizeof full,
		"command too long: %s", cmd);

	FILE *p = popen(full, "r");
	TEST_ASSERT(p != NULL, "popen failed: %s", cmd);

	size_t got;
	while (r.len + 1 < sizeof r.buf
	       && (got = fread(r.buf + r.len, 1,
				sizeof r.buf - 1 - r.len, p)) > 0)
		r.len += got;
	r.buf[r.len] = '\0';

	const int rc = pclose(p);
	if (rc == -1)
		r.exit_code = -1;
	else if (WIFEXITED(rc))
		r.exit_code = WEXITSTATUS(rc);
	else if (WIFSIGNALED(rc))
		r.exit_code = 128 + WTERMSIG(rc);
	else
		r.exit_code = -1;
	return r;
}

/*
 * Substring-match helper with a descriptive failure
 * message.  Embeds (a truncated head of) the captured
 * output so a failure makes the divergence obvious.
 */
static void must_contain(const struct run_out *r, const char *needle,
	const char *what)
{
	if (strstr(r->buf, needle) != NULL)
		return;
	fprintf(stderr, "expected %s (%s) in runner output, got:\n%.2000s\n",
		what, needle, r->buf);
	TEST_FAIL("substring %s not found", needle);
}

static void must_not_contain(const struct run_out *r, const char *needle,
	const char *what)
{
	if (strstr(r->buf, needle) == NULL)
		return;
	fprintf(stderr, "unexpected %s (%s) in runner output:\n%.2000s\n",
		what, needle, r->buf);
	TEST_FAIL("substring %s should not appear", needle);
}


/* -- scenarios ------------------------------------------------------- */

/*
 * Single passing test -> exit 0, "ok" verdict, "1
 * passed; 0 failed" totals.
 */
static void t_pass_one(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color -- " FIX_PASS);
	TEST_EQ(r.exit_code, 0);
	must_contain(&r, "test t-pass", "status line for pass fixture");
	must_contain(&r, " ok ", "ok verdict");
	must_contain(&r, "1 passed; 0 failed", "summary totals");
}

/*
 * Single failing test -> exit 1, "FAILED" verdict, "0
 * passed; 1 failed" totals.
 */
static void t_fail_one(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color -- " FIX_FAIL);
	TEST_EQ(r.exit_code, 1);
	must_contain(&r, "test t-fail", "status line for fail fixture");
	must_contain(&r, "FAILED", "FAILED verdict");
	must_contain(&r, "0 passed; 1 failed", "summary totals");
	must_contain(&r, "(exit 1)", "exit code annotation");
}

/*
 * Mixed run: pass + fail simultaneously -> exit 1, the
 * fail is itemised under the "failures:" block.
 */
static void t_mixed(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color --jobs=2 -- "
		FIX_PASS " " FIX_FAIL);
	TEST_EQ(r.exit_code, 1);
	must_contain(&r, "1 passed; 1 failed", "mixed summary totals");
	must_contain(&r, "failures:", "failures block header");
	must_contain(&r, "    t-fail", "fail listed in failures block");
}

/*
 * Signal-killed test (SIGABRT) -> runner reports exit
 * code 128+6 = 134.
 */
static void t_abort_signal(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color -- " FIX_ABORT);
	TEST_EQ(r.exit_code, 1);
	must_contain(&r, "FAILED", "abort yields FAILED verdict");
	must_contain(&r, "(exit 134)", "SIGABRT maps to exit 134");
}

/*
 * Failing test that wrote to both streams -> runner
 * dumps the captured output under labelled headers.
 */
static void t_failure_dump(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color -- " FIX_BANNER);
	TEST_EQ(r.exit_code, 1);
	must_contain(&r, "---- t-banner stdout ----", "stdout header");
	must_contain(&r, "BANNER_STDOUT_LINE", "stdout body captured");
	must_contain(&r, "---- t-banner stderr ----", "stderr header");
	must_contain(&r, "BANNER_STDERR_LINE", "stderr body captured");
}

/*
 * --filter narrows the set by effective name.  Pass +
 * fail supplied, filter='pass*' -> only pass runs ->
 * exit 0.
 */
static void t_filter(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color --filter='pass*' -- "
		FIX_PASS " " FIX_FAIL);
	TEST_EQ(r.exit_code, 0);
	must_contain(&r, "running 1 test", "filter narrows to 1");
	must_contain(&r, "1 passed; 0 failed", "filter summary totals");
	must_not_contain(&r, "t-fail", "fail filtered out");
}

/*
 * --jobs=4 with 4 sleep(250 ms) tests should finish
 * well under the 1.0 s serial bound.  Bound is loose
 * (0.8 s) to tolerate fork/exec overhead and CI load,
 * but tight enough to detect the runner falling back
 * to serial execution.
 */
static void t_parallel_speedup(void)
{
	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC, &t0);
	struct run_out r = run_cmd(RUN_BIN " --no-color --jobs=4 -- "
		FIX_SLEEP " " FIX_SLEEP " " FIX_SLEEP " " FIX_SLEEP);
	clock_gettime(CLOCK_MONOTONIC, &t1);

	TEST_EQ(r.exit_code, 0);
	must_contain(&r, "4 passed; 0 failed", "all four passed");

	const double elapsed = (double)(t1.tv_sec - t0.tv_sec)
		+ (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
	TEST_ASSERT(elapsed < 0.8,
		"expected parallel speedup (4x usleep(250ms) "
		"with --jobs=4), got %.2fs wall", elapsed);
}

/*
 * --jobs=1 with the same workload should take roughly
 * the serial time (>= 4 x 250 ms).  This is the
 * negative of t_parallel_speedup -- proves the bound
 * above is actually load-bearing.
 */
static void t_serial_no_speedup(void)
{
	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC, &t0);
	struct run_out r = run_cmd(RUN_BIN " --no-color --jobs=1 -- "
		FIX_SLEEP " " FIX_SLEEP " " FIX_SLEEP " " FIX_SLEEP);
	clock_gettime(CLOCK_MONOTONIC, &t1);

	TEST_EQ(r.exit_code, 0);
	const double elapsed = (double)(t1.tv_sec - t0.tv_sec)
		+ (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
	TEST_ASSERT(elapsed >= 0.9,
		"--jobs=1 should run serially, got %.2fs wall "
		"(expected >= 0.9s for 4x usleep(250ms))",
		elapsed);
}

/*
 * --no-color must not emit ANSI escape introducers.
 * The runner's status line is already plain ASCII --
 * the assertion guards against accidental introduction
 * of escapes in some future edit.
 */
static void t_no_color(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color -- " FIX_PASS);
	TEST_EQ(r.exit_code, 0);
	must_not_contain(&r, "\x1b[", "ANSI escape under --no-color");
}

/*
 * FORCE_COLOR=1 should re-enable colour even when the
 * stdout side of popen is a pipe (i.e. isatty false).
 * Confirms the runner honours the env override.
 */
static void t_force_color(void)
{
	struct run_out r = run_cmd("FORCE_COLOR=1 " RUN_BIN " -- " FIX_PASS);
	TEST_EQ(r.exit_code, 0);
	must_contain(&r, "\x1b[32m", "green ANSI under FORCE_COLOR");
}

/*
 * No paths after `--` is a usage error (exit 2).  Guard
 * against the runner silently exiting 0 on an empty
 * suite, which would mask CI invocation bugs.
 */
static void t_no_paths(void)
{
	struct run_out r = run_cmd(RUN_BIN " --no-color");
	TEST_EQ(r.exit_code, 2);
	must_contain(&r, "no test paths supplied", "empty-list diagnostic");
}


/* -- main ------------------------------------------------------------ */

/*
 * Detect a wrapping mpirun and skip: this meta-test is
 * single-process and would otherwise have every rank
 * popen() its own copy of ./test/run.  The runner's MPI
 * mode does not selectively skip native tests, so the
 * guard lives here.  OMPI_COMM_WORLD_SIZE is set by
 * Open MPI's mpirun in every spawned child.
 */
static bool under_mpirun(void)
{
	return getenv("OMPI_COMM_WORLD_SIZE") != NULL
		|| getenv("PMI_SIZE") != NULL;
}

int main(void)
{
	if (under_mpirun())
		return 0;

	t_pass_one();
	t_fail_one();
	t_mixed();
	t_abort_signal();
	t_failure_dump();
	t_filter();
	t_parallel_speedup();
	t_serial_no_speedup();
	t_no_color();
	t_force_color();
	t_no_paths();
	return 0;
}
