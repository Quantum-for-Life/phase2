/*
 * test/run -- parallel test harness for phase2.
 *
 * The runner takes a list of test paths on the command
 * line, spawns up to N children at once (N defaults to
 * the number of online cores), captures each test's
 * stdout and stderr to per-test tempfiles, and reports
 * results in a cargo-test-style stream:
 *
 *     running 31 tests
 *     test t-bitstring_index ............ ok (0.04s)
 *     test t-data_mpi ................... FAILED (1.23s)
 *     ...
 *
 *     failures:
 *     ---- t-data_mpi stdout ----
 *     <captured>
 *     ---- t-data_mpi stderr ----
 *     <captured>
 *
 *     test result: FAILED. 30 passed; 1 failed; \
 *         finished in 6.32s
 *
 * Architectural notes
 * -------------------
 *
 *  - Tests are *not* discovered here.  The caller (the
 *    Makefile, or a developer invoking the binary
 *    directly) passes the list of test paths after `--`.
 *    The runner stays decoupled from the TESTS /
 *    TESTS_SLOW manifest.
 *
 *  - Each child is a plain fork+execvp.  stdout and
 *    stderr are redirected to mkstemp(3) tempfiles in a
 *    per-run mkdtemp(3) directory under /tmp.  Parent
 *    waits with waitpid(-1, &status, 0) so completions
 *    are processed in arrival order, then dispatches
 *    the next pending test into the freed slot.
 *
 *  - Python tests (suffix `.py`) are recognised and
 *    executed via `python3 <path>` instead of running
 *    the path as a native binary.  Under MPI mode
 *    (--mpiranks=N) Python tests are skipped: they are
 *    rank-0 cross-validators and don't make sense under
 *    mpirun.
 *
 *  - Colour is on iff stdout is a TTY (or FORCE_COLOR is
 *    set).  --no-color forces plain output.  The plain
 *    format is byte-identical apart from the ANSI
 *    escapes, so logs are diffable.
 *
 *  - SIGINT propagates naturally: the children share
 *    the parent's process group, so Ctrl-C at the
 *    terminal kills all live children and the parent
 *    sees them as terminated.  No custom handler.
 *
 * No phase2 / MPI / HDF5 dependencies; libc + POSIX
 * only.  Build line:
 *
 *     cc -std=c11 -Wall -Wextra -O3 test/run.c -o test/run
 */

#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <fnmatch.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>


/* -- diagnostics ----------------------------------------------------- */

/*
 * Print a fatal error to stderr (with errno text if
 * non-zero) and exit non-zero.  Used for unrecoverable
 * setup failures (mkstemp, fork, ...).  Test failures
 * proper go through the regular pass/fail path; this is
 * for the harness itself dying.
 */
__attribute__((noreturn, format(printf, 1, 2)))
static void die(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fputs("test/run: ", stderr);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	if (errno)
		fprintf(stderr, ": %s", strerror(errno));
	fputc('\n', stderr);
	exit(2);
}


/* -- options --------------------------------------------------------- */

struct opts {
	int jobs;		/* parallel slots; 1..N */
	int mpiranks;		/* 0 = no MPI wrap */
	const char *filter;	/* fnmatch glob; NULL = no filter */
	bool verbose;		/* don't capture stdout; forces jobs=1 */
	bool color;		/* emit ANSI escapes? */
};

/*
 * sysconf-derived default for --jobs.  Falls back to 1
 * if the platform refuses to answer (some sandboxes).
 */
static int default_jobs(void)
{
	const long n = sysconf(_SC_NPROCESSORS_ONLN);
	return (n > 0) ? (int)n : 1;
}

/*
 * Colour-on-TTY policy: the user's --no-color flag wins
 * absolutely.  Otherwise honour FORCE_COLOR (set in CI
 * configs that wrap colour-aware tools), and finally
 * fall back to isatty(stdout).
 */
static bool default_color(void)
{
	if (getenv("FORCE_COLOR"))
		return true;
	return isatty(STDOUT_FILENO) != 0;
}

static void usage(FILE *out)
{
	fputs(
		"usage: test/run [OPTIONS] -- TEST...\n"
		"\n"
		"  --jobs=N         parallel job count "
		"(default: number of cores)\n"
		"  --mpiranks=N     wrap each test with "
		"`mpirun -n N` (skips .py tests)\n"
		"  --filter=GLOB    fnmatch against the "
		"effective test name (no\n"
		"                   `t-` prefix, no `.py` "
		"suffix)\n"
		"  --verbose, -v    pass each test's stdout "
		"through; forces --jobs=1\n"
		"  --no-color       force plain output "
		"(default: detect TTY)\n"
		"  --help, -h       this help\n"
		"\n"
		"  TEST...          list of test paths (native "
		"binaries and/or .py)\n",
		out);
}

static void parse_opts(int argc, char **argv, struct opts *o, int *first_pos)
{
	*o = (struct opts){
		.jobs = default_jobs(),
		.mpiranks = 0,
		.filter = NULL,
		.verbose = false,
		.color = default_color(),
	};

	static const struct option long_opts[] = {
		{ "jobs",     required_argument, NULL, 'j' },
		{ "mpiranks", required_argument, NULL, 'm' },
		{ "filter",   required_argument, NULL, 'f' },
		{ "verbose",  no_argument,       NULL, 'v' },
		{ "no-color", no_argument,       NULL, 'C' },
		{ "help",     no_argument,       NULL, 'h' },
		{ 0 },
	};
	int c;
	while ((c = getopt_long(argc, argv, "vh", long_opts, NULL)) != -1) {
		switch (c) {
		case 'j': o->jobs = atoi(optarg); break;
		case 'm': o->mpiranks = atoi(optarg); break;
		case 'f': o->filter = optarg; break;
		case 'v': o->verbose = true; break;
		case 'C': o->color = false; break;
		case 'h': usage(stdout); exit(0);
		default:  usage(stderr); exit(2);
		}
	}
	if (o->verbose)
		o->jobs = 1;	/* keep per-test output legible */
	if (o->jobs < 1)
		o->jobs = 1;
	*first_pos = optind;
}


/* -- test_result ----------------------------------------------------- */

/*
 * Per-test state carried through the run.  Populated at
 * setup time (path/name/subsys) and filled in as the
 * test runs (pid -> exit_code, elapsed, captured-output
 * paths).
 */
struct test_result {
	const char *path;	/* exactly as supplied on argv */
	char *name;		/* basename(path) minus t- prefix
				 * and .py suffix; owned */
	bool is_python;		/* path ends in `.py` */
	bool dispatched;	/* spawned? */
	bool reaped;		/* completed and recorded? */
	bool passed;		/* exit_code == 0 */
	pid_t pid;
	int exit_code;		/* WEXITSTATUS, or 128+sig on signal */
	double elapsed;		/* wall-clock seconds */
	struct timespec t0;	/* spawn time */
	char *out_path;		/* captured stdout (owned) */
	char *err_path;		/* captured stderr (owned) */
	int out_fd, err_fd;	/* open fds during spawn; closed in parent */
};

/*
 * Compute a malloc'd "effective name" -- basename minus
 * the conventional `t-` prefix and the `.py` suffix if
 * present.  Examples:
 *
 *     test/t-data_mpi             -> "data_mpi"
 *     ./test/t-ref-coeff_matrix.py -> "ref-coeff_matrix"
 *
 * Caller owns the result.
 */
static char *make_name(const char *path)
{
	const char *slash = strrchr(path, '/');
	const char *base = slash ? slash + 1 : path;
	if (strncmp(base, "t-", 2) == 0)
		base += 2;
	const size_t len = strlen(base);
	const size_t cut = (len >= 3 && strcmp(base + len - 3, ".py") == 0)
		? len - 3 : len;
	char *out = malloc(cut + 1);
	if (!out)
		die("malloc");
	memcpy(out, base, cut);
	out[cut] = '\0';
	return out;
}

/* -- filtering ------------------------------------------------------- */

/*
 * Drop tests whose effective name doesn't match `glob`
 * (fnmatch with default flags -- POSIX glob syntax,
 * `*` does not cross `/` but our names are flat).
 *
 * `*n_tests` shrinks in place; freed entries are *not*
 * reclaimed (the run is short-lived).
 */
static void filter_tests(struct test_result *tests, size_t *n_tests,
	const char *glob)
{
	if (!glob)
		return;
	size_t w = 0;
	for (size_t r = 0; r < *n_tests; r++) {
		if (fnmatch(glob, tests[r].name, 0) == 0)
			tests[w++] = tests[r];
	}
	*n_tests = w;
}


/* -- child spawn ----------------------------------------------------- */

/*
 * Prepare per-test stdout/stderr tempfiles, fork, and
 * execvp into the test (or python3, or mpirun).  Records
 * the pid + spawn time on `r`; bubbles fatal errno up
 * via die().
 *
 * The tempfile fds are opened in the parent before fork
 * so both processes see the same inode.  The child
 * dup2's them onto fds 1 and 2, then closes the
 * originals.  The parent closes its copies immediately
 * after fork -- the file remains open via the child's
 * stdio handles, and we'll reopen it for reading once
 * the child exits.
 */
static void spawn_one(struct test_result *r, const struct opts *o,
	const char *tmpdir)
{
	/* Tempfile names: <tmpdir>/<name>.out and .err.
	 * The name is already filename-safe (alnum + `_` +
	 * `-`); no escaping needed. */
	const size_t pathsz = strlen(tmpdir) + strlen(r->name) + 8;
	r->out_path = malloc(pathsz);
	r->err_path = malloc(pathsz);
	if (!r->out_path || !r->err_path)
		die("malloc");
	snprintf(r->out_path, pathsz, "%s/%s.out", tmpdir, r->name);
	snprintf(r->err_path, pathsz, "%s/%s.err", tmpdir, r->name);
	r->out_fd = open(r->out_path,
		O_CREAT | O_TRUNC | O_RDWR, 0600);
	r->err_fd = open(r->err_path,
		O_CREAT | O_TRUNC | O_RDWR, 0600);
	if (r->out_fd < 0 || r->err_fd < 0)
		die("open(%s/%s)", tmpdir, r->name);

	clock_gettime(CLOCK_MONOTONIC, &r->t0);

	const pid_t pid = fork();
	if (pid < 0)
		die("fork");

	if (pid == 0) {
		/* child */
		if (o->verbose) {
			/* Stream through to the parent's TTY for
			 * live observation.  No capture. */
			close(r->out_fd);
			close(r->err_fd);
		} else {
			if (dup2(r->out_fd, STDOUT_FILENO) < 0
				|| dup2(r->err_fd, STDERR_FILENO) < 0)
				_exit(126);
			close(r->out_fd);
			close(r->err_fd);
		}

		/* argv assembly: at most six slots (mpirun -n
		 * NN python3 path NULL).  Static buffer is
		 * fine -- it lives until exec consumes it. */
		const char *argv[8] = { 0 };
		char ranks_buf[16];
		int ai = 0;
		if (o->mpiranks > 0 && !r->is_python) {
			snprintf(ranks_buf, sizeof ranks_buf, "%d",
				o->mpiranks);
			argv[ai++] = "mpirun";
			argv[ai++] = "-n";
			argv[ai++] = ranks_buf;
		}
		if (r->is_python)
			argv[ai++] = "python3";
		argv[ai++] = r->path;
		argv[ai] = NULL;

		execvp(argv[0], (char *const *)argv);
		_exit(127);
	}

	/* parent */
	close(r->out_fd);
	close(r->err_fd);
	r->out_fd = r->err_fd = -1;
	r->pid = pid;
	r->dispatched = true;
}


/* -- parent reap ----------------------------------------------------- */

/*
 * Wait for any one child to complete.  Find its
 * test_result by pid, record exit code + elapsed, mark
 * reaped.  Returns the index of the completed test, or
 * -1 if no child was outstanding.
 */
static ssize_t reap_one(struct test_result *tests, size_t n_tests)
{
	int status;
	const pid_t pid = waitpid(-1, &status, 0);
	if (pid < 0) {
		if (errno == ECHILD)
			return -1;
		die("waitpid");
	}
	struct timespec t1;
	clock_gettime(CLOCK_MONOTONIC, &t1);

	for (size_t i = 0; i < n_tests; i++) {
		if (!tests[i].dispatched || tests[i].reaped)
			continue;
		if (tests[i].pid != pid)
			continue;
		tests[i].elapsed =
			(double)(t1.tv_sec - tests[i].t0.tv_sec)
			+ (double)(t1.tv_nsec - tests[i].t0.tv_nsec) * 1e-9;
		if (WIFEXITED(status))
			tests[i].exit_code = WEXITSTATUS(status);
		else if (WIFSIGNALED(status))
			tests[i].exit_code = 128 + WTERMSIG(status);
		else
			tests[i].exit_code = -1;
		tests[i].passed = (tests[i].exit_code == 0);
		tests[i].reaped = true;
		return (ssize_t)i;
	}
	/* Should not happen: kernel handed us a pid we
	 * never spawned.  Surface it loudly. */
	die("waitpid returned unknown pid %d", (int)pid);
}


/* -- output ---------------------------------------------------------- */

/* ANSI sequences; emitted only when opts->color is on. */
#define ANSI_GREEN "\x1b[32m"
#define ANSI_RED   "\x1b[31m"
#define ANSI_DIM   "\x1b[2m"
#define ANSI_RESET "\x1b[0m"

/*
 * Width of the longest effective name across `tests`,
 * used to align the dot padding so the status column
 * lines up visually.
 */
static size_t name_column_width(const struct test_result *tests, size_t n)
{
	size_t w = 0;
	for (size_t i = 0; i < n; i++) {
		const size_t len = strlen(tests[i].name) + 2; /* `t-` */
		if (len > w)
			w = len;
	}
	return w;
}

/*
 * Print one test's status line.  The format is
 *
 *     test <name> <dots> <verdict> (<elapsed>s)
 *
 * with the dot run padding `name` out to `name_col`
 * columns + 4 minimum dots.
 */
static void print_status(const struct test_result *r, const struct opts *o,
	size_t name_col)
{
	const char *g_on = o->color ? ANSI_GREEN : "";
	const char *r_on = o->color ? ANSI_RED   : "";
	const char *d_on = o->color ? ANSI_DIM   : "";
	const char *off  = o->color ? ANSI_RESET : "";

	const size_t name_len = strlen(r->name) + 2; /* "t-" */
	const size_t dots = (name_col >= name_len)
		? (name_col - name_len) + 4 : 4;

	printf("test t-%s ", r->name);
	fputs(d_on, stdout);
	for (size_t i = 0; i < dots; i++)
		fputc('.', stdout);
	fputs(off, stdout);
	if (r->passed)
		printf(" %sok%s", g_on, off);
	else
		printf(" %sFAILED%s", r_on, off);
	printf(" %s(%.3fs)%s\n", d_on, r->elapsed, off);
	fflush(stdout);
}

/*
 * Dump a file by path with a labelled header.  Used to
 * surface captured stdout / stderr of a failed test.
 * Silently skips empty files.
 */
static void dump_file(const char *label, const char *path)
{
	struct stat st;
	if (stat(path, &st) < 0 || st.st_size == 0)
		return;
	printf("---- %s ----\n", label);
	FILE *f = fopen(path, "r");
	if (!f) {
		printf("(could not open %s: %s)\n", path, strerror(errno));
		return;
	}
	char buf[4096];
	size_t r, total = 0;
	int last_ch = '\n';
	while ((r = fread(buf, 1, sizeof buf, f)) > 0) {
		fwrite(buf, 1, r, stdout);
		total += r;
		last_ch = (unsigned char)buf[r - 1];
	}
	fclose(f);
	if (total > 0 && last_ch != '\n')
		fputc('\n', stdout);
}

/*
 * After the run completes, print captured output of
 * each failed test and the final pass/fail summary
 * line.
 */
static int print_summary(const struct test_result *tests, size_t n_tests,
	double total, const struct opts *o)
{
	const char *g_on = o->color ? ANSI_GREEN : "";
	const char *r_on = o->color ? ANSI_RED   : "";
	const char *d_on = o->color ? ANSI_DIM   : "";
	const char *off  = o->color ? ANSI_RESET : "";

	size_t passed = 0, failed = 0;
	for (size_t i = 0; i < n_tests; i++) {
		if (tests[i].passed)
			passed++;
		else
			failed++;
	}

	if (failed > 0) {
		printf("\nfailures:\n\n");
		for (size_t i = 0; i < n_tests; i++) {
			if (tests[i].passed)
				continue;
			char buf[128];
			snprintf(buf, sizeof buf, "t-%s stdout",
				tests[i].name);
			dump_file(buf, tests[i].out_path);
			snprintf(buf, sizeof buf, "t-%s stderr",
				tests[i].name);
			dump_file(buf, tests[i].err_path);
			printf("(exit %d)\n\n", tests[i].exit_code);
		}
		printf("failures:\n");
		for (size_t i = 0; i < n_tests; i++) {
			if (!tests[i].passed)
				printf("    t-%s\n", tests[i].name);
		}
	}

	const char *verdict_col = (failed == 0) ? g_on : r_on;
	const char *verdict = (failed == 0) ? "ok" : "FAILED";
	printf("\ntest result: %s%s%s. %zu passed; %zu failed;"
	       " %sfinished in %.2fs%s\n",
		verdict_col, verdict, off, passed, failed,
		d_on, total, off);

	return (failed == 0) ? 0 : 1;
}


/* -- tempdir cleanup ------------------------------------------------- */

/*
 * Remove the per-run tempdir and any files inside it.
 * Best-effort: any errors are reported but don't escalate
 * to a non-zero exit, since the test results matter more
 * than tempdir hygiene.  The tempfiles named in `tests`
 * are the only contents we ever created, so we don't need
 * a directory walk.
 */
static void cleanup_tmpdir(const char *tmpdir,
	const struct test_result *tests, size_t n_tests)
{
	for (size_t i = 0; i < n_tests; i++) {
		if (tests[i].out_path)
			unlink(tests[i].out_path);
		if (tests[i].err_path)
			unlink(tests[i].err_path);
	}
	if (rmdir(tmpdir) < 0)
		perror("test/run: rmdir tmpdir");
}


/* -- main loop ------------------------------------------------------- */

int main(int argc, char **argv)
{
	struct opts o;
	int first_pos;
	parse_opts(argc, argv, &o, &first_pos);

	/*
	 * HDF5 1.10+ takes an advisory flock(2) on every
	 * file it opens; when N parallel children each
	 * H5Fopen the same committed test fixture under
	 * test/data/, the second-and-later opens
	 * intermittently fail.  All our tests open
	 * fixtures read-only, so disabling the locking
	 * is safe and stops the parallel runs from
	 * flaking.  Setting the env var here means every
	 * fork+exec'd child inherits it.
	 */
	setenv("HDF5_USE_FILE_LOCKING", "FALSE", 0);

	const int n_paths = argc - first_pos;
	if (n_paths <= 0) {
		fputs("test/run: no test paths supplied (after `--`)\n",
			stderr);
		usage(stderr);
		exit(2);
	}

	/* Build the test_result array.  In MPI mode, Python
	 * tests get filtered out at this stage so the user-
	 * supplied list can be common between MPI and
	 * non-MPI invocations. */
	struct test_result *tests = calloc(n_paths, sizeof *tests);
	if (!tests)
		die("calloc");
	size_t n_tests = 0;
	for (int i = 0; i < n_paths; i++) {
		const char *p = argv[first_pos + i];
		const size_t plen = strlen(p);
		const bool py = (plen > 3 && strcmp(p + plen - 3, ".py") == 0);
		if (py && o.mpiranks > 0)
			continue;
		tests[n_tests] = (struct test_result){
			.path = p,
			.is_python = py,
		};
		tests[n_tests].name = make_name(p);
		n_tests++;
	}
	filter_tests(tests, &n_tests, o.filter);
	if (n_tests == 0) {
		fputs("test/run: no tests after filter\n", stderr);
		exit(0);
	}

	/* Per-run tempdir for captured output. */
	char tmpdir[] = "/tmp/phase2-test-run-XXXXXX";
	if (!mkdtemp(tmpdir))
		die("mkdtemp");

	const size_t name_col = name_column_width(tests, n_tests);
	printf("\nrunning %zu test%s\n\n", n_tests, n_tests == 1 ? "" : "s");

	struct timespec wall0;
	clock_gettime(CLOCK_MONOTONIC, &wall0);

	/* Single bounded loop: dispatch up to `jobs` at a
	 * time, reap as they finish, repeat until every
	 * test has both dispatched and reaped. */
	size_t next_pending = 0;
	size_t running = 0;
	while (next_pending < n_tests || running > 0) {
		while (running < (size_t)o.jobs && next_pending < n_tests) {
			spawn_one(&tests[next_pending], &o, tmpdir);
			next_pending++;
			running++;
		}
		if (running > 0) {
			const ssize_t done = reap_one(tests, n_tests);
			if (done < 0)
				break;
			print_status(&tests[done], &o, name_col);
			running--;
		}
	}

	struct timespec wall1;
	clock_gettime(CLOCK_MONOTONIC, &wall1);
	const double total =
		(double)(wall1.tv_sec - wall0.tv_sec)
		+ (double)(wall1.tv_nsec - wall0.tv_nsec) * 1e-9;

	const int rc = print_summary(tests, n_tests, total, &o);

	cleanup_tmpdir(tmpdir, tests, n_tests);
	for (size_t i = 0; i < n_tests; i++) {
		free(tests[i].name);
		free(tests[i].out_path);
		free(tests[i].err_path);
	}
	free(tests);

	return rc;
}
