/*
 * t-log -- unit tests for the phase2 logging subsystem.
 *
 * Built with -DDEBUG (the default for every test target),
 * so log_trace and log_debug compile to real emits.  Tests
 * capture stdout / stderr through a self-pipe, drain the
 * pipe, and assert against the captured bytes.
 *
 * The complementary release-build coverage lives in
 * test/t-log_release.c -- that binary is compiled with
 * -UDEBUG and asserts the trace/debug strip.
 */

/* pipe(), dup2(), read(), close() need POSIX. */
#define _POSIX_C_SOURCE 200809L

#define LOG_SUBSYS "t-log"

#include "c23_compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "log.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xa4f5c91d3e7b8062)
#define CAP_BYTES (16 * 1024)

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
	c->pipe_w = -1;

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

static void env_reset(void)
{
	unsetenv(PHASE2_LOG_ENVVAR);
	unsetenv(PHASE2_LOG_ALL_ENVVAR);
	log_fini();
}

static int digits(const char *s, int n)
{
	for (int i = 0; i < n; i++)
		if (s[i] < '0' || s[i] > '9')
			return 0;
	return 1;
}

/* -- tests ------------------------------------------------------------- */

static void t_default_level(void)
{
	env_reset();
	TEST_EQ(log_init(), 0);
	TEST_EQ(log_threshold, LOG_INFO);
}

static void t_parse_levels(void)
{
	static const struct {
		const char *env;
		int level;
	} cases[] = {
		{ "trace", LOG_TRACE },
		{ "Debug", LOG_DEBUG },
		{ "INFO", LOG_INFO },
		{ "w", LOG_WARN },
		{ "error", LOG_ERROR },
		{ "F", LOG_FATAL },
	};
	for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
		env_reset();
		setenv(PHASE2_LOG_ENVVAR, cases[i].env, 1);
		TEST_EQ(log_init(), 0);
		TEST_EQ(log_threshold, cases[i].level);
	}
}

static void t_parse_garbage(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "xyz", 1);
	TEST_EQ(log_init(), 0);
	/* Unparseable -> falls back to the default. */
	TEST_EQ(log_threshold, LOG_INFO);
}

static void t_threshold_filters(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "warn", 1);
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	log_info("hidden info");
	capture_end(&cap, buf, sizeof buf);
	TEST_ASSERT(strstr(buf, "hidden info") == NULL,
		"info should be suppressed when threshold=warn");

	TEST_EQ(capture_start(&cap, STDERR_FILENO), 0);
	log_warn("visible warn");
	capture_end(&cap, buf, sizeof buf);
	TEST_ASSERT(strstr(buf, "visible warn") != NULL,
		"warn must pass the warn threshold");
}

static void t_stream_split(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "trace", 1);
	log_init();

	char out[CAP_BYTES], err[CAP_BYTES];
	struct capture co, ce;
	TEST_EQ(capture_start(&co, STDOUT_FILENO), 0);
	TEST_EQ(capture_start(&ce, STDERR_FILENO), 0);
	log_trace("t-msg");
	log_debug("d-msg");
	log_info("i-msg");
	log_warn("w-msg");
	log_error("e-msg");
	log_fatal("f-msg");
	capture_end(&co, out, sizeof out);
	capture_end(&ce, err, sizeof err);

	TEST_ASSERT(strstr(out, "t-msg") != NULL, "trace on stdout");
	TEST_ASSERT(strstr(out, "d-msg") != NULL, "debug on stdout");
	TEST_ASSERT(strstr(out, "i-msg") != NULL, "info on stdout");
	TEST_ASSERT(strstr(out, "w-msg") == NULL, "warn NOT on stdout");
	TEST_ASSERT(strstr(out, "e-msg") == NULL, "error NOT on stdout");
	TEST_ASSERT(strstr(out, "f-msg") == NULL, "fatal NOT on stdout");

	TEST_ASSERT(strstr(err, "w-msg") != NULL, "warn on stderr");
	TEST_ASSERT(strstr(err, "e-msg") != NULL, "error on stderr");
	TEST_ASSERT(strstr(err, "f-msg") != NULL, "fatal on stderr");
	TEST_ASSERT(strstr(err, "i-msg") == NULL, "info NOT on stderr");
	TEST_ASSERT(strstr(err, "d-msg") == NULL, "debug NOT on stderr");
}

static void t_format_shape(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "info", 1);
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	log_info("payload-XYZ");
	capture_end(&cap, buf, sizeof buf);

	/* Expect: YYYY-MM-DD HH:MM:SS.mmm INFO  [t-log] payload-XYZ\n */
	TEST_ASSERT(digits(buf + 0, 4), "year digits");
	TEST_ASSERT(buf[4] == '-', "dash after year");
	TEST_ASSERT(digits(buf + 5, 2), "month digits");
	TEST_ASSERT(buf[7] == '-', "dash after month");
	TEST_ASSERT(digits(buf + 8, 2), "day digits");
	TEST_ASSERT(buf[10] == ' ', "space after date");
	TEST_ASSERT(digits(buf + 11, 2), "hour digits");
	TEST_ASSERT(buf[13] == ':', "colon");
	TEST_ASSERT(digits(buf + 14, 2), "minute digits");
	TEST_ASSERT(buf[16] == ':', "colon");
	TEST_ASSERT(digits(buf + 17, 2), "second digits");
	TEST_ASSERT(buf[19] == '.', "dot before ms");
	TEST_ASSERT(digits(buf + 20, 3), "millisecond digits");
	TEST_ASSERT(buf[23] == ' ', "space after ms");
	TEST_ASSERT(strstr(buf, " INFO  ") != NULL, "INFO level pad");
	TEST_ASSERT(strstr(buf, "[t-log]") != NULL, "subsystem tag");
	TEST_ASSERT(strstr(buf, "payload-XYZ\n") != NULL, "message + newline");
}

static void t_subsystem_explicit(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "info", 1);
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	/* Bypass the macro to install a custom subsystem. */
	log_emit(LOG_INFO, "myunit", __FILE__, __LINE__, "hello %s", "world");
	capture_end(&cap, buf, sizeof buf);

	TEST_ASSERT(strstr(buf, "[myunit]") != NULL,
		"explicit subsystem appears verbatim");
	TEST_ASSERT(strstr(buf, "hello world") != NULL, "printf passthrough");
}

static void t_subsystem_missing(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "info", 1);
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	/* LOG_SUBSYS is defined as "?" by the header when no
	 * file-level override is given.  Simulate by passing
	 * "?" directly. */
	log_emit(LOG_INFO, "?", __FILE__, __LINE__, "anon");
	capture_end(&cap, buf, sizeof buf);
	TEST_ASSERT(strstr(buf, "[?]") != NULL, "missing LOG_SUBSYS -> [?]");
}

static void t_printf_passthrough(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "info", 1);
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	log_info("d=%d s=%s z=%zu g=%g p=%%", -42, "abc",
		(size_t)1234567890u, 3.14);
	capture_end(&cap, buf, sizeof buf);
	TEST_ASSERT(strstr(buf, "d=-42 s=abc z=1234567890 g=3.14 p=%") != NULL,
		"printf format spec round-trip");
}

static void t_long_message(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "info", 1);
	log_init();

	char big[8192];
	memset(big, 'x', sizeof big - 1);
	big[sizeof big - 1] = '\0';

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	log_info("%s", big);
	int n = capture_end(&cap, buf, sizeof buf);
	TEST_ASSERT(n > 0, "long message produced output");
	TEST_ASSERT(buf[n - 1] == '\n', "output ends with newline despite trunc");
}

static void t_rank_prefix_off(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "info", 1);
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	log_info("no-prefix");
	capture_end(&cap, buf, sizeof buf);
	/* Rank 0, ALL off: no [r=...] prefix; line starts
	 * with a digit (the year). */
	TEST_ASSERT(buf[0] >= '0' && buf[0] <= '9',
		"rank-0 line starts with year, no rank prefix");
}

static void t_init_idempotent(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "info", 1);
	TEST_EQ(log_init(), 0);
	/* Calling again is a no-op; threshold unchanged. */
	const int before = log_threshold;
	setenv(PHASE2_LOG_ENVVAR, "warn", 1);
	TEST_EQ(log_init(), 0);
	TEST_EQ(log_threshold, before);

	/* After fini, init reloads env. */
	log_fini();
	TEST_EQ(log_init(), 0);
	TEST_EQ(log_threshold, LOG_WARN);
}

static void t_filtered_perf(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "warn", 1);
	log_init();

	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC, &t0);
	for (int i = 0; i < 1000000; i++)
		log_trace("filtered: %d", i);
	clock_gettime(CLOCK_MONOTONIC, &t1);

	const double dt = (t1.tv_sec - t0.tv_sec) +
		(t1.tv_nsec - t0.tv_nsec) * 1e-9;
	TEST_ASSERT(dt < 0.25,
		"1e6 filtered log_trace calls in %.3f s (budget <0.25 s)",
		dt);
}

static void t_trace_debug_emit_in_debug_build(void)
{
	env_reset();
	setenv(PHASE2_LOG_ENVVAR, "trace", 1);
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	TEST_EQ(capture_start(&cap, STDOUT_FILENO), 0);
	log_trace("trace-emit");
	log_debug("debug-emit");
	capture_end(&cap, buf, sizeof buf);

#ifdef DEBUG
	TEST_ASSERT(strstr(buf, "trace-emit") != NULL,
		"DEBUG build emits trace");
	TEST_ASSERT(strstr(buf, "debug-emit") != NULL,
		"DEBUG build emits debug");
#else
	TEST_ASSERT(strstr(buf, "trace-emit") == NULL,
		"non-DEBUG build strips trace");
	TEST_ASSERT(strstr(buf, "debug-emit") == NULL,
		"non-DEBUG build strips debug");
#endif
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_default_level();
	t_parse_levels();
	t_parse_garbage();
	t_threshold_filters();
	t_stream_split();
	t_format_shape();
	t_subsystem_explicit();
	t_subsystem_missing();
	t_printf_passthrough();
	t_long_message();
	t_rank_prefix_off();
	t_init_idempotent();
	t_filtered_perf();
	t_trace_debug_emit_in_debug_build();

	world_free();
	return 0;
}
