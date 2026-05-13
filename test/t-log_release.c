/*
 * t-log_release -- verifies that log_trace and log_debug
 * vanish from a build that does NOT define DEBUG.
 *
 * The Makefile builds this single test target with
 * -UDEBUG cancelling the project-wide -DDEBUG that
 * applies to every other test.  Inside this translation
 * unit, log_trace and log_debug must expand to ((void)0)
 * and therefore produce no bytes on either sink even with
 * log_threshold = LOG_TRACE.
 */

#define _POSIX_C_SOURCE 200809L

#define LOG_SUBSYS "t-log_release"

#include "c23_compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "log.h"
#include "phase2/world.h"

#include "test.h"

#define WD_SEED UINT64_C(0xa4f5c91d3e7b8062)
#define CAP_BYTES (4096)

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
	if (c->saved_fd < 0 || dup2(c->pipe_w, target_fd) < 0) {
		close(c->pipe_r);
		close(c->pipe_w);
		return -1;
	}
	return 0;
}

static int capture_end(struct capture *c, char *buf, size_t cap)
{
	if (c->target_fd == STDOUT_FILENO)
		fflush(stdout);
	else if (c->target_fd == STDERR_FILENO)
		fflush(stderr);
	dup2(c->saved_fd, c->target_fd);
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

static void t_trace_debug_stripped(void)
{
#ifdef DEBUG
	TEST_FAIL("t-log_release MUST be built with -UDEBUG");
#endif

	unsetenv(PHASE2_LOG_ENVVAR);
	unsetenv(PHASE2_LOG_ALL_ENVVAR);
	log_fini();
	log_init();
	log_threshold = LOG_TRACE; /* lowest gate -- everything in */

	char out[CAP_BYTES], err[CAP_BYTES];
	struct capture co, ce;
	capture_start(&co, STDOUT_FILENO);
	capture_start(&ce, STDERR_FILENO);
	log_trace("MUST-NOT-APPEAR-trace");
	log_debug("MUST-NOT-APPEAR-debug");
	capture_end(&co, out, sizeof out);
	capture_end(&ce, err, sizeof err);

	TEST_ASSERT(out[0] == '\0',
		"log_trace/log_debug must produce zero bytes on stdout"
		" in release builds (got: %s)", out);
	TEST_ASSERT(err[0] == '\0',
		"log_trace/log_debug must produce zero bytes on stderr"
		" in release builds (got: %s)", err);
}

static void t_info_still_works(void)
{
	unsetenv(PHASE2_LOG_ENVVAR);
	log_fini();
	log_init();

	char buf[CAP_BYTES];
	struct capture cap;
	capture_start(&cap, STDOUT_FILENO);
	log_info("info-survives-release");
	capture_end(&cap, buf, sizeof buf);
	TEST_ASSERT(strstr(buf, "info-survives-release") != NULL,
		"log_info must still emit in release builds");
}

int main(void)
{
	world_init(nullptr, nullptr, WD_SEED);

	t_trace_debug_stripped();
	t_info_still_works();

	world_free();
	return 0;
}
