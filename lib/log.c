/*
 * phase2 logging — implementation.
 *
 * See include/log.h for the public surface and the format
 * contract.  No dependencies beyond libc + POSIX time
 * (clock_gettime, localtime_r).  MPI is consulted only via
 * world_info() — a cheap struct copy.
 */

/* clock_gettime + localtime_r need POSIX.1-2008. */
#define _POSIX_C_SOURCE 200809L

#define LOG_SUBSYS "log"

#include "c23_compat.h"
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "log.h"
#include "phase2/world.h"

int log_threshold = LOG_INFO;

#define LINE_BUF_BYTES (4096)

static const char *const LEVEL_NAMES[] = {
	"TRACE",
	"DEBUG",
	"INFO",
	"WARN",
	"ERROR",
	"FATAL",
};

static struct log_state {
	int initialized;
	int all_ranks;
} state = { 0, 0 };

static int parse_level(const char *name, int *out)
{
	if (!name || !*name)
		return -1;
	switch (name[0]) {
	case 't':
	case 'T':
		*out = LOG_TRACE;
		return 0;
	case 'd':
	case 'D':
		*out = LOG_DEBUG;
		return 0;
	case 'i':
	case 'I':
		*out = LOG_INFO;
		return 0;
	case 'w':
	case 'W':
		*out = LOG_WARN;
		return 0;
	case 'e':
	case 'E':
		*out = LOG_ERROR;
		return 0;
	case 'f':
	case 'F':
		*out = LOG_FATAL;
		return 0;
	}
	return -1;
}

void log_emit(const int level, const char *subsys, const char *file,
	const int line, const char *fmt, ...)
{
	(void)file;
	(void)line;

	struct world_info wd;
	const int wstat = world_info(&wd);
	const int rank = (wstat >= 0) ? wd.rank : 0;

	/* Rank gating: non-zero ranks emit only when
	 * PHASE2_LOG_ALL is set.  Before log_init() the
	 * all_ranks flag is 0; rank-0 still emits. */
	if (rank > 0 && !state.all_ranks)
		return;

	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	struct tm tm;
	localtime_r(&ts.tv_sec, &tm);

	char hdr[128];
	int hlen;
	if (state.all_ranks && rank > 0)
		hlen = snprintf(hdr, sizeof hdr,
			"[r=%d] %04d-%02d-%02d %02d:%02d:%02d.%03ld "
			"%-5s [%s] ",
			rank, 1900 + tm.tm_year, 1 + tm.tm_mon, tm.tm_mday,
			tm.tm_hour, tm.tm_min, tm.tm_sec,
			ts.tv_nsec / 1000000, LEVEL_NAMES[level], subsys);
	else
		hlen = snprintf(hdr, sizeof hdr,
			"%04d-%02d-%02d %02d:%02d:%02d.%03ld "
			"%-5s [%s] ",
			1900 + tm.tm_year, 1 + tm.tm_mon, tm.tm_mday,
			tm.tm_hour, tm.tm_min, tm.tm_sec,
			ts.tv_nsec / 1000000, LEVEL_NAMES[level], subsys);
	if (hlen < 0)
		return;
	if ((size_t)hlen >= sizeof hdr)
		hlen = (int)sizeof hdr - 1;

	char body[LINE_BUF_BYTES];
	va_list ap;
	va_start(ap, fmt);
	int blen = vsnprintf(body, sizeof body, fmt, ap);
	va_end(ap);
	if (blen < 0)
		return;
	if ((size_t)blen >= sizeof body)
		blen = (int)sizeof body - 1;

	FILE *sink = (level >= LOG_WARN) ? stderr : stdout;
	fwrite(hdr, 1, (size_t)hlen, sink);
	fwrite(body, 1, (size_t)blen, sink);
	fputc('\n', sink);
	fflush(sink);
}

int log_init(void)
{
	if (state.initialized)
		return 0;

	int lvl = LOG_INFO;
	const char *env_lvl = getenv(PHASE2_LOG_ENVVAR);
	if (env_lvl) {
		int parsed;
		if (parse_level(env_lvl, &parsed) == 0)
			lvl = parsed;
	}
	log_threshold = lvl;

	const char *env_all = getenv(PHASE2_LOG_ALL_ENVVAR);
	state.all_ranks = (env_all && *env_all) ? 1 : 0;

	setvbuf(stdout, NULL, _IOLBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);

	state.initialized = 1;
	return 0;
}

void log_fini(void)
{
	fflush(stdout);
	fflush(stderr);
	state.initialized = 0;
}
