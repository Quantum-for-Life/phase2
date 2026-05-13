#ifndef LOG_H
#define LOG_H

/*
 * phase2 logging.
 *
 * Format (left to right):
 *
 *   [r=R] YYYY-MM-DD HH:MM:SS.mmm LEVEL [subsys] message
 *
 * The rank prefix [r=R] appears only when PHASE2_LOG_ALL=1
 * AND the calling rank is non-zero.  Rank 0 stays unprefixed
 * in every mode.
 *
 * Sinks: TRACE/DEBUG/INFO go to stdout; WARN/ERROR/FATAL go
 * to stderr.  Both are fflushed after every emit so SLURM
 * `tail -f` shows output as it is produced.
 *
 * Each translation unit that emits log_* should set its own
 * subsystem tag before including this header:
 *
 *   #define LOG_SUBSYS "qreg"
 *   #include "log.h"
 *
 * Files without LOG_SUBSYS fall back to "?".  Reviewers
 * should reject patches that leave that fallback in place.
 *
 * Build-time gating: log_trace and log_debug expand to a
 * no-op when DEBUG is undefined (the default for `make
 * build`).  Release binaries pay zero cost for trace/debug
 * call sites — no load, no branch, no format string.
 * Tests and `make debug` builds keep them.
 */

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"
#define PHASE2_LOG_ALL_ENVVAR "PHASE2_LOG_ALL"

#ifndef LOG_SUBSYS
#define LOG_SUBSYS "?"
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum log_level {
	LOG_TRACE = 0,
	LOG_DEBUG = 1,
	LOG_INFO = 2,
	LOG_WARN = 3,
	LOG_ERROR = 4,
	LOG_FATAL = 5,
};

/* Runtime threshold.  Set by log_init() from PHASE2_LOG.
 * The log_at() macro reads it directly to short-circuit
 * filtered-out calls before any function is invoked. */
extern int log_threshold;

void log_emit(int level, const char *subsys, const char *file, int line,
	const char *fmt, ...) __attribute__((format(printf, 5, 6)));

#define log_at(lvl, ...)                                                       \
	do {                                                                   \
		if (__builtin_expect((lvl) >= log_threshold, 0))               \
			log_emit((lvl), LOG_SUBSYS, __FILE__, __LINE__,        \
				__VA_ARGS__);                                  \
	} while (0)

#ifdef DEBUG
#define log_trace(...) log_at(LOG_TRACE, __VA_ARGS__)
#define log_debug(...) log_at(LOG_DEBUG, __VA_ARGS__)
#else
#define log_trace(...) ((void)0)
#define log_debug(...) ((void)0)
#endif
#define log_info(...) log_at(LOG_INFO, __VA_ARGS__)
#define log_warn(...) log_at(LOG_WARN, __VA_ARGS__)
#define log_error(...) log_at(LOG_ERROR, __VA_ARGS__)
#define log_fatal(...) log_at(LOG_FATAL, __VA_ARGS__)

/* Read PHASE2_LOG and PHASE2_LOG_ALL; set stdout to
 * line-buffered, stderr to unbuffered.  Idempotent.
 * Returns 0 on success, -1 on error. */
int log_init(void);

/* Flush both sinks.  Safe to call without log_init(). */
void log_fini(void);

#ifdef __cplusplus
}
#endif

#endif /* LOG_H */
