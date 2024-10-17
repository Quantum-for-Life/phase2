#ifndef LOG_H
#define LOG_H

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"

#ifdef __cplusplus
extern "C" {
#endif

enum log_level {
	LOG_TRACE,
	LOG_DEBUG,
	LOG_INFO,
	LOG_WARN,
	LOG_ERROR,
	LOG_FATAL
};

/* Call world_init() before using this. */
void log_log(int level, const char *fmt, ...);

#define log_trace(...) log_log(LOG_TRACE, __VA_ARGS__)
#define log_debug(...) log_log(LOG_DEBUG, __VA_ARGS__)
#define log_info(...) log_log(LOG_INFO, __VA_ARGS__)
#define log_warn(...) log_log(LOG_WARN, __VA_ARGS__)
#define log_error(...) log_log(LOG_ERROR, __VA_ARGS__)
#define log_fatal(...) log_log(LOG_FATAL, __VA_ARGS__)

int log_init(void);

#ifdef __cplusplus
}
#endif

#endif /* LOG_H */
