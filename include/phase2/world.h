#ifndef WORLD_H
#define WORLD_H

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"

enum log_level {
	LOG_TRACE,
	LOG_DEBUG,
	LOG_INFO,
	LOG_WARN,
	LOG_ERROR,
	LOG_FATAL
};

void log_log(int level, const char *fmt, ...);

#define log_trace(...) log_log(LOG_TRACE, __VA_ARGS__)
#define log_debug(...) log_log(LOG_DEBUG, __VA_ARGS__)
#define log_info(...) log_log(LOG_INFO, __VA_ARGS__)
#define log_warn(...) log_log(LOG_WARN, __VA_ARGS__)
#define log_error(...) log_log(LOG_ERROR, __VA_ARGS__)
#define log_fatal(...) log_log(LOG_FATAL, __VA_ARGS__)

enum world_stat {
	WORLD_UNDEF	= -1,
	WORLD_READY	=  0,
	WORLD_DONE	=  1,
	WORLD_ERR	=  2,
};

/*
 * Global world info stored as static data inside the module.
 * Users can keep a local copy of this when populated by world_info().
 */
struct world {
	enum world_stat stat;
	int size;
	int rank;
};

int world_init(int *argc, char ***argv);
int world_fin(void);

int world_info(struct world *wd);

#endif /* WORLD_H */
