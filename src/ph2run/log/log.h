/**
 * Copyright (c) 2020 rxi
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the MIT license. See `log.c` for details.
 */

#ifndef LOG_H
#define LOG_H

#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <time.h>

#define LOG_VERSION "0.1.0"

struct log_event {
	va_list ap;
	const char *fmt;
	const char *file;
	struct tm *time;
	void *data;
	int line;
	int level;
};

typedef void (*log_logfn)(struct log_event *ev);

typedef void (*log_lockfn)(bool lock, void *data);

enum log_level {
	LOG_TRACE,
	LOG_DEBUG,
	LOG_INFO,
	LOG_WARN,
	LOG_ERROR,
	LOG_FATAL
};

#define log_trace(...) log_log(LOG_TRACE, __VA_ARGS__)
#define log_debug(...) log_log(LOG_DEBUG, __VA_ARGS__)
#define log_info(...) log_log(LOG_INFO, __VA_ARGS__)
#define log_warn(...) log_log(LOG_WARN, __VA_ARGS__)
#define log_error(...) log_log(LOG_ERROR, __VA_ARGS__)
#define log_fatal(...) log_log(LOG_FATAL, __VA_ARGS__)

const char *log_level_string(int level);

int log_level_from_lowercase(enum log_level *level, const char *name);

void log_set_lock(log_lockfn fn, void *data);

void log_set_level(int level);

int log_add_callback(log_logfn fn, void *data, int level);

void log_vlog(const int level, const char *fmt, va_list ap);

void log_log(int level, const char *fmt, ...);

#endif
