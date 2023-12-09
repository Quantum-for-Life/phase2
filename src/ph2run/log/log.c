/*
 * Copyright (c) 2020 rxi
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifdef DISTRIBUTED
#include <string.h>
#include "mpi.h"
#endif

#include "log.h"

#define MAX_CALLBACKS (32)

struct callback {
	log_logfn fn;
	void *data;
	int level;
};

static struct {
	void *data;
	log_lockfn lock;
	int level;
	struct callback cb[MAX_CALLBACKS];
} L;

static const char *level_strings[] = { "TRACE", "DEBUG", "INFO",
				       "WARN",	"ERROR", "FATAL" };

#ifdef LOG_USE_COLOR
static const char *level_colors[] = { "\x1b[94m", "\x1b[36m", "\x1b[32m",
				      "\x1b[33m", "\x1b[31m", "\x1b[35m" };
#endif

static void lock(void)
{
	if (L.lock) {
		L.lock(true, L.data);
	}
}

static void unlock(void)
{
	if (L.lock) {
		L.lock(false, L.data);
	}
}

const char *log_level_string(const int level)
{
	return level_strings[level];
}

int log_level_from_lowercase(enum log_level *level, const char *name)
{
	if (!name) {
		return -1;
	}
	if (strncmp(name, "trace", 5) == 0) {
		*level = LOG_TRACE;
		return 0;
	}
	if (strncmp(name, "debug", 5) == 0) {
		*level = LOG_DEBUG;
		return 0;
	}
	if (strncmp(name, "info", 4) == 0) {
		*level = LOG_INFO;
		return 0;
	}
	if (strncmp(name, "warn", 4) == 0) {
		*level = LOG_WARN;
		return 0;
	}
	if (strncmp(name, "error", 5) == 0) {
		*level = LOG_ERROR;
		return 0;
	}
	if (strncmp(name, "fatal", 5) == 0) {
		*level = LOG_FATAL;
		return 0;
	}

	return -1;
}

void log_set_lock(const log_lockfn fn, void *data)
{
	L.lock = fn;
	L.data = data;
}

void log_set_level(const int level)
{
	L.level = level;
}

int log_add_callback(const log_logfn fn, void *data, const int level)
{
	for (int i = 0; i < MAX_CALLBACKS; i++) {
		if (!L.cb[i].fn) {
			L.cb[i] = (struct callback){ fn, data, level };
			return 0;
		}
	}
	return -1;
}

static void init_event(struct log_event *ev, void *data)
{
	if (!ev->time) {
		const time_t t = time(NULL);
		ev->time = localtime(&t);
	}
	ev->data = data;
}

void log_vlog(const int level, const char *fmt, va_list ap)
{
	struct log_event ev = {
		.fmt = fmt,
		.level = level,
	};

	lock();

	for (int i = 0; i < MAX_CALLBACKS && L.cb[i].fn; i++) {
		const struct callback *cb = &L.cb[i];
		if (level >= cb->level) {
			init_event(&ev, cb->data);
			va_copy(ev.ap, ap);
			cb->fn(&ev);
			va_end(ev.ap);
		}
	}

	unlock();
}

void log_log(const int level, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	log_vlog(level, fmt, ap);
	va_end(ap);
}