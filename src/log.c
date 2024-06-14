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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mpi.h"

#include "log.h"

#define MAX_CALLBACKS (32)

struct log_event {
	int	   level;
	struct tm *time;

	va_list	    ap;
	const char *fmt;

	void *data;
};

struct callback {
	void  (*fn)(struct log_event *ev);
	void *data;
	int   level;
};

static struct {
	void	       *data;
	int		level;
	struct callback cb[MAX_CALLBACKS];
} L;

static const char *level_strings[] = { "TRACE", "DEBUG", "INFO", "WARN",
	"ERROR", "FATAL" };

const char *log_level_string(const int level)
{
	return level_strings[level];
}

int log_level_from_lowercase(enum log_level *level, const char *name)
{
	if (!name)
		return -1;

	int rt = 0;
	switch (name[0]) {
	case 't':
	case 'T':
		*level = LOG_TRACE;
		break;
	case 'd':
	case 'D':
		*level = LOG_DEBUG;
		break;
	case 'i':
	case 'I':
		*level = LOG_INFO;
		break;
	case 'w':
	case 'W':
		*level = LOG_WARN;
		break;
	case 'e':
	case 'E':
		*level = LOG_ERROR;
		break;
	case 'f':
	case 'F':
		*level = LOG_FATAL;
		break;
	default:
		rt = -1;
	}

	return rt;
}

void log_set_level(const int level)
{
	L.level = level;
}

int log_add_callback(
	void (*fn)(struct log_event *ev), void *data, const int level)
{
	for (int i = 0; i < MAX_CALLBACKS; i++)
		if (!L.cb[i].fn) {
			L.cb[i] = (struct callback){ fn, data, level };
			return 0;
		}

	return -1;
}

static void init_event(struct log_event *ev, void *data)
{
	if (!ev->time) {
		const time_t t = time(NULL);
		ev->time       = localtime(&t);
	}
	ev->data = data;
}

void log_vlog(const int level, const char *fmt, va_list ap)
{
	struct log_event ev = {
		.fmt   = fmt,
		.level = level,
	};

	for (int i = 0; i < MAX_CALLBACKS && L.cb[i].fn; i++) {
		const struct callback *cb = &L.cb[i];
		if (level >= cb->level) {
			init_event(&ev, cb->data);
			va_copy(ev.ap, ap);
			cb->fn(&ev);
			va_end(ev.ap);
		}
	}
}

void log_log(const int level, const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	log_vlog(level, fmt, ap);
	va_end(ap);
}

void log_callback(struct log_event *ev)
{
	int initialized, finalized;
	MPI_Initialized(&initialized);
	MPI_Finalized(&finalized);
	if (initialized && !finalized) {
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank > 0)
			return;
	}

	char  buf[64];
	FILE *fd = ev->data;
	buf[strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", ev->time)] = '\0';
	fprintf(fd, "%s %-5s ", buf, log_level_string(ev->level));
	vfprintf(fd, ev->fmt, ev->ap);
	fprintf(fd, "\n");
	fflush(fd);
}

int log_init(void)
{
	enum log_level lvl;
	const char    *lvl_str = getenv(PHASE2_LOG_ENVVAR);
	if (!lvl_str || log_level_from_lowercase(&lvl, lvl_str) < 0)
		lvl = LOG_ERROR;

	log_set_level(lvl);
	log_add_callback(log_callback, stderr, lvl);

	return 0;
}
