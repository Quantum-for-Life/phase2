/* -------------------------------------------------------------------------- *
 * bench.c: Simple benchmarking framework.                                    *
 *                                                                            *
 * Copyright 2024 ⧉⧉⧉                                                         *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify it    *
 * under the terms of the GNU General Public License as published by the      *
 * Free Software Foundation, either version 3 of the License, or (at your     *
 * option) any later version.                                                 *
 *                                                                            *
 * This program is distributed in the hope that it will be useful, but        *
 * WITHOUT ANY WARRANTY* without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with this program.  If not, see <https://www.gnu.org/licenses/>.           *
 * -------------------------------------------------------------------------- */
#define _XOPEN_SOURCE 700

#include <stddef.h>
#include <time.h>

#include "bench.h"

#define REPS_INIT (99)

int bench_mark(struct bench *b, size_t reps, int (*op)(void *), void *data)
{
	int rt = 0;

	struct timespec t1, t2;
	volatile size_t r;

	b->nanos = 0;
	b->reps = 0;
	/* First, warm the cache up */
	for (volatile size_t i = 0; i < REPS_INIT; i++)
		if ((rt = op(data)) != 0)
			return rt;

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
	for (r = 0; r < reps; r++)
		if ((rt = op(data)) != 0)
			break;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t2);

	b->nanos = 1000000000UL * (t2.tv_sec - t1.tv_sec) +
		   (t2.tv_nsec - t1.tv_nsec);
	b->reps = r;

	return rt;
}

double bench_msrep(struct bench b)
{
	return (double)b.nanos / b.reps / 1.0e6;
}
