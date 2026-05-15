/*
 * Synthetic test fixture: sleeps 250 ms, then exits 0.
 * Four copies run in parallel exercise the runner's
 * --jobs concurrency: serial would take ~1.0 s wall,
 * parallel ~0.25 s.
 */
#define _POSIX_C_SOURCE 200809L

#include <time.h>

int main(void)
{
	const struct timespec ts = { 0, 250000000 };
	nanosleep(&ts, NULL);
	return 0;
}
