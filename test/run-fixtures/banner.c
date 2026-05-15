/*
 * Synthetic test fixture: writes a distinctive banner
 * to both stdout and stderr, then exits 1.  Exercises
 * the runner's failure-dump output (captured streams
 * surfaced under `---- t-banner stdout ----` /
 * `---- t-banner stderr ----` headers).
 */
#include <stdio.h>

int main(void)
{
	fputs("BANNER_STDOUT_LINE\n", stdout);
	fputs("BANNER_STDERR_LINE\n", stderr);
	return 1;
}
