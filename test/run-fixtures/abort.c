/*
 * Synthetic test fixture: raises SIGABRT.  The runner
 * reports a signal-killed child as `128 + signum`,
 * i.e. exit code 134 here (SIGABRT == 6).
 */
#include <stdlib.h>

int main(void)
{
	abort();
}
