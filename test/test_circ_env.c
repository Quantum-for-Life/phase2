#include <threads.h>

#include "circ.h"

#include "test.h"

#define NUM_THREADS (8)

int
init_par(void *)
{
	return circ_initialize();
}

int
main(void)
{
	volatile int bar;
	thrd_t	     th[NUM_THREADS];
	int	     th_ret[NUM_THREADS];

	for (size_t i = 0; i < NUM_THREADS; i++) {
		bar	  = 0;
		th_ret[i] = -1;
		if (thrd_create(&th[i], init_par, NULL) != thrd_success) {
			TEST_FAIL("creating thread %zu", i);
			goto error;
		}
	}

	for (size_t i = 0; i < NUM_THREADS; i++) {
		bar = 0;
		if (thrd_join(th[i], th_ret + i) != thrd_success) {
			TEST_FAIL("joining thread %zu", i);
			goto error;
		}
	}

	size_t num_zeroes = 0, num_ones = 0;
	for (size_t i = 0; i < NUM_THREADS; i++) {
		bar = 0;
		if (th_ret[i] < 0) {
			TEST_FAIL("thread[%zu] returned %d", i, th_ret[i]);
			goto error;
		}

		if (th_ret[i] == 0)
			num_zeroes++;
		if (th_ret[i] == 1)
			num_ones++;
	}

	if (num_zeroes != 1) {
		TEST_FAIL("init should happen exactly once");
		goto error;
	}

	if (num_ones != NUM_THREADS - 1) {
		TEST_FAIL("the rest should be passes");
		goto error;
	}

	return 0;
error:
	return -1;
}