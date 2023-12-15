#include <threads.h>

#include "circ.h"

#include "test.h"

#define NUM_THREADS (8)

int init_par(void *)
{
	return circ_initialize();
}

TEST(init_parallel, void)
{
	volatile int bar;
	thrd_t th[NUM_THREADS];
	int th_ret[NUM_THREADS];

	for (size_t i = 0; i < NUM_THREADS; i++) {
		bar = 0;
		th_ret[i] = -1;
		TEST_ASSERT(thrd_create(&th[i], init_par, NULL) == thrd_success,
			    "creating thread %zu", i);
	}

	for (size_t i = 0; i < NUM_THREADS; i++) {
		bar = 0;
		TEST_ASSERT(thrd_join(th[i], th_ret + i) == thrd_success,
			    "joining thread %zu", i);
	}

	size_t num_zeroes = 0, num_ones = 0;
	for (size_t i = 0; i < NUM_THREADS; i++) {
		bar = 0;
		TEST_ASSERT(th_ret[i] >= 0, "thread[%zu] returned %d", i,
			    th_ret[i]);

		if (th_ret[i] == 0)
			num_zeroes++;
		if (th_ret[i] == 1)
			num_ones++;
	}

	TEST_ASSERT(num_zeroes == 1, "init should happen exactly once");
	TEST_ASSERT(num_ones == NUM_THREADS - 1, "the rest should be passes");

	TEST_FIN(circ_shutdown())
}

int main(void)
{
	return init_parallel();
}