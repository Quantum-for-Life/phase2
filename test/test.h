#ifndef TEST_H
#define TEST_H

#include <stdbool.h>
#include <stdio.h>

static _Atomic bool TEST_RT = true;

/*
 * __VA_OPT__ has been added to C23. It is also part of GNU C and is
 * supported by gcc.
 */
#define TEST_FAIL(fmt, ...)                                                    \
	({                                                                     \
		fprintf(stderr, "%s:%d FAIL \"" fmt "\"\n", __FILE__,          \
			__LINE__ __VA_OPT__(, ) __VA_ARGS__);                  \
		TEST_RT = false;                                               \
	})

#define TEST_ASSERT(expr, fmt, ...)                                            \
	({                                                                     \
		if (!(expr))                                                   \
			TEST_FAIL(fmt __VA_OPT__(, ) __VA_ARGS__);             \
	})

#define TEST_STR(a) #a
#define TEST_XSTR(a) TEST_STR(a)

#define TEST_EQ(a, b)                                                          \
	({                                                                     \
		TEST_ASSERT((a) == (b), "expected %s == %s", TEST_XSTR(a),     \
			TEST_XSTR(b));                                         \
	})

static void TEST_MAIN(void);

int main(int argc, char **argv)
{
	(void)argc;
	(void)argv;

	TEST_MAIN();

	return TEST_RT ? 0 : -1;
}

#endif /* TEST_H */
