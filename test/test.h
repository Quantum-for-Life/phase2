#ifndef TEST_H
#define TEST_H

#include <stdio.h>

#define TEST(name, ...)                        \
	static int name(__VA_ARGS__)           \
	{                                      \
		const char *test_name = #name; \
		int test_rc = 0;

#define TEST_ASSERT(exp, ...)                                         \
	if (!(exp)) {                                                 \
		fprintf(stderr, "FAIL: %s\n", test_name);             \
		fprintf(stderr, "-- %s:%d\n-- ", __FILE__, __LINE__); \
		fprintf(stderr, __VA_ARGS__);                         \
		fprintf(stderr, "\n");                                \
		goto test_error;                                      \
	}

#define TEST_FIN(SMNT)  \
	test_rc = 0;    \
	goto test_exit; \
test_error:             \
	test_rc = -1;   \
test_exit:              \
	SMNT;           \
	return test_rc; \
	}

#define TEST_CASE(exp) TEST_ASSERT(exp == 0, "%s", #exp)

#endif // TEST_H
