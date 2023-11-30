#ifndef TEST_H
#define TEST_H

#include <stdio.h>

#define TEST(name, ...)	\
int name(__VA_ARGS__) { \
	const char* __test_name = #name; \
	int __test_rc;

#define TEST_ASSERT(cond, ...)	\
	if (!(cond)) { \
		fprintf(stderr, "FAIL: %s\n", __test_name); \
		fprintf(stderr, "-- %s:%d\n-- ",  __FILE__, __LINE__);       \
		fprintf(stderr, __VA_ARGS__);	\
		fprintf(stderr, "\n");	\
		goto __test_error;	\
	}

#define TEST_FINALIZE	\
	__test_rc = 0;	\
	goto __test_exit;	\
	__test_error:	\
	__test_rc = -1;	\
	__test_exit:

#define TEST_END	\
	return	__test_rc; \
}

#define TEST_CASE(exp) TEST_ASSERT(exp == 0, "%s", #exp)

#endif // TEST_H
