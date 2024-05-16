#ifndef TEST_H
#define TEST_H

#include <stdio.h>

#define TEST_FAIL(...)                                                         \
	do {                                                                   \
		fprintf(stderr, "FAILED %s:%d \"", __FILE__, __LINE__);        \
		fprintf(stderr, __VA_ARGS__);                                  \
		fprintf(stderr, "\"\n");                                       \
	} while (0)

#endif // TEST_H
