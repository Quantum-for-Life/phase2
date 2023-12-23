#include "data.h"

#include "test.h"

#ifndef TEST_DATA_SRC
#error
#endif

#include "test_data.h"

int
test_data_open(void);

int
test_data_multidet(void);

int
main(void)
{
	if (test_data_open() < 0) {
		TEST_FAIL("data_open");
		goto error;
	}

	if (test_data_multidet() < 0) {
		TEST_FAIL("data_multidet");
		goto error;
	}

	return 0;
error:
	return -1;
}