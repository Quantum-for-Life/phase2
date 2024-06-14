#include "test.h"

int test_data_open(void);

int test_data_multidet(void);

int test_data_hamil(void);

int test_data_trott_steps(void);

int main(void)
{
	if (test_data_open() < 0) {
		TEST_FAIL("data_open");
		goto err;
	}
	if (test_data_multidet() < 0) {
		TEST_FAIL("data_multidet");
		goto err;
	}
	if (test_data_hamil() < 0) {
		TEST_FAIL("data_hamil");
		goto err;
	}
	if (test_data_trott_steps() < 0) {
		TEST_FAIL("data_trott_steps");
		goto err;
	}

	return 0;
err:
	return -1;
}
