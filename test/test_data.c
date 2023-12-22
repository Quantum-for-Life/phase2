#ifndef TEST_DATA_FILE
#error
#endif

int
test_data_open(void);

int
test_data_open_nonexist(void);

int
main(void)
{
	if (test_data_open() < 0)
		goto error;
	if (test_data_open_nonexist() < 0)
		goto error;

	return 0;
error:
	return -1;
}