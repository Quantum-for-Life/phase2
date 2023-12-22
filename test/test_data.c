
int
test_data_open(void);

int
main(void)
{
	if (test_data_open() < 0)
		goto error;

	return 0;
error:
	return -1;
}