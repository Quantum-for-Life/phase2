#ifndef TEST_DATA_H
#define TEST_DATA_H

#define NUM_TEST_FILES (2)
static struct test_data {
	const char *filename;
	size_t	    num_qubits;
	size_t	    num_terms;
	size_t	    num_dets;
	size_t	    num_steps;
	double	    norm;
} TEST_DATA[NUM_TEST_FILES] = {
	{ .filename	    = TEST_DATA_SRC "/H2O_CAS56.h5",
		.num_qubits = 10,
		.num_terms  = 251,
		.num_dets   = 1,
		.num_steps  = 4,
		.norm	    = 0.07170948 },
	{ .filename	    = TEST_DATA_SRC "/case-d9f603dc.h5_solved",
		.num_qubits = 3,
		.num_terms  = 10,
		.num_dets   = 3,
		.num_steps  = 16,
		.norm	    = 0.18535287 }
};

#endif // TEST_DATA_H
