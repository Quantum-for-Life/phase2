#ifndef TEST_DATA_H
#define TEST_DATA_H

#include <stddef.h>
#include <stdio.h>

#include "phase2/world.h"

#include "test.h"

/*
 * Copy a committed fixture to /tmp so per-step writes
 * issued by trott_simul / qdrift_simul / etc. do not
 * mutate the on-disk reference file.  Rank 0 does the
 * copy; followers no-op.  Caller deletes the temporary
 * after use.
 */
static inline void test_fixture_copy(const char *src, const char *dst)
{
	struct world_info wd;
	world_info(&wd);
	if (wd.rank != 0)
		return;
	FILE *in = fopen(src, "rb");
	TEST_ASSERT(in != NULL, "test_fixture_copy: open src %s", src);
	FILE *out = fopen(dst, "wb");
	TEST_ASSERT(out != NULL, "test_fixture_copy: open dst %s", dst);
	char buf[4096];
	size_t n;
	while ((n = fread(buf, 1, sizeof buf, in)) > 0)
		TEST_ASSERT(fwrite(buf, 1, n, out) == n,
			"test_fixture_copy: short write to %s", dst);
	fclose(in);
	fclose(out);
}

#define NUM_TEST_FILES (2)
[[maybe_unused]] static struct test_data {
	const char *filename;
	size_t num_qubits;
	size_t num_terms;
	size_t num_dets;
	size_t num_steps;
	double norm;
} TEST_DATA[NUM_TEST_FILES] = {
	{ .filename = PH2_TESTDIR "/data/H2O_CAS56.h5",
		.num_qubits = 10,
		.num_terms = 251,
		.num_dets = 1,
		.num_steps = 4,
		.norm = 0.07170948 },
	{ .filename = PH2_TESTDIR "/data/case-d9f603dc.h5_solved",
		.num_qubits = 3,
		.num_terms = 10,
		.num_dets = 3,
		.num_steps = 16,
		.norm = 0.18535287 }
};

#endif // TEST_DATA_H
