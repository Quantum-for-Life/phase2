#include "algos/linen.h"
#include "algos/rayon.h"
#include "algos/silk.h"

int run_linen(void)
{
	return linen_simulate();
}

int run_rayon(data2_id fid)
{
	int rc = 0;

	struct rayon_data rd;
	rayon_data_init(&rd);
	if (rayon_data_from_data(&rd, fid) < 0)
		goto error;
	if (rayon_simulate(&rd) < 0)
		goto error;
	rayon_data_times_write(fid, &rd.times);
	goto cleanup;
error:
	rc = -1;
cleanup:
	rayon_data_destroy(&rd);

	return rc;
}

int run_silk(data2_id fid, size_t num_steps)
{
	int rc = 0;

	struct silk_data rd;
	silk_data_init(&rd, num_steps);
	if (silk_data_from_data(&rd, fid) < 0)
		goto error;
	if (silk_simulate(&rd) < 0)
		goto error;
	data2_trotter_write_values(fid, rd.trotter_steps, num_steps);
	goto cleanup;
error:
	rc = -1;
cleanup:
	silk_data_destroy(&rd);

	return rc;
}
