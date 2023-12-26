#include "algos/linen.h"
#include "algos/rayon.h"

int
run_linen(void)
{
	return linen_simulate();
}

int
run_rayon(data_id fid)
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
