#include "circ.h"
#include "rayon.h"
#include "linen.h"

int run_linen(struct circ_env *env, struct data *dat)
{
	(void)(dat);

	return linen_simulate(env);
}

int run_rayon(struct circ_env *env, struct data *dat)
{
	int rc = 0;

	struct rayon_data rd;
	rayon_data_init(&rd);
	if (rayon_data_from_data(&rd, dat) < 0)
		goto error;
	if (rayon_simulate(env, &rd) < 0)
		goto error;
	rayon_data_write_times(&dat->time_series, &rd.times);
	goto cleanup;
error:
	rc = -1;
cleanup:
	rayon_data_destroy(&rd);

	return rc;
}
