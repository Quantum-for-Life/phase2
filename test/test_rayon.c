#include <stdio.h>
#include <math.h>

#include "circ.h"
#include "rayon.h"
#include "data.h"

#define MARGIN (0.001)

static double ref_values[][2] = { { 0.8237180013918515, 0.056881902063548984 },
				  { 0.9188847808451465, 0.041193087387450644 },
				  { 0.9792943197305308, 0.02162383108455912 },
				  { 0.8752221697187768, 0.04962192917736108 },
				  { 0.20553290304802077, 0.056359153263639945 },
				  { 0.27247599116749827, 0.06202025540798886 },
				  { 0.9947965503140022, 0.010943075593945473 },
				  { 0.7014724677600823, 0.06729436930500987 },
				  { 0.20553290304802077, 0.056359153263639945 },
				  { 0.34276716474828706, 0.06651455873389611 },
				  { 0.7014724677600823, 0.06729436930500987 },
				  { 0.9188847808451465, 0.041193087387450644 },
				  { 0.4152558581177978, 0.06968668379016273 },
				  { 0.4152558581177978, 0.06968668379016273 },
				  { 0.4152558581177978, 0.06968668379016273 },
				  { 0.953814854174465, 0.03178731460121187 } };

int main(void)
{
	int rc = 0;

	struct data dat;
	data_init(&dat);
	data_id fid = data_file_open(DATA_DIR "/case-rand/case-26bcd92a.h5");
	data_parse(&dat, fid);
	data_file_close(fid);

	struct circ_env env;
	circ_env_init(&env);
	struct rayon_data ct_dat;
	rayon_data_init(&ct_dat);
	rayon_data_from_data(&ct_dat, &dat);
	rayon_simulate(env, &ct_dat, &dat.time_series);
	circ_env_destroy(&env);

	for (size_t i = 0; i < dat.time_series.num_steps; i++) {
		if (fabs(dat.time_series.values[2 * i] - ref_values[i][0]) >
		    MARGIN) {
			return -1;
		}
		if (fabs(dat.time_series.values[2 * i + 1] - ref_values[i][1]) >
		    MARGIN) {
			return -1;
		}
	}

	rayon_data_destroy(&ct_dat);
	return rc;
}
