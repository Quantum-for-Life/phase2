#include <math.h>

#include "circ.h"
#include "rayon.h"
#include "data.h"

#define MARGIN (0.005)
#define DATA_FILE DATA_DIR "/case-rand/case-6669b85f.h5"

static double ref_values[][2] = {
	{ 0.12518097361769823, -0.10170409677999773 },
	{ 0.17670425340124754, -0.017778613702941905 },
	{ -0.2884917311702638, -0.04210083472168025 },
	{ 0.11701076895595565, -0.10895444476284752 },
	{ -0.20747966038006477, -0.040093232106182985 },
	{ 0.3388401889521972, 0.020876330290626258 },
	{ -0.3218496848447874, -0.04402145588914208 },
	{ 0.3388401889521972, 0.020876330290626258 }
};

int main(void)
{
	struct data dat;
	data_init(&dat);
	{
		data_id fid = data_file_open(DATA_FILE);
		data_parse(&dat, fid);
		data_file_close(fid);
	}

	{
		struct circ_env env;
		circ_env_init(&env);
		{
			struct rayon_data rd;
			rayon_data_init(&rd);
			rayon_data_from_data(&rd, &dat);
			rayon_simulate(&env, &rd);
			rayon_data_write_times(&dat.time_series, &rd.times);
			rayon_data_destroy(&rd);
		}
		circ_env_destroy(&env);
	}

	for (size_t i = 0; i < dat.time_series.num_steps; i++) {
		double val;
		val = dat.time_series.values[2 * i];
		if (isnan(val))
			return -1;

		double diff = fabs(val - ref_values[i][0]);
		if (diff > MARGIN)
			return -1;

		val = dat.time_series.values[2 * i + 1];
		if (isnan(val))
			return -1;
		diff = fabs(val - ref_values[i][1]);
		if (diff > MARGIN)
			return -1;
	}
	data_destroy(&dat);

	return 0;
}
