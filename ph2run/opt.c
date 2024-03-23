#include <stdio.h>

#include "opt.h"

#include <stdlib.h>
#include <string.h>

void opt_help_page(int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5] [SIMUL_ARGS]\n\n",
		argv[0]);
	fprintf(stderr, "where CIRCUIT is must be one of:\n"
			"\n"
			"    linen\n"
			"    rayon\n"
			"    silk\n"
			"\n"
			"If no simulation input file (HDF5) is specified,\n"
			"the default is " PH2RUN_DEFAULT_H5FILE " in the "
			"current directory.\n");
	fprintf(stderr, "\n"
			"To set logging level, set environment variable:\n"
			"\n    " PH2RUN_LOG_ENVVAR "={trace, debug, info, "
			"warn, error, fatal}"
			"\n\n");
}

int opt_parse(struct opt *o, int argc, char **argv)
{
	/* Parse command line arguments. */
	if (argc < 2) {
		opt_help_page(argc, argv);
		return -1;
	}

	if (strncmp(argv[1], "silk", 4) == 0) {
		o->cicuit = OPT_CICUIT_SILK;
	} else {
		fprintf(stderr, "No circ named %s\n", argv[1]);
		return -1;
	}

	o->dat_filename = PH2RUN_DEFAULT_H5FILE;
	if (argc >= 3)
		o->dat_filename = argv[2];

	if (o->cicuit == OPT_CICUIT_SILK) {
		if (argc < 4)
			o->circuit_args.silk.num_steps =
				PH2RUN_SILK_DEFAULT_NUM_STEPS;
		else {
			unsigned long long num_steps =
				strtoull(argv[3], NULL, 10);
			if (num_steps == 0) {
				fprintf(stderr,
					"Wrong number of Trotter steps\n");
				return -1;
			}
			o->circuit_args.silk.num_steps = num_steps;
		}
	}

	return 0;
}
