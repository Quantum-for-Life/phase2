#include <stdio.h>

#include "argums.h"

void print_help_page(int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5]\n\n", argv[0]);
	fprintf(stderr, "where CIRCUIT is must be one of:\n"
			"\n"
			"    linen\n"
			"    rayon\n"
			"\n"
			"If no simulation input file (HDF5) is specified,\n"
			"the default is " PH2RUN_DEFAULT_H5FILE
			" in the current directory.\n");
	fprintf(stderr, "\n"
			"To set logging level, set environment variable:\n"
			"\n    " PH2RUN_LOG_ENVVAR
			"={trace, debug, info, warn, error, fatal}"
			"\n\n");
}

int argums_parse(struct argums *a, int argc, char **argv)
{
	int rc = 0;

	/* Parse command line arguments. */
	if (argc < 2) {
		print_help_page(argc, argv);
		return -1;
	}
	a->dat_filename = PH2RUN_DEFAULT_H5FILE;
	if (argc >= 3)
		a->dat_filename = argv[2];

	return rc;
}
