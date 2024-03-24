#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "circ.h"
#include "data2.h"
#include "log.h"

#define PH2RUN_DEFAULT_H5FILE "simul.h5"

#define PH2RUN_SILK_DEFAULT_NUM_STEPS (8)

enum opt_cicuit {
	OPT_CICUIT_SILK,
};

struct args_silk {
	size_t num_steps;
};

struct opt {
	enum opt_cicuit cicuit;
	const char     *dat_filename;
	union {
		struct args_silk silk;
	} circuit_args;
};

static struct opt opt;

void opt_help_page(int argc, char **argv);
int  opt_parse(struct opt *o, int argc, char **argv);

/* Runners */
int run_silk(data2_id fid, size_t num_steps);

static int MAIN_RET = 0;

#define ABORT_ON_ERROR(...)                                                    \
	{                                                                      \
		log_error(__VA_ARGS__);                                        \
		MAIN_RET = EXIT_FAILURE;                                       \
		goto error;                                                    \
	}

int main(int argc, char **argv)
{
	int initialized;
	MPI_Initialized(&initialized);
	if (!initialized && MPI_Init(&argc, &argv) != MPI_SUCCESS)
		return -1;

	/* Initiallize logging */
	if (log_init() < 0)
		exit(EXIT_FAILURE);

	int rank, num_ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if ((num_ranks & (num_ranks - 1)) != 0)
		ABORT_ON_ERROR("number of MPI ranks must be a power of two");

	log_info("*** Init ***");
	log_info("MPI num_ranks: %d", num_ranks);
	log_info("This is rank no. %d", rank);

	/* Parse command line arguments. */
	if (opt_parse(&opt, argc, argv) < 0)
		exit(EXIT_FAILURE);

	circ_initialize();

	data2_id fid = data2_open(opt.dat_filename);
	if (fid == DATA2_INVALID_FID)
		ABORT_ON_ERROR("cannot process input data");

	log_info("*** Circuit ***");
	switch (opt.cicuit) {
	case OPT_CICUIT_SILK:
		log_info("Circuit: silk");
		log_info("Num_steps: %zu", opt.circuit_args.silk.num_steps);
		if (run_silk(fid, opt.circuit_args.silk.num_steps) < 0) {
			log_error("Failure: simulation error");
			goto error;
		}
		break;
	}

	data2_close(fid);

	log_info("Shut down simulation environment");
	return EXIT_SUCCESS;

error:
	log_error("Shut down simulation environment");

	return MAIN_RET;
}

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
