#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "circ.h"
#include "data2.h"
#include "log.h"

struct opt {
	const char *dat_filename;
	size_t	    num_steps;
};

static struct opt opt;

void opt_help_page(int argc, char **argv);
int  opt_parse(struct opt *o, int argc, char **argv);

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
	log_info("Num_steps: %zu", opt.num_steps);
	if (circuit_run(fid, opt.num_steps) < 0) {
		log_error("Failure: simulation error");
		goto error;
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
	fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5] [NUM_STEPS]\n\n",
		argv[0]);
}

int opt_parse(struct opt *o, int argc, char **argv)
{
	/* Parse command line arguments. */
	if (argc < 3) {
		opt_help_page(argc, argv);
		return -1;
	}

	o->dat_filename		     = argv[1];
	unsigned long long num_steps = strtoull(argv[2], NULL, 10);
	if (num_steps == 0) {
		fprintf(stderr, "Wrong number of Trotter steps\n");
		return -1;
	}
	o->num_steps = num_steps;

	return 0;
}
