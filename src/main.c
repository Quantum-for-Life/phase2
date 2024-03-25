#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "circ.h"
#include "data2.h"
#include "log.h"

static struct opt {
	const char *filename;
	size_t	    num_steps;
} OPT;

void opt_help_page(int argc, char **argv);
int  opt_parse(int argc, char **argv);

static int MAIN_RET = 0;

#define ABORT_ON_ERROR(...)                                                    \
	{                                                                      \
		log_error(__VA_ARGS__);                                        \
		MAIN_RET = EXIT_FAILURE;                                       \
		goto error;                                                    \
	}

int run_circuit(data2_id fid, size_t num_steps);

int main(int argc, char **argv)
{
	int initialized;
	MPI_Initialized(&initialized);
	if (!initialized && MPI_Init(&argc, &argv) != MPI_SUCCESS)
		exit(EXIT_FAILURE);

	if (log_init() < 0)
		exit(EXIT_FAILURE);

	int rank, num_ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if ((num_ranks & (num_ranks - 1)) != 0)
		ABORT_ON_ERROR("number of MPI ranks must be a power of two")
	log_info("*** Init ***");
	log_info("MPI num_ranks: %d", num_ranks);
	log_info("This is rank no. %d", rank);

	if (opt_parse(argc, argv) < 0)
		exit(EXIT_FAILURE);

	data2_id fid = data2_open(OPT.filename);
	if (fid == DATA2_INVALID_FID)
		ABORT_ON_ERROR("cannot process input data")

	log_info("*** Circuit ***");
	log_info("Num_steps: %zu", OPT.num_steps);
	if (run_circuit(fid, OPT.num_steps) < 0) {
		log_error("Failure: simulation error");
		goto error;
	}

	data2_close(fid);
	goto cleanup;

error:
	MAIN_RET = EXIT_FAILURE;
cleanup:

	log_info("Shut down simulation environment");
	int finalized;
	MPI_Finalized(&finalized);
	if (!finalized && MPI_Finalize() != MPI_SUCCESS)
		MAIN_RET = EXIT_FAILURE;

	return MAIN_RET;
}

void opt_help_page(int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5] [NUM_STEPS]\n\n",
		argv[0]);
}

int opt_parse(int argc, char **argv)
{
	if (argc < 3) {
		opt_help_page(argc, argv);
		return -1;
	}

	OPT.filename	 = argv[1];
	size_t num_steps = strtoull(argv[2], NULL, 10);
	if (num_steps == 0) {
		fprintf(stderr, "Wrong number of Trotter steps\n");
		return -1;
	}
	OPT.num_steps = num_steps;

	return 0;
}

int run_circuit(data2_id fid, size_t num_steps)
{
	int rc = 0;

	struct circuit_data rd;
	circuit_data_init(&rd, num_steps);
	if (circuit_data_from_file(&rd, fid) < 0)
		goto error;
	if (circuit_simulate(&rd) < 0)
		goto error;
	data2_trotter_write_values(fid, rd.trotter_steps, num_steps);
	goto cleanup;
error:
	rc = -1;
cleanup:
	circuit_data_destroy(&rd);

	return rc;
}
