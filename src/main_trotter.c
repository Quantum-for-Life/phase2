#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "circ_trotter.h"
#include "data.h"
#include "log.h"
#include "qreg.h"

static struct opt {
	const char *filename;
	size_t	    num_steps;
} OPT;

void
opt_help_page(int argc, char **argv);
int
opt_parse(int argc, char **argv);

static int MAIN_RET = 0;

#define ABORT_ON_ERROR(...)                                                    \
	{                                                                      \
		log_error(__VA_ARGS__);                                        \
		MAIN_RET = EXIT_FAILURE;                                       \
		goto error;                                                    \
	}

int
run_circuit(data_id fid, size_t num_steps);

int
main(int argc, char **argv)
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

	data_id fid = data_open(OPT.filename);
	if (fid == DATA_INVALID_FID)
		ABORT_ON_ERROR("cannot process input data")

	log_info("*** Circuit ***");
	log_info("Floating point precision: %d", QREG_PREC);
	log_info("Num_steps: %zu", OPT.num_steps);
	if (run_circuit(fid, OPT.num_steps) < 0) {
		log_error("Failure: simulation error");
		goto error;
	}

	data_close(fid);
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

void
opt_help_page(int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s [SIMUL_FILE] [NUM_STEPS]\n\n", argv[0]);
}

int
opt_parse(int argc, char **argv)
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

int
run_circuit(data_id fid, size_t num_steps)
{
	int rc = 0;

	struct circ_trotter_data rd;
	circ_trotter_data_init(&rd, num_steps);
	if (circ_trotter_data_from_file(&rd, fid) < 0)
		goto error;
	if (circ_trotter_simulate(&rd) < 0)
		goto error;
	data_circ_trotter_write_values(fid, rd.trott_steps, num_steps);
	goto cleanup;
error:
	rc = -1;
cleanup:
	circ_trotter_data_destroy(&rd);

	return rc;
}
