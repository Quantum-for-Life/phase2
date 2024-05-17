#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "circ.h"
#include "data.h"
#include "log.h"
#include "qreg.h"

static struct opt {
	const char *filename;
} OPT;

void
opt_help_page(int argc, char **argv);
int
opt_parse(int argc, char **argv);

static int MAIN_RET = 0;

#define ABORT_ON_ERROR(...)                                                    \
	do {                                                                   \
		log_error(__VA_ARGS__);                                        \
		MAIN_RET = EXIT_FAILURE;                                       \
		goto error;                                                    \
	} while (0)

int
run_circuit(data_id fid);

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
		ABORT_ON_ERROR("number of MPI ranks must be a power of two");
	log_info("*** Init ***");
	log_info("MPI num_ranks: %d", num_ranks);
	log_info("This is rank no. %d", rank);

	if (opt_parse(argc, argv) < 0)
		exit(EXIT_FAILURE);

	data_id fid = data_open(OPT.filename);
	if (fid == DATA_INVALID_FID)
		ABORT_ON_ERROR("cannot process input data");

	log_info("*** Circuit ***");
	log_info("Floating point precision: %d", QREG_PREC);
	log_info("QDRIFT >>>");
	if (run_circuit(fid) < 0) {
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
	fprintf(stderr,
		"usage: %s SIMUL_FILE NUM_SAMPLES STEP_SIZE DEPTH"
		"\n\n",
		argv[0]);
}

int
opt_parse(int argc, char **argv)
{
	if (argc < 2) {
		opt_help_page(argc, argv);
		return -1;
	}
	OPT.filename = argv[1];

	return 0;
}

int
run_circuit(data_id fid)
{
	int rc = 0;

	struct circ_data rd;
	if (circ_data_init(&rd, fid) < 0)
		goto error;
	log_info("num_samples: %zu", rd.num_samples);
	log_info("step_size: %f", rd.step_size);
	log_info("depth: %zu", rd.depth);
	if (circ_simulate(&rd) < 0)
		goto error;
	data_trotter_write_values(fid, rd.samples, rd.num_samples);
	goto cleanup;
error:
	rc = -1;
cleanup:
	circ_data_destroy(&rd);

	return rc;
}
