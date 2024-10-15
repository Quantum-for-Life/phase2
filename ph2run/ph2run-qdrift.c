#define _XOPEN_SOURCE 600
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

#include "phase2/circ.h"
#include "phase2/world.h"

#define WD_SEED UINT64_C(0xd871e5d39fc0222d)
static struct world WD;

static struct opt {
	const char *filename;
} OPT;

void opt_help_page(int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s SIMUL_FILE\n\n", argv[0]);
}

int opt_parse(int argc, char **argv)
{
	if (argc < 2) {
		opt_help_page(argc, argv);
		return -1;
	}
	OPT.filename = argv[1];

	return 0;
}

int run_circuit(data_id fid)
{
	int rt = -1; /* Return value */

	struct timespec t1, t2;
	double t_tot;

	struct circ_qdrift_data rd;
	if (circ_qdrift_data_init(&rd, fid) < 0)
		goto exit_qdrift_data;
	log_info("num_samples: %zu", rd.num_samples);
	log_info("step_size: %f", rd.step_size);
	log_info("depth: %zu", rd.depth);

	clock_gettime(CLOCK_REALTIME, &t1);
	if (circ_qdrift_simulate(&rd) < 0)
		goto exit_qdrift_simulate;
	clock_gettime(CLOCK_REALTIME, &t2);
	t_tot = (double)(t2.tv_sec - t1.tv_sec) +
		(double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	if (data_circ_qdrift_write_values(fid, rd.samples, rd.num_samples) < 0)
		goto exit_qdrift_write;

	rt = 0; /* Success. */

exit_qdrift_write:
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_samples,step_size,depth,"
		 "n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%zu,%.3f,%d,%.3f", rd.hamil.num_qubits,
		rd.hamil.num_terms, rd.multidet.num_dets, rd.num_samples,
		rd.step_size, rd.depth, WD.size, t_tot);

exit_qdrift_simulate:
	circ_qdrift_data_destroy(&rd);
exit_qdrift_data:

	return rt;
}

int main(int argc, char **argv)
{
	int rt = -1; /* Return value. */

	if (world_init(&argc, &argv, WD_SEED) != WORLD_READY)
		goto exit_world_init;
	world_info(&WD);

	unsigned int num_ranks = WD.size;
	if (num_ranks == 0 || (num_ranks & (num_ranks - 1)) != 0) {
		log_error("number of MPI ranks (%u) "
			  "must be a power of two.",
			num_ranks);
		goto exit_num_ranks;
	}

	log_info("*** Init ***");
	log_info("MPI num_ranks: %d", num_ranks);
	log_info("This is rank no. %d", WD.rank);

	if (opt_parse(argc, argv) < 0) {
		log_error("can't parse program arguments");
		goto exit_opt_parse;
	}

	data_id fid = data_open(OPT.filename);
	if (fid == DATA_INVALID_FID) {
		log_error("cannot process input data");
		goto exit_data_open;
	}

	log_info("*** Circuit: qDRIFT >>> ***");
	if (run_circuit(fid) < 0) {
		log_error("Failure: simulation error");
		goto exit_run_circuit;
	}

	rt = 0; /* Success. */

exit_run_circuit:
	log_info("Shut down simulation environment");
	data_close(fid);
exit_data_open:
exit_opt_parse:
exit_num_ranks:
exit_world_init:
	world_destroy();

	return rt;
}
