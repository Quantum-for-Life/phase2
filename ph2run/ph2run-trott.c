#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

#include "phase2/circ.h"
#include "phase2/world.h"

static struct world WD;

static struct opt {
	const char *filename;
	size_t num_steps;
} OPT;

void opt_help_page(int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s [SIMUL_FILE] [NUM_STEPS]\n\n", argv[0]);
}

int opt_parse(int argc, char **argv)
{
	if (argc < 3) {
		opt_help_page(argc, argv);
		return -1;
	}

	OPT.filename = argv[1];
	size_t num_steps = strtoull(argv[2], NULL, 10);
	if (num_steps == 0) {
		fprintf(stderr, "Wrong number of Trotter steps\n");
		return -1;
	}
	OPT.num_steps = num_steps;

	return 0;
}

int run_circuit(data_id fid, size_t num_steps)
{
	int rt = -1; /* Return value */

	struct timespec t1, t2;
	double t_tot;

	struct circ_trott_data rd;
	if (circ_trott_data_init_from_file(&rd, num_steps, fid) < 0)
		goto exit_trott_data;

	clock_gettime(CLOCK_REALTIME, &t1);
	if (circ_trott_simulate(&rd) < 0)
		goto exit_trott_simulate;

	clock_gettime(CLOCK_REALTIME, &t2);
	t_tot = (double)(t2.tv_sec - t1.tv_sec) +
		(double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	if (data_circ_trott_write_values(fid, rd.trott_steps, num_steps) < 0)
		goto exit_trott_write;

	rt = 0; /* Success. */

exit_trott_write:
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%d,%.3f",
		rd.hamil.num_qubits, rd.hamil.num_terms,
		rd.multidet.num_dets, rd.num_trott_steps, WD.size, t_tot);

exit_trott_simulate:
	circ_trott_data_destroy(&rd);

exit_trott_data:

	return rt;
}

int main(int argc, char **argv)
{
	int rt = -1; /* Return value. */

	if (world_init(&argc, &argv) != WORLD_READY)
		goto exit_world_init;
	world_info(&WD);

	unsigned int num_ranks = WD.size;
	if (num_ranks == 0 || (num_ranks & (num_ranks - 1)) != 0) {
		log_error("number of MPI ranks (%u) "
			"must be a power of two.", num_ranks);
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

	log_info("*** Circuit: trott ***");
	log_info("Num_steps: %zu", OPT.num_steps);
	if (run_circuit(fid, OPT.num_steps) < 0) {
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
	world_fin();

	return rt;
}
