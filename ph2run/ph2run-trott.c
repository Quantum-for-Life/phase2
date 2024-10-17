#define _XOPEN_SOURCE 600
#include "c23_compat.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

#include "phase2/circ.h"
#include "phase2/world.h"

#define WD_SEED UINT64_C(0xd326119d4859ebb2)
static struct world WD;

static struct opt {
	const char *filename;
	size_t nsteps;
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
	size_t nsteps = strtoull(argv[2], NULL, 10);
	if (nsteps == 0) {
		fprintf(stderr, "Wrong number of Trotter steps\n");
		return -1;
	}
	OPT.nsteps = nsteps;

	return 0;
}

int run_circuit(data_id fid, size_t nsteps)
{
	int rt = -1; /* Return value */

	struct timespec t1, t2;
	double t_tot;

	struct circ_trott_data cd;
	if (circ_trott_data_init_from_file(&cd, nsteps, fid) < 0)
		goto exit_trott_data;

	clock_gettime(CLOCK_REALTIME, &t1);
	if (circ_trott_simulate(&cd) < 0)
		goto exit_trott_simulate;

	clock_gettime(CLOCK_REALTIME, &t2);
	t_tot = (double)(t2.tv_sec - t1.tv_sec) +
		(double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	if (data_circ_trott_write_values(fid, cd.tsteps, nsteps) < 0)
		goto exit_trott_write;

	rt = 0; /* Success. */

exit_trott_write:
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%d,%.3f", cd.hamil.nqb, cd.hamil.nterms,
		cd.multidet.ndets, cd.ntsteps, WD.size, t_tot);

exit_trott_simulate:
	circ_trott_data_destroy(&cd);

exit_trott_data:

	return rt;
}

int main(int argc, char **argv)
{
	int rt = -1; /* Return value. */

	if (world_init(&argc, &argv, WD_SEED) != WORLD_READY)
		goto exit_world_init;
	world_info(&WD);

	unsigned int nranks = WD.size;
	if (nranks == 0 || (nranks & (nranks - 1)) != 0) {
		log_error("number of MPI ranks (%u) "
			  "must be a power of two.",
			nranks);
		goto exit_nranks;
	}

	log_info("*** Init ***");
	log_info("MPI nranks: %d", nranks);
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
	log_info("Num_steps: %zu", OPT.nsteps);
	if (run_circuit(fid, OPT.nsteps) < 0) {
		log_error("Failure: simulation error");
		goto exit_run_circuit;
	}

	rt = 0; /* Success. */

exit_run_circuit:
	log_info("Shut down simulation environment");
	data_close(fid);
exit_data_open:
exit_opt_parse:
exit_nranks:
exit_world_init:
	world_destroy();

	return rt;
}
