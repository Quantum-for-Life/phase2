#define _XOPEN_SOURCE 600
#include "c23_compat.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "circ/trott.h"
#include "phase2.h"

#define WD_SEED UINT64_C(0xd326119d4859ebb2)
static struct world WD;

static struct args {
	const char *progname;
	bool opt_help;
	bool opt_version;
	double delta;
	size_t steps;
	const char *filename;
} ARGS = {
	.progname = nullptr,
	.opt_help = false,
	.opt_version = false,
	.delta = 1.0,
	.steps = 1,
	.filename = nullptr,
};

static void print_help(const char *progname)
{
	fprintf(stderr, "%s: Simulate \"trott\" circuit.\n\n", progname);
	fprintf(stderr, "  usage: %s [OPTIONS] FILENAME\n", progname);
	fprintf(stderr, "\nOptions:\n"
			"  -h, --help          Show this help.\n"
			"  -v, --version       Print version number.\n"
			"  --delta=1.0         Time scale.\n"
			"  --steps=1           Number of Trotter steps.\n"
			"\n");
	fprintf(stderr, "FILENAME is a HDF5 simulation worksheet.\n");
}

static void print_version(const char *progname)
{
	fprintf(stderr, "%s-v%d.%d.%d\n", progname, PHASE2_VER_MAJOR,
		PHASE2_VER_MINOR, PHASE2_VER_PATCH);
}

static int args_parse_shortopt(const int *argc, char ***argv)
{
	(void)argc;

	char *o = **argv + 1;
	while (*o != '\0') {
		switch (*o) {
		case ('h'):
			ARGS.opt_help = true;
			break;
		case ('v'):
			ARGS.opt_version = true;
			break;
		default:
			fprintf(stderr, "Unrecognized option: -%c\n", *o);
			return -1;
		}
		o++;
	}

	return 0;
}

static int args_parse_longopt(const int *argc, char ***argv)
{
	(void)argc;

	char *o = **argv;
	if (strncmp(o, "--delta=", 8) == 0) {
		double delta = strtod(o + 8, nullptr);
		if (delta == 0) {
			fprintf(stderr, "Option: --delta=, wrong value.\n");
			return -1;
		}
		ARGS.delta = delta;

		return 0;
	}
	if (strncmp(o, "--help", 6) == 0) {
		ARGS.opt_help = true;
		return 0;
	}
	if (strncmp(o, "--steps=", 8) == 0) {
		unsigned long n = strtoull(o + 8, nullptr, 10);
		if (n == 0) {
			fprintf(stderr,
				"Option: --steps=n, wrong no. of steps.\n");
			return -1;
		}
		ARGS.steps = n;

		return 0;
	}
	if (strncmp(o, "--version", 9) == 0) {
		ARGS.opt_version = true;
		return 0;
	}

	fprintf(stderr, "Unrecognized option: %s (no option argument?)\n", o);
	return -1;
}

static int args_parse(int argc, char **argv)
{
	const char *prog;
	ARGS.progname = prog = *argv;
	while (*prog != '\0') {
		if (*prog == '/')
			ARGS.progname = prog + 1;
		prog++;
	}

	while (++argv, --argc) {
		switch (argv[0][0]) {
		case ('-'):
			if (argv[0][1] != '-' &&
				args_parse_shortopt(&argc, &argv) < 0)
				return -1;
			if (argv[0][1] == '-' &&
				args_parse_longopt(&argc, &argv) < 0)
				return -1;
			break;
		default:
			if (!ARGS.filename)
				ARGS.filename = argv[0];
			else {
				fprintf(stderr, "Unrecognized argument: %s\n",
					argv[0]);
				return -1;
			}
		}
	}

	return 0;
}

static void args_validate(void)
{
	if (ARGS.opt_help) {
		print_help(ARGS.progname);
		exit(0);
	}
	if (ARGS.opt_version) {
		print_version(ARGS.progname);
		exit(0);
	}
	if (!ARGS.filename) {
		fprintf(stderr,
			"No simulation file specified.  "
			"See '%s -h' for more detail.\n",
			ARGS.progname);
		exit(-1);
	}
}

static int run_circuit(const struct args *args)
{
	int rt = -1; /* Return value */

	data_id fid;
	struct timespec t1, t2;

	struct trott tt;
	struct trott_data tt_dt;
	tt_dt.delta = args->delta;
	tt_dt.steps = args->steps;

	log_info("open data file: %s", args->filename);
	if ((fid = data_open(args->filename)) == DATA_INVALID_FID) {
		log_error("open file: %s", args->filename);
		goto ex_circ_init;
	}
	if (trott_init(&tt, &tt_dt, fid) < 0)
		goto ex_circ_init;
	log_info("close data file: %s", args->filename);
	data_close(fid);

	clock_gettime(CLOCK_REALTIME, &t1);
	if (circ_simul(&tt.ct) < 0)
		goto ex_circ_simulate;
	clock_gettime(CLOCK_REALTIME, &t2);
	const double t_tot = (double)(t2.tv_sec - t1.tv_sec) +
			     (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	log_info("open data file: %s", args->filename);
	if ((fid = data_open(args->filename)) == DATA_INVALID_FID) {
		log_error("open file: %s", args->filename);
		goto ex_circ_write_res;
	}
	if (trott_write_res(&tt, fid) < 0)
		goto ex_circ_write_res;
	log_info("close data file: %s", args->filename);
	data_close(fid);

	rt = 0; /* Success. */
ex_circ_write_res:
	trott_destroy(&tt);
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,delta,n_steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%f,%zu,%d,%.3f", tt.ct.hamil.qb,
		tt.ct.hamil.len, tt.ct.muldet.len, args->delta,
		args->steps, WD.size, t_tot);
ex_circ_simulate:
ex_circ_init:
	return rt;
}

int main(int argc, char **argv)
{
	int rt = -1; /* Return value. */
	if (args_parse(argc, argv) < 0)
		return -1;
	args_validate();

	if (world_init(&argc, &argv, WD_SEED) != WORLD_READY)
		goto ex_world_init;
	world_info(&WD);

	const unsigned int nranks = WD.size;
	if (nranks == 0 || (nranks & (nranks - 1)) != 0) {
		log_error("number of MPI ranks (%u) "
			  "must be a power of two.",
			nranks);
		goto ex_nranks;
	}

	log_info("*** Init ***");
	log_info("MPI nranks: %d", nranks);
	log_info("world_backend: %s", WORLD_BACKEND);

	log_info("*** Circuit: trott ***");
	log_info("delta: %f", ARGS.delta);
	log_info("num_steps: %zu", ARGS.steps);
	if (run_circuit(&ARGS) < 0) {
		log_error("Failure: simulation error");
		goto ex_run_circuit;
	}

	rt = 0; /* Success. */

ex_run_circuit:
	log_info("Shut down simulation environment");
ex_nranks:
ex_world_init:
	world_destroy();

	return rt;
}
