#define _XOPEN_SOURCE 600
#include "c23_compat.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "circ/cmpsit.h"
#include "ph2run.h"
#include "phase2.h"

#define WD_SEED UINT64_C(0xd871e5d39fc0222d)
static struct world WD;

static struct args {
	const char *progname;
	bool opt_help;
	bool opt_version;
	size_t depth;
	size_t length;
	size_t samples;
	size_t steps;
	double step_size;
	const char *filename;
} ARGS = {
	.progname = nullptr,
	.opt_help = false,
	.opt_version = false,
	.depth = 64,
	.length = 1,
	.samples = 1,
	.step_size = 1.0,
	.steps = 1,
	.filename = nullptr,
};
static void print_help(const char *progname)
{
	fprintf(stderr, "%s: Simulate \"composite\" circuit.\n\n", progname);
	fprintf(stderr, "  usage: %s [OPTIONS] FILENAME\n", progname);
	fprintf(stderr,
		"\nOptions:\n"
		"  -h, --help          Show this help.\n"
		"  -v, --version       Print version number.\n"
		"  --depth=64          Depth of the sampled sub-circuit.\n"
		"  --length=1          Length of the deterministic Hamiltonian.\n"
		"  --samples=1         Number of samples.\n"
		"  --step-size=1.0     Time evolution step size.\n"
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
	if (strncmp(o, "--depth=", 8) == 0) {
		size_t depth = strtoul(o + 8, nullptr, 10);
		if (depth == 0) {
			fprintf(stderr, "Option: --depth=, wrong value.\n");
			return -1;
		}
		ARGS.depth = depth;

		return 0;
	}
	if (strncmp(o, "--help", 6) == 0) {
		ARGS.opt_help = true;
		return 0;
	}
	if (strncmp(o, "--length=", 9) == 0) {
		size_t length = strtoul(o + 9, nullptr, 10);
		if (length == 0) {
			fprintf(stderr, "Option: --length=, wrong value.\n");
			return -1;
		}
		ARGS.length = length;

		return 0;
	}
	if (strncmp(o, "--samples=", 10) == 0) {
		size_t samples = strtoul(o + 10, nullptr, 10);
		if (samples == 0) {
			fprintf(stderr, "Option: --samples=, wrong value.\n");
			return -1;
		}
		ARGS.samples = samples;

		return 0;
	}
	if (strncmp(o, "--step-size=", 12) == 0) {
		double ss = strtod(o + 12, nullptr);
		if (ss == 0) {
			fprintf(stderr, "Option: --step-size=, wrong value.\n");
			return -1;
		}
		ARGS.step_size = ss;

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

int run_circuit(const struct args *args)
{
	int rt = -1; /* Return value */

	data_id fid;
	struct timespec t1, t2;

	struct cmpsit ct;
	struct cmpsit_data data = { /* */
		.depth = args->depth,
		.length = args->length,
		.samples = args->samples,
		.step_size = args->step_size,
		.steps = args->steps
	};

	log_info("open data file: %s", args->filename);
	if ((fid = data_open(args->filename)) == DATA_INVALID_FID) {
		log_error("open file: %s", args->filename);
		goto ex_circ_init;
	}
	if (cmpsit_init(&ct, &data, fid) < 0)
		goto ex_circ_init;
	log_info("close data file: %s", args->filename);
	data_close(fid);

	clock_gettime(CLOCK_REALTIME, &t1);
	if (circ_simulate(&ct.ct) < 0)
		goto ex_circ_simulate;
	clock_gettime(CLOCK_REALTIME, &t2);
	const double t_tot = (double)(t2.tv_sec - t1.tv_sec) +
			     (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	log_info("open data file: %s", args->filename);
	if ((fid = data_open(args->filename)) == DATA_INVALID_FID) {
		log_error("open file: %s", args->filename);
		goto ex_circ_res_write;
	}
	if (circ_write_res(&ct.ct, fid) < 0)
		goto ex_circ_res_write;
	log_info("close data file: %s", args->filename);
	data_close(fid);

	rt = 0; /* Success. */
ex_circ_res_write:
	cmpsit_destroy(&ct);
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_samples,depth,length,step_size,steps,"
		 "n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%zu,%zu,%.6f,%zu,%d,%.3f",
		ct.ct.hamil.nqb, ct.ct.hamil.nterms, ct.ct.muldet.ndets,
		data.samples, data.depth, data.length, data.step_size,
		data.steps, WD.size, t_tot);
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
	log_info("MPI num_ranks: %d", nranks);
	log_info("world_backend: %s", WORLD_BACKEND);
	log_info("*** Circuit: cmpsit >>> ***");
	log_info("depth: %zu", ARGS.depth);
	log_info("length: %zu", ARGS.length);
	log_info("samples: %zu", ARGS.samples);
	log_info("step_size: %f", ARGS.step_size);
	log_info("steps: %zu", ARGS.steps);

	if (run_circuit(&ARGS) < 0) {
		log_error("Failure: simulation error");
		goto exit_run_circuit;
	}

	rt = 0; /* Success. */

exit_run_circuit:
	log_info("Shut down simulation environment");
exit_nranks:
exit_world_init:
	world_destroy();

	return rt;
}
