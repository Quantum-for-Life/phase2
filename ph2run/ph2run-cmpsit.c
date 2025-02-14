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
	size_t length;
	size_t depth;
	size_t steps;
	double angle_det;
	double angle_rand;
	size_t samples;
	uint64_t seed;
	const char *filename;
} ARGS = {
	.progname = nullptr,
	.opt_help = false,
	.opt_version = false,
	.seed = 0,
	.length = 1,
	.depth = 64,
	.steps = 1,
	.angle_det = 1.0,
	.angle_rand = 1.0,
	.samples = 1,
	.filename = nullptr,
};

static void print_help(const char *progname)
{
	fprintf(stderr, "%s-%d.%d.%d: Simulate \"trott\" circuit.\n\n",
		progname, PHASE2_VER_MAJOR, PHASE2_VER_MINOR, PHASE2_VER_PATCH);
	fprintf(stderr, "  usage: %s [OPTIONS] FILENAME\n", progname);
	fprintf(stderr,
		"\nOptions:\n"
		"  -h, --help          Show this help.\n"
		"  -v, --version       Print version number.\n"
		"  --seed=N            Seed of the PRNG (default value if not specified).\n"
		"  --length=1          Length of the deterministic Hamiltonian.\n"
		"  --depth=64          Depth of the sampled circuit.\n"
		"  --steps=1           Number of Trotter steps.\n"
		"  --angle-det=1.0     Deterministic time evolution angle (delta).\n"
		"  --angle-rand=1.0    Randomized evolution angle (tau).\n"
		"  --samples=1         Number of samples.\n"
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
	if (strncmp(o, "--help", 6) == 0) {
		ARGS.opt_help = true;
		return 0;
	}
	if (strncmp(o, "--depth=", 8) == 0) {
		size_t depth = strtoul(o + 8, nullptr, 10);
		if (depth == 0) {
			fprintf(stderr, "Option: --depth=, wrong value.\n");
			return -1;
		}
		ARGS.depth = depth;

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
	if (strncmp(o, "--angle-det=", 12) == 0) {
		double ss = strtod(o + 12, nullptr);
		if (ss == 0) {
			fprintf(stderr, "Option: --angle-det=, wrong value.\n");
			return -1;
		}
		ARGS.angle_det = ss;

		return 0;
	}
	if (strncmp(o, "--angle-rand=", 13) == 0) {
		double tau = strtod(o + 13, nullptr);
		if (tau == 0) {
			fprintf(stderr, "Option: --anglerand=, wrong value.\n");
			return -1;
		}
		ARGS.angle_rand = tau;

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
	if (strncmp(o, "--seed=", 7) == 0) {
		uint64_t seed = strtoul(o + 7, nullptr, 10);
		if (seed != 0) {
			ARGS.seed = seed;
		}

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
	struct cmpsit_data data = {
		/* */
		.depth = args->depth,
		.length = args->length,
		.samples = args->samples,
		.angle_det = args->angle_det,
		.steps = args->steps,
		.seed = args->seed,
		.angle_rand = args->angle_rand,
	};

	struct cmpsit cp;

	log_info("open data file: %s", args->filename);
	if ((fid = data_open(args->filename)) == DATA_INVALID_FID) {
		log_error("open file: %s", args->filename);
		goto ex_circ_init;
	}
	if (cmpsit_init(&cp, &data, fid) < 0)
		goto ex_circ_init;
	log_info("close data file: %s", args->filename);
	data_close(fid);

	clock_gettime(CLOCK_REALTIME, &t1);
	if (cmpsit_simul(&cp) < 0)
		goto ex_circ_simulate;
	clock_gettime(CLOCK_REALTIME, &t2);
	const double t_tot = (double)(t2.tv_sec - t1.tv_sec) +
			     (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	log_info("open data file: %s", args->filename);
	if ((fid = data_open(args->filename)) == DATA_INVALID_FID) {
		log_error("open file: %s", args->filename);
		goto ex_circ_res_write;
	}
	if (cmpsit_write_res(&cp, fid) < 0)
		goto ex_circ_res_write;
	log_info("close data file: %s", args->filename);
	data_close(fid);

	rt = 0; /* Success. */
ex_circ_res_write:
	cmpsit_free(&cp);
	log_info("> Simulation summary (CSV):");
	log_info(
		"> n_qb,n_terms,n_dets,samples,length,depth,angle_det,angle_rand,steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%zu,%zu,%.16f,%.16f,%zu,%d,%.3f",
		cp.ct.hm.qb, cp.ct.hm.len, cp.ct.md.len, data.samples,
		data.length, data.depth, data.angle_det, data.angle_rand,
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
	log_info("seed: %lu", ARGS.seed);
	log_info("length: %zu", ARGS.length);
	log_info("depth: %zu", ARGS.depth);
	log_info("steps: %zu", ARGS.steps);
	log_info("angle_det: %.16f", ARGS.angle_det);
	log_info("angle_rand: %.16f", ARGS.angle_rand);
	log_info("samples: %zu", ARGS.samples);

	if (run_circuit(&ARGS) < 0) {
		log_error("Failure: simulation error");
		goto exit_run_circuit;
	}

	rt = 0; /* Success. */

exit_run_circuit:
	log_info("Shut down simulation environment");
exit_nranks:
exit_world_init:
	world_free();

	return rt;
}
