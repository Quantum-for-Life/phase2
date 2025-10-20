#define _XOPEN_SOURCE 600
#include "c23_compat.h"

#include <argp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "circ/cmpsit.h"
#include "circ/qdrift.h"
#include "circ/trott.h"
#include "phase2.h"

#define WD_SEED UINT64_C(0xd326119d4859ebb2)
static struct world wd;

#define xstr(s) str(s)
#define str(s) #s

const char *argp_program_version = PHASE2_VERSION;
const char *argp_program_bug_address = "Marek Miller <mlm@math.ku.dk>";

#define doc "Run phase2 sumilations."
#define args_doc "CMD [CMD_OPT]"

static struct args {
	bool verbose;
	char *simul;
	char *cmd;
	unsigned cmd_num;
} args = { .verbose = false,
	.simul = "./simul.h5",
	.cmd = nullptr,
	.cmd_num = 0 };

static struct argp_option opts[] = {
	{ "verbose", 'v', 0, 0, "Print verbose output", 0 },
	{ "simul", 'S', "FILE", 0,
		"Simulation HDF5 data file (default: ./simul.h5)", 0 },
	{ 0 }
};

static error_t opts_parser(int key, char *arg, struct argp_state *state)
{
	struct args *args = state->input;

	switch (key) {
	case 'v':
		args->verbose = true;
		break;
	case 'S':
		args->simul = arg;
		break;

	case ARGP_KEY_ARG:
		args->cmd = arg;
		args->cmd_num = state->next - 1;
		state->next = state->argc;
		break;

	case ARGP_KEY_NO_ARGS:
		argp_usage(state);
		break;

	default:
		return ARGP_ERR_UNKNOWN;
	}

	return 0;
}

static struct argp argp = { opts, opts_parser, args_doc, doc, 0, 0, 0 };

/* Commands */
enum {
	CMD_TROTT = 0,
	CMD_QDRIFT = 1,
	CMD_CMPSIT = 2,
};

#define CMD_TROTT_STR "trott"
#define CMD_QDRIFT_STR "qdrift"
#define CMD_CMPSIT_STR "cmpsit"

/* Command: "trott" */
#define doc_trott "Run trotter algo"
#define argv0_trott "ph2run [OPTS] trott"
#define args_doc_trott ""

static struct args_trott {
	double delta;
	size_t steps;
} args_trott = { .delta = 1.0, .steps = 1 };

static struct argp_option opts_trott[] = {
	{ "delta", 'D', "VAL", 0, "Floating point number (default: 1.0)", 0 },
	{ "steps", 's', "N", 0, "Number of Trotter steps", 0 },
	{ 0 }
};

static error_t opts_parser_trott(int key, char *arg, struct argp_state *state)
{
	struct args_trott *args = state->input;

	switch (key) {
	case 'D':
		args->delta = strtod(arg, nullptr);
		break;
	case 's':
		args->steps = strtoull(arg, nullptr, 10);
		break;

	case ARGP_KEY_ARG:
		argp_usage(state);
		break;

	case ARGP_KEY_NO_ARGS:
		break;

	default:
		return ARGP_ERR_UNKNOWN;
	}

	return 0;
}

static struct argp argp_trott = { opts_trott, opts_parser_trott, args_doc_trott,
	doc_trott, 0, 0, 0 };

int cmd_trott(void)
{
	int rt = -1;

	log_info("*** Circuit: trott ***");
	log_info("delta: %f", args_trott.delta);
	log_info("num_steps: %zu", args_trott.steps);

	data_id fid;
	struct timespec t1, t2;

	struct trott tt;
	struct trott_data tt_dt;
	tt_dt.delta = args_trott.delta;
	tt_dt.steps = args_trott.steps;

	log_info("open data file: %s", args.simul);
	if ((fid = data_open(args.simul)) == DATA_INVALID_FID) {
		log_error("open file: %s", args.simul);
		goto ex_circ_init;
	}
	if (trott_init(&tt, &tt_dt, fid) < 0)
		goto ex_circ_init;
	log_info("close data file: %s", args.simul);
	data_close(fid);

	clock_gettime(CLOCK_REALTIME, &t1);
	if (trott_simul(&tt) < 0)
		goto ex_circ_simulate;
	clock_gettime(CLOCK_REALTIME, &t2);
	const double t_tot = (double)(t2.tv_sec - t1.tv_sec) +
			     (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	log_info("open data file: %s", args.simul);
	if ((fid = data_open(args.simul)) == DATA_INVALID_FID) {
		log_error("open file: %s", args.simul);
		goto ex_circ_write_res;
	}
	if (trott_write_res(&tt, fid) < 0)
		goto ex_circ_write_res;
	log_info("close data file: %s", args.simul);
	data_close(fid);

	rt = 0; /* Success. */

ex_circ_write_res:
	trott_free(&tt);
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,delta,n_steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%f,%zu,%d,%.3f", tt.ct.hm.qb, tt.ct.hm.len,
		tt.ct.md.len, args_trott.delta, args_trott.steps, wd.size,
		t_tot);
ex_circ_simulate:
ex_circ_init:
	log_info("Shut down simulation environment");

	return rt;
}

/* Command: "qdrift" */
#define doc_qdrift "Run qDrift algorithm"
#define argv0_qdrift "ph2run [OPTS] qdrift"
#define args_doc_qdrift ""

static struct args_qdrift {
	double delta;
	size_t depth;
	size_t nsamples;
	double step_size;
	uint64_t seed;
} args_qdrift = {
	.delta = 1.0,
	.depth = 64,
	.nsamples = 1,
	.seed = 0,
};

static struct argp_option opts_qdrift[] = {
	{ "delta", 'D', "VAL", 0, "Floating point number (default: 1.0)", 0 },
	{ "depth", 'd', "VAL", 0, "Random depth (default: 64)", 0 },
	{ "samples", 'n', "N", 0, "Number of independent samples (default: 1)",
		0 },
	{ "seed", 'x', "", 0, "Seed to the pseudo random number generator", 0 },
	{ 0 }
};

static error_t opts_parser_qdrift(int key, char *arg, struct argp_state *state)
{
	struct args_qdrift *args = state->input;

	switch (key) {
	case 'D':
		args->delta = strtod(arg, nullptr);
		break;
	case 'd':
		args->depth = strtoull(arg, nullptr, 10);
		break;
	case 'n':
		args->nsamples = strtoull(arg, nullptr, 10);
		break;
	case 's':
		args->step_size = strtod(arg, nullptr);
		break;
	case 'x':
		args->seed = strtoull(arg, nullptr, 10);
		break;

	case ARGP_KEY_ARG:
		argp_usage(state);
		break;

	case ARGP_KEY_NO_ARGS:
		break;

	default:
		return ARGP_ERR_UNKNOWN;
	}

	return 0;
}

static struct argp argp_qdrift = { opts_qdrift, opts_parser_qdrift,
	args_doc_qdrift, doc_qdrift, 0, 0, 0 };

int cmd_qdrift(void)
{
	int rt = -1;

	log_info("*** Circuit: qDRIFT >>> ***");
	log_info("delta: %f", args_qdrift.delta);
	log_info("depth: %zu", args_qdrift.depth);
	log_info("samples: %zu", args_qdrift.nsamples);
	log_info("seed: %lu", args_qdrift.seed);

	data_id fid;
	struct timespec t1, t2;

	struct qdrift qd;
	struct qdrift_data data = { .depth = args_qdrift.depth,
		.samples = args_qdrift.nsamples,
		.step_size = args_qdrift.delta,
		.seed = args_qdrift.seed };

	log_info("open data file: %s", args.simul);
	if ((fid = data_open(args.simul)) == DATA_INVALID_FID) {
		log_error("open file: %s", args.simul);
		goto ex_circ_init;
	}
	if (qdrift_init(&qd, &data, fid) < 0)
		goto ex_circ_init;
	log_info("close data file: %s", args.simul);
	data_close(fid);

	clock_gettime(CLOCK_REALTIME, &t1);
	if (qdrift_simul(&qd) < 0)
		goto ex_circ_simulate;
	clock_gettime(CLOCK_REALTIME, &t2);
	const double t_tot = (double)(t2.tv_sec - t1.tv_sec) +
			     (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	log_info("open data file: %s", args.simul);
	if ((fid = data_open(args.simul)) == DATA_INVALID_FID) {
		log_error("open file: %s", args.simul);
		goto ex_circ_res_write;
	}
	if (qdrift_write_res(&qd, fid) < 0)
		goto ex_circ_res_write;
	log_info("close data file: %s", args.simul);
	data_close(fid);

	rt = 0; /* Success. */

ex_circ_res_write:
	qdrift_free(&qd);
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_samples,step_size,depth,"
		 "n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%zu,%.3f,%d,%.3f", qd.ct.hm.qb,
		qd.ct.hm.len, qd.ct.md.len, data.samples, data.step_size,
		data.depth, wd.size, t_tot);
ex_circ_simulate:
ex_circ_init:
	log_info("Shut down simulation environment");

	return rt;
}

/* Command: "cmpsit" */

int main(int argc, char **argv)
{
	int rt = -1;

	argp_parse(&argp, argc, argv, ARGP_IN_ORDER, nullptr, &args);

	/* Parse subcommands. */
	argc -= args.cmd_num;
	argv += args.cmd_num;
	int cmd = -1;

#define cmd_parse(name, str, val)                                              \
	({                                                                     \
		if (strncmp(args.cmd, str, strlen(str)) == 0) {                \
			argv[0] = argv0_##name;                                \
			argp_parse(&argp_##name, argc, argv, ARGP_IN_ORDER,    \
				nullptr, &args_##name);                        \
			cmd = val;                                             \
		}                                                              \
	})

	cmd_parse(trott, CMD_TROTT_STR, CMD_TROTT);
	cmd_parse(qdrift, CMD_QDRIFT_STR, CMD_QDRIFT);

	if (world_init(nullptr, nullptr, WD_SEED) != WORLD_READY)
		exit(-1);

	switch (cmd) {
	case CMD_TROTT:
		rt = cmd_trott();
		break;
	case CMD_QDRIFT:
		rt = cmd_qdrift();
		break;
	default:
		fprintf(stderr, "Unrecognized command.\n");
	}

	world_free();
	exit(rt);
}
