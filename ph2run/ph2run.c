#define _XOPEN_SOURCE 600
#include "c23_compat.h"

#include <argp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "circ/trott.h"
#include "phase2.h"

#define WD_SEED UINT64_C(0xd326119d4859ebb2)
static struct world WD;

#define xstr(s) str(s)
#define str(s) #s

const char *argp_program_version = PHASE2_VERSION;
const char *argp_program_bug_address = "Marek Miller <mlm@math.ku.dk>";

#define doc "Run phase2 sumilations."
#define args_doc "CMD [CMD_OPT]"

static struct args {
	bool verbose;
	char *simul;
	double delta;
	char *cmd;
	unsigned cmd_num;
} args = {
	.verbose	= false,
	.simul		= "./simul.h5",
	.delta		= 1.0,
	.cmd		= nullptr,
	.cmd_num	= 0
};

static struct argp_option opts[] = {
	{ "verbose", 'v', 0, 0,
		"Print verbose output", 0 },
	{ "simul", 'S', "FILE", 0,
		"Simulation HDF5 data file (default: ./simul.h5)", 0 }, 
	{ "delta", 'D', "VALUE", 0,
		"Floating point number (default: 1.0)", 0},
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
	case 'D':
		args->delta = strtod(arg, nullptr);
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

#define CMD_TROTT "trott"
#define CMD_QDRIFT "qdrift"
#define CMD_CMPSIT "cmpsit"

#define doc_trott "Run trotter algo"
#define argv0_trott "ph2run [OPTS] trott"
#define args_doc_trott ""

static struct args_trott {
	size_t steps;
} args_trott = {
	.steps = 1
};

static struct argp_option opts_trott[] = {
	{ "steps", 's', "N", 0, 
		"Number of Trotter steps", 0},
	{ 0 }
};

static error_t opts_parser_trott(int key, char *arg, struct argp_state *state)
{
	struct args_trott *args = state->input;

	switch (key) {
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

static struct argp argp_trott = { opts_trott, opts_parser_trott, 
	args_doc_trott, doc_trott, 0, 0, 0 };

static int trott_run_circuit(void)
{
	int rt = -1; /* Return value */

	data_id fid;
	struct timespec t1, t2;

	struct trott tt;
	struct trott_data tt_dt;
	tt_dt.delta = args.delta;
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
		tt.ct.md.len, args.delta, args_trott.steps, WD.size, t_tot);
ex_circ_simulate:
ex_circ_init:
	return rt;
}

int cmd_trott(void)
{
	int rt = -1;

	if (world_init(nullptr, nullptr, WD_SEED) != WORLD_READY)
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
	log_info("delta: %f", args.delta);
	log_info("num_steps: %zu", args_trott.steps);
	if (trott_run_circuit() < 0) {
		log_error("Failure: simulation error");
		goto ex_run_circuit;
	}

	rt = 0; /* Success. */

ex_run_circuit:
	log_info("Shut down simulation environment");
ex_nranks:
ex_world_init:
	world_free();

	return rt;
}

int main(int argc, char **argv)
{
	argp_parse(&argp, argc, argv, ARGP_IN_ORDER, nullptr, &args);

	/* Parse subcommands. */
	argc -= args.cmd_num;
	argv += args.cmd_num;
	if (strncmp(args.cmd, CMD_TROTT, strlen(CMD_TROTT)) == 0) {
		argv[0] = argv0_trott;
		argp_parse(&argp_trott, argc, argv, ARGP_IN_ORDER, nullptr,
			&args_trott);
		return cmd_trott();
	}


	exit(0);
}
