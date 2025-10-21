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

#define str(s) #s
#define xstr(s) str(s)

int timeit(int (*fn)(void *), void *data, double *t)
{
	int rt = -1;
	struct timespec t1, t2;

	clock_gettime(CLOCK_REALTIME, &t1);
	rt = fn(data);
	clock_gettime(CLOCK_REALTIME, &t2);

	if (t)
		*t = (double)(t2.tv_sec - t1.tv_sec) +
		     (double)(t2.tv_nsec - t1.tv_nsec) * 1.0e-9;

	return rt;
}

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

static int data_exec(int (*fn)(data_id, void *), void *data)
{
	int rt = -1;

	data_id fid;
	log_info("open data file: %s", args.simul);
	if ((fid = data_open(args.simul)) == DATA_INVALID_FID) {
		log_error("open file: %s", args.simul);
		return -1;
	}

	rt = fn(fid, data);

	log_info("close data file: %s", args.simul);
	data_close(fid);

	return rt;
}

/* Commands */
enum {
	CMD_TROTT = 0,
	CMD_QDRIFT = 1,
	CMD_CMPSIT = 2,
};

/* Command: "trott" */
static struct cmd_trott_dt {
	struct trott tt;
	struct trott_data tt_dt;
} cmd_trott_dt = {
	.tt_dt = { .delta = 1.0, .steps = 1 }
};

static int cmd_trott_init(data_id fid, void *data)
{
	struct cmd_trott_dt *const dt = data;

	log_info("*** Circuit: trott ***");
	log_info("delta: %f", dt->tt_dt.delta);
	log_info("num_steps: %zu", dt->tt_dt.steps);

	return trott_init(&dt->tt, &dt->tt_dt, fid);
}

static int cmd_trott_write(data_id fid, void *data)
{
	int rt = -1;
	struct cmd_trott_dt *const dt = data;

	rt = trott_write_res(&dt->tt, fid);
	trott_free(&dt->tt);

	return rt;
}

static int cmd_trott_run(void *data)
{
	struct cmd_trott_dt *dt = data;

	return trott_simul(&dt->tt);
}

#define doc_trott "Run trotter algo"
#define argv0_trott "ph2run [OPTS] trott"
#define args_doc_trott ""

/* Alias for the option parser. */
static struct trott_data *const args_trott = &cmd_trott_dt.tt_dt;

static struct argp_option opts_trott[] = {
	{ "delta", 'D', "VAL", 0, "Floating point number (default: 1.0)", 0 },
	{ "steps", 's', "N", 0, "Number of Trotter steps", 0 },
	{ 0 }
};

static error_t opts_parser_trott(int key, char *arg, struct argp_state *state)
{
	struct trott_data *dt = state->input;

	switch (key) {
	case 'D':
		dt->delta = strtod(arg, nullptr);
		break;
	case 's':
		dt->steps = strtoull(arg, nullptr, 10);
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

static int cmd_trott(void)
{
	int rt = -1;
	double t_tot;

	if (data_exec(cmd_trott_init, &cmd_trott_dt) < 0)
		goto ex;
	if (timeit(cmd_trott_run, &cmd_trott_dt, &t_tot) < 0) {
		log_error("Simulation error");
		goto ex;
	}
	if (data_exec(cmd_trott_write, &cmd_trott_dt) < 0)
		goto ex;

	rt = 0; /* Success. */

	struct trott *const tt = &cmd_trott_dt.tt;
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,delta,n_steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%f,%zu,%d,%.3f", tt->ct.hm.qb, tt->ct.hm.len,
		tt->ct.md.len, args_trott->delta, args_trott->steps, wd.size,
		t_tot);

ex:
	log_info("Shut down simulation environment");

	return rt;
}

/* Command: "qdrift" */
static struct cmd_qdrift_dt {
	struct qdrift qd;
	struct qdrift_data qd_dt;
} cmd_qdrift_dt = {
	.qd_dt = { .step_size = 1.0, .depth = 64, .samples = 1, .seed = 23 }
};

static int cmd_qdrift_init(data_id fid, void *data)
{
	struct cmd_qdrift_dt *const dt = data;

	log_info("*** Circuit: qDRIFT >>> ***");
	log_info("delta: %f", dt->qd_dt.step_size);
	log_info("depth: %zu", dt->qd_dt.depth);
	log_info("samples: %zu", dt->qd_dt.samples);
	log_info("seed: %lu", dt->qd_dt.seed);

	return qdrift_init(&dt->qd, &dt->qd_dt, fid);
}

static int cmd_qdrift_write(data_id fid, void *data)
{
	int rt = -1;
	struct cmd_qdrift_dt *const dt = data;

	rt = qdrift_write_res(&dt->qd, fid);
	qdrift_free(&dt->qd);

	return rt;
}

static int cmd_qdrift_run(void *data)
{
	struct cmd_qdrift_dt *dt = data;

	return qdrift_simul(&dt->qd);
}

#define doc_qdrift "Run qDrift algorithm"
#define argv0_qdrift "ph2run [OPTS] qdrift"
#define args_doc_qdrift ""

static struct qdrift_data *const args_qdrift = &cmd_qdrift_dt.qd_dt;

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
	struct qdrift_data *dt = state->input;

	switch (key) {
	case 'D':
		dt->step_size = strtod(arg, nullptr);
		break;
	case 'd':
		dt->depth = strtoull(arg, nullptr, 10);
		break;
	case 'n':
		dt->samples = strtoull(arg, nullptr, 10);
		break;
	case 'x':
		dt->seed = strtoull(arg, nullptr, 10);
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
	double t_tot;

	if (data_exec(cmd_qdrift_init, &cmd_qdrift_dt) < 0)
		goto ex;
	if (timeit(cmd_qdrift_run, &cmd_qdrift_dt, &t_tot) < 0) {
		log_error("Simulation error");
		goto ex;
	}
	if (data_exec(cmd_qdrift_write, &cmd_qdrift_dt) < 0)
		goto ex;

	rt = 0; /* Success. */

	const struct qdrift *qd = &cmd_qdrift_dt.qd;
	const struct qdrift_data *qd_dt = &cmd_qdrift_dt.qd_dt;
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_samples,step_size,depth,"
		 "n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%zu,%.3f,%d,%.3f", qd->ct.hm.qb,
		qd->ct.hm.len, qd->ct.md.len, qd_dt->samples, qd_dt->step_size,
		qd_dt->depth, wd.size, t_tot);
ex:
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

#define cmd_parse(name, val)                                                   \
	({                                                                     \
		if (strncmp(args.cmd, xstr(name), strlen(xstr(name))) == 0) {  \
			argv[0] = argv0_##name;                                \
			argp_parse(&argp_##name, argc, argv, ARGP_IN_ORDER,    \
				nullptr, args_##name);                         \
			cmd = val;                                             \
			break;                                                 \
		}                                                              \
	})

	while (1) {
		cmd_parse(trott, CMD_TROTT);
		cmd_parse(qdrift, CMD_QDRIFT);
	}

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
