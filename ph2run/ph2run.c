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
#include "log.h"
#include "phase2.h"

#define WD_SEED UINT64_C(0xd326119d4859ebb2)
static struct world_info wd;

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

#define doc                                                                    \
	"Run phase2 sumilations.  CMD can be one of the algorithms:\n"         \
	"  trott\n"                                                            \
	"  qdrift\n"                                                           \
	"  cmpsit\n"                                                           \
	"\nRun CMD --help for more information"
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
	double t_tot;
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

	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,delta,n_steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%f,%zu,%d,%.3f", dt->tt.ct.hm.qb,
		dt->tt.ct.hm.len, dt->tt.ct.md.len, dt->tt_dt.delta,
		dt->tt_dt.steps, wd.size, dt->t_tot);

	return rt;
}

static int cmd_trott_run(void *data)
{
	struct cmd_trott_dt *dt = data;

	return trott_simul(&dt->tt);
}

#define doc_trott "Run deterministic Trotter product formula (1st order)."
#define argv0_trott "ph2run [OPTS] trott"
#define args_doc_trott ""

/* Alias for the option parser. */
static struct trott_data *const args_trott = &cmd_trott_dt.tt_dt;

static struct argp_option opts_trott[] = {
	{ "delta", 'D', "VAL", 0,
		"Rotation angle, a floating point number (default: 1.0)", 0 },
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

	if (data_exec(cmd_trott_init, &cmd_trott_dt) < 0)
		goto ex;
	if (timeit(cmd_trott_run, &cmd_trott_dt, &cmd_trott_dt.t_tot) < 0) {
		log_error("Simulation error");
		goto ex;
	}
	if (data_exec(cmd_trott_write, &cmd_trott_dt) < 0)
		goto ex;

	rt = 0; /* Success. */
ex:
	log_info("Shut down simulation environment");

	return rt;
}

/* Command: "qdrift" */
static struct cmd_qdrift_dt {
	struct qdrift qd;
	struct qdrift_data qd_dt;
	double t_tot;
} cmd_qdrift_dt = {
	.qd_dt = { .step_size = 1.0, .depth = 64, .samples = 1, .seed = 1 }
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

	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_samples,step_size,depth,"
		 "n_ranks,t_tot");
	log_info("> %u,%zu,%zu,%zu,%.6f,%zu,%d,%.3f", dt->qd.ct.hm.qb,
		dt->qd.ct.hm.len, dt->qd.ct.md.len, dt->qd_dt.samples,
		dt->qd_dt.step_size, dt->qd_dt.depth, wd.size, dt->t_tot);

	return rt;
}

static int cmd_qdrift_run(void *data)
{
	struct cmd_qdrift_dt *dt = data;

	return qdrift_simul(&dt->qd);
}

#define doc_qdrift "Run qDRIFT randomized algorithm."
#define argv0_qdrift "ph2run [OPTS] qdrift"
#define args_doc_qdrift ""

static struct qdrift_data *const args_qdrift = &cmd_qdrift_dt.qd_dt;

static struct argp_option opts_qdrift[] = {
	{ "delta", 'D', "VAL", 0,
		"Rotation angle, a floating point number (default: 1.0)", 0 },
	{ "depth", 'd', "VAL", 0, "Random depth (default: 64)", 0 },
	{ "samples", 'n', "N", 0, "Number of independent samples (default: 1)",
		0 },
	{ "seed", 'x', "N", 0,
		"Seed of the pseudo random number generator (default: 1, must differ from zero)",
		0 },
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

	if (data_exec(cmd_qdrift_init, &cmd_qdrift_dt) < 0)
		goto ex;
	if (timeit(cmd_qdrift_run, &cmd_qdrift_dt, &cmd_qdrift_dt.t_tot) < 0) {
		log_error("Simulation error");
		goto ex;
	}
	if (data_exec(cmd_qdrift_write, &cmd_qdrift_dt) < 0)
		goto ex;

	rt = 0; /* Success. */
ex:
	log_info("Shut down simulation environment");

	return rt;
}

/* Command: "cmpsit" */
static struct cmd_cmpsit_dt {
	struct cmpsit cp;
	struct cmpsit_data cp_dt;
	double t_tot;
} cmd_cmpsit_dt = { .cp_dt = {
	/* */
	.seed = 1,
	.length = 1,
	.depth = 64,
	.steps = 1,
	.angle_det = 1.0,
	.angle_rand = 1.0,
	.samples = 1,
	}
};

static int cmd_cmpsit_init(data_id fid, void *data)
{
	struct cmd_cmpsit_dt *const dt = data;

	log_info("*** Circuit: cmpsit ***");
	log_info("seed: %lu", dt->cp_dt.seed);
	log_info("length: %zu", dt->cp_dt.length);
	log_info("depth: %zu", dt->cp_dt.depth);
	log_info("steps: %zu", dt->cp_dt.steps);
	log_info("angle_det: %.16f", dt->cp_dt.angle_det);
	log_info("angle_rand: %.16f", dt->cp_dt.angle_rand);
	log_info("samples: %zu", dt->cp_dt.samples);

	return cmpsit_init(&dt->cp, &dt->cp_dt, fid);
}

static int cmd_cmpsit_write(data_id fid, void *data)
{
	int rt = -1;
	struct cmd_cmpsit_dt *const dt = data;

	rt = cmpsit_write_res(&dt->cp, fid);
	cmpsit_free(&dt->cp);

	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,samples,length,depth,angle_det"
		 ",angle_rand,steps,n_ranks,t_tot");
	log_info("> %zu,%zu,%zu,%zu,%zu,%zu,%.16f,%.16f,%zu,%d,%.3f",
		dt->cp.ct.hm.qb, dt->cp.ct.hm.len, dt->cp.ct.md.len,
		dt->cp_dt.samples, dt->cp_dt.length, dt->cp_dt.depth,
		dt->cp_dt.angle_det, dt->cp_dt.angle_rand, dt->cp_dt.steps,
		wd.size, dt->t_tot);

	return rt;
}

static int cmd_cmpsit_run(void *data)
{
	struct cmd_cmpsit_dt *dt = data;

	return cmpsit_simul(&dt->cp);
}

#define doc_cmpsit                                                             \
	"Run composite algorithm (partially randomized), 2nd order Trotter."
#define argv0_cmpsit "ph2run [OPTS] cmpsit"
#define args_doc_cmpsit ""

static struct cmpsit_data *const args_cmpsit = &cmd_cmpsit_dt.cp_dt;

static struct argp_option opts_cmpsit[] = {
	{ "length", 'l', "VAL", 0, "Deterministic legth (default: 1)", 0 },
	{ "depth", 'd', "VAL", 0, "Random depth (default: 64)", 0 },
	{ "steps", 's', "N", 0, "Number of Trotter steps (default: 1)", 0 },
	{ "angle-det", 'D', "VAL", 0,
		"Deterministic angle, a floating point number (default: 1.0)", 0 },
	{ "angle-rand", 'R', "VAL", 0,
		"Angle for the randomized part, a floating point number (default: 1.0)",
		0 },
	{ "samples", 'n', "N", 0, "Number of independent samples (default: 1)",
		0 },
	{ "seed", 'x', "N", 0,
		"Seed of the pseudo random number generator (default: 1, must differ from zero)",
		0 },
	{ 0 }
};

static error_t opts_parser_cmpsit(int key, char *arg, struct argp_state *state)
{
	struct cmpsit_data *dt = state->input;

	switch (key) {
	case 's':
		dt->steps = strtoull(arg, nullptr, 10);
		break;
	case 'D':
		dt->angle_det = strtod(arg, nullptr);
		break;
	case 'R':
		dt->angle_rand = strtod(arg, nullptr);
		break;
	case 'd':
		dt->depth = strtoull(arg, nullptr, 10);
		break;
	case 'l':
		dt->length = strtoull(arg, nullptr, 10);
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

static struct argp argp_cmpsit = { opts_cmpsit, opts_parser_cmpsit,
	args_doc_cmpsit, doc_cmpsit, 0, 0, 0 };

int cmd_cmpsit(void)
{
	int rt = -1;

	if (data_exec(cmd_cmpsit_init, &cmd_cmpsit_dt) < 0)
		goto ex;
	if (timeit(cmd_cmpsit_run, &cmd_cmpsit_dt, &cmd_cmpsit_dt.t_tot) < 0) {
		log_error("Simulation error");
		goto ex;
	}
	if (data_exec(cmd_cmpsit_write, &cmd_cmpsit_dt) < 0)
		goto ex;

	rt = 0; /* Success. */
ex:
	log_info("Shut down simulation environment");

	return rt;
}

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
		cmd_parse(cmpsit, CMD_CMPSIT);
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
	case CMD_CMPSIT:
		rt = cmd_cmpsit();
		break;
	default:
		fprintf(stderr, "Unrecognized command.\n");
	}

	world_free();
	exit(rt);
}
