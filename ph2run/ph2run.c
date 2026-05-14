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
#include "circ/trott2.h"
#define LOG_SUBSYS "ph2run"
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
	"Run phase2 simulations.  CMD can be one of the algorithms:\n"         \
	"  trott    1st-order Trotter product formula\n"                       \
	"  trott2   2nd-order symmetric (Strang) Trotter\n"                    \
	"  qdrift   qDRIFT randomised product formula\n"                       \
	"  cmpsit   composite (deterministic + randomised, 2nd order)\n"       \
	"\nRun CMD --help for more information."
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
	CMD_TROTT2 = 3,
};

/* Command: "trott" */
static struct cmd_trott_dt {
	struct trott tt;
	struct trott_data tt_dt;
	double t_tot;
} cmd_trott_dt = {
	.tt_dt = { .delta = 1.0, .steps = 1 }
};

static int cmd_trott_run(void *data)
{
	struct cmd_trott_dt *dt = data;

	log_info("running simulation: %zu Trotter steps", dt->tt_dt.steps);
	return trott_simul(&dt->tt);
}

static void cmd_trott_summary(const struct cmd_trott_dt *dt)
{
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,delta,n_steps,n_ranks,t_tot");
	log_info("> %u,%zu,%zu,%f,%zu,%d,%.3f", dt->tt.ct.hm.qb,
		dt->tt.ct.hm.len, dt->tt.ct.md.len, dt->tt_dt.delta,
		dt->tt_dt.steps, wd.size, dt->t_tot);
}

#define doc_trott "Run 1st-order Trotter product formula."
#define argv0_trott "ph2run [OPTS] trott"
#define args_doc_trott ""

/* Alias for the option parser. */
static struct trott_data *const args_trott = &cmd_trott_dt.tt_dt;

static struct argp_option opts_trott[] = {
	{ "delta", 'D', "VAL", 0,
		"Trotter step size (default: 1.0)", 0 },
	{ "steps", 's', "N", 0,
		"Number of Trotter steps (default: 1)", 0 },
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
	log_info("open data file: %s", args.simul);
	const data_id fid = data_open(args.simul);
	if (fid == DATA_INVALID_FID)
		goto ex_open;

	log_info("*** Circuit: trott ***");
	log_info("delta: %f", cmd_trott_dt.tt_dt.delta);
	log_info("num_steps: %zu", cmd_trott_dt.tt_dt.steps);
	if (trott_init(&cmd_trott_dt.tt, &cmd_trott_dt.tt_dt, fid) < 0) {
		log_error("trott: init failed");
		goto ex_close;
	}
	if (timeit(cmd_trott_run, &cmd_trott_dt, &cmd_trott_dt.t_tot) < 0) {
		log_error("trott: simulation failed after %.3f s",
			cmd_trott_dt.t_tot);
		goto ex_free;
	}
	log_info("simulation finished in %.3f s", cmd_trott_dt.t_tot);
	cmd_trott_summary(&cmd_trott_dt);

	rt = 0;
ex_free:
	trott_free(&cmd_trott_dt.tt);
ex_close:
	log_info("close data file: %s", args.simul);
	data_close(fid);
ex_open:
	log_info("Shut down simulation environment");
	return rt;
}

/* Command: "trott2" */
static struct cmd_trott2_dt {
	struct trott2 t2;
	struct trott2_data t2_dt;
	double t_tot;
} cmd_trott2_dt = {
	.t2_dt = { .delta = 1.0, .steps = 1 }
};

static int cmd_trott2_run(void *data)
{
	struct cmd_trott2_dt *dt = data;

	log_info("running simulation: %zu symmetric Trotter steps",
		dt->t2_dt.steps);
	return trott2_simul(&dt->t2);
}

static void cmd_trott2_summary(const struct cmd_trott2_dt *dt)
{
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,delta,n_steps,n_ranks,t_tot");
	log_info("> %u,%zu,%zu,%f,%zu,%d,%.3f", dt->t2.ct.hm.qb,
		dt->t2.ct.hm.len, dt->t2.ct.md.len, dt->t2_dt.delta,
		dt->t2_dt.steps, wd.size, dt->t_tot);
}

#define doc_trott2 "Run 2nd-order symmetric (Strang) Trotter product formula."
#define argv0_trott2 "ph2run [OPTS] trott2"
#define args_doc_trott2 ""

static struct trott2_data *const args_trott2 = &cmd_trott2_dt.t2_dt;

static struct argp_option opts_trott2[] = {
	{ "delta", 'D', "VAL", 0,
		"Trotter step size (default: 1.0)", 0 },
	{ "steps", 's', "N", 0,
		"Number of Trotter steps (default: 1)", 0 },
	{ 0 }
};

static error_t opts_parser_trott2(int key, char *arg, struct argp_state *state)
{
	struct trott2_data *dt = state->input;

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

static struct argp argp_trott2 = { opts_trott2, opts_parser_trott2,
	args_doc_trott2, doc_trott2, 0, 0, 0 };

static int cmd_trott2(void)
{
	int rt = -1;
	log_info("open data file: %s", args.simul);
	const data_id fid = data_open(args.simul);
	if (fid == DATA_INVALID_FID)
		goto ex_open;

	log_info("*** Circuit: trott2 (Strang 2nd-order) ***");
	log_info("delta: %f", cmd_trott2_dt.t2_dt.delta);
	log_info("num_steps: %zu", cmd_trott2_dt.t2_dt.steps);
	if (trott2_init(&cmd_trott2_dt.t2, &cmd_trott2_dt.t2_dt, fid) < 0) {
		log_error("trott2: init failed");
		goto ex_close;
	}
	if (timeit(cmd_trott2_run, &cmd_trott2_dt, &cmd_trott2_dt.t_tot) < 0) {
		log_error("trott2: simulation failed after %.3f s",
			cmd_trott2_dt.t_tot);
		goto ex_free;
	}
	log_info("simulation finished in %.3f s", cmd_trott2_dt.t_tot);
	cmd_trott2_summary(&cmd_trott2_dt);

	rt = 0;
ex_free:
	trott2_free(&cmd_trott2_dt.t2);
ex_close:
	log_info("close data file: %s", args.simul);
	data_close(fid);
ex_open:
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

static int cmd_qdrift_run(void *data)
{
	struct cmd_qdrift_dt *dt = data;

	log_info("running simulation: %zu samples, depth %zu",
		dt->qd_dt.samples, dt->qd_dt.depth);
	return qdrift_simul(&dt->qd);
}

static void cmd_qdrift_summary(const struct cmd_qdrift_dt *dt)
{
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,n_samples,step_size,depth,"
		 "n_ranks,t_tot");
	log_info("> %u,%zu,%zu,%zu,%.6f,%zu,%d,%.3f", dt->qd.ct.hm.qb,
		dt->qd.ct.hm.len, dt->qd.ct.md.len, dt->qd_dt.samples,
		dt->qd_dt.step_size, dt->qd_dt.depth, wd.size, dt->t_tot);
}

#define doc_qdrift "Run qDRIFT randomised product formula."
#define argv0_qdrift "ph2run [OPTS] qdrift"
#define args_doc_qdrift ""

static struct qdrift_data *const args_qdrift = &cmd_qdrift_dt.qd_dt;

static struct argp_option opts_qdrift[] = {
	{ "step-size", 'D', "VAL", 0,
		"qDRIFT step size (default: 1.0)", 0 },
	{ "depth", 'd', "N", 0,
		"Number of randomised terms per sample (default: 64)", 0 },
	{ "samples", 'n', "N", 0,
		"Number of independent samples (default: 1)", 0 },
	{ "seed", 'x', "N", 0,
		"PRNG seed; must be non-zero (default: 1)", 0 },
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
	log_info("open data file: %s", args.simul);
	const data_id fid = data_open(args.simul);
	if (fid == DATA_INVALID_FID)
		goto ex_open;

	log_info("*** Circuit: qdrift ***");
	log_info("step_size: %f", cmd_qdrift_dt.qd_dt.step_size);
	log_info("depth: %zu", cmd_qdrift_dt.qd_dt.depth);
	log_info("samples: %zu", cmd_qdrift_dt.qd_dt.samples);
	log_info("seed: %lu", cmd_qdrift_dt.qd_dt.seed);
	if (qdrift_init(&cmd_qdrift_dt.qd, &cmd_qdrift_dt.qd_dt, fid) < 0) {
		log_error("qdrift: init failed");
		goto ex_close;
	}
	if (timeit(cmd_qdrift_run, &cmd_qdrift_dt, &cmd_qdrift_dt.t_tot) < 0) {
		log_error("qdrift: simulation failed after %.3f s",
			cmd_qdrift_dt.t_tot);
		goto ex_free;
	}
	log_info("simulation finished in %.3f s", cmd_qdrift_dt.t_tot);
	cmd_qdrift_summary(&cmd_qdrift_dt);

	rt = 0;
ex_free:
	qdrift_free(&cmd_qdrift_dt.qd);
ex_close:
	log_info("close data file: %s", args.simul);
	data_close(fid);
ex_open:
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

static int cmd_cmpsit_run(void *data)
{
	struct cmd_cmpsit_dt *dt = data;

	log_info("running simulation: %zu samples x %zu steps,"
		 " L=%zu depth=%zu",
		dt->cp_dt.samples, dt->cp_dt.steps, dt->cp_dt.length,
		dt->cp_dt.depth);
	return cmpsit_simul(&dt->cp);
}

static void cmd_cmpsit_summary(const struct cmd_cmpsit_dt *dt)
{
	log_info("> Simulation summary (CSV):");
	log_info("> n_qb,n_terms,n_dets,samples,length,depth,angle_det"
		 ",angle_rand,steps,n_ranks,t_tot");
	log_info("> %u,%zu,%zu,%zu,%zu,%zu,%.16f,%.16f,%zu,%d,%.3f",
		dt->cp.ct.hm.qb, dt->cp.ct.hm.len, dt->cp.ct.md.len,
		dt->cp_dt.samples, dt->cp_dt.length, dt->cp_dt.depth,
		dt->cp_dt.angle_det, dt->cp_dt.angle_rand, dt->cp_dt.steps,
		wd.size, dt->t_tot);
}

#define doc_cmpsit                                                             \
	"Run composite (partially randomised) 2nd-order Trotter."
#define argv0_cmpsit "ph2run [OPTS] cmpsit"
#define args_doc_cmpsit ""

static struct cmpsit_data *const args_cmpsit = &cmd_cmpsit_dt.cp_dt;

static struct argp_option opts_cmpsit[] = {
	{ "length", 'l', "N", 0,
		"Number of deterministic top-|c_k| terms (default: 1)", 0 },
	{ "depth", 'd', "N", 0,
		"Number of randomised terms per step (default: 64)", 0 },
	{ "steps", 's', "N", 0,
		"Number of Trotter steps (default: 1)", 0 },
	{ "angle-det", 'D', "VAL", 0,
		"Step size for the deterministic part (default: 1.0)", 0 },
	{ "angle-rand", 'R', "VAL", 0,
		"Step size for the randomised part (default: 1.0)", 0 },
	{ "samples", 'n', "N", 0,
		"Number of independent samples (default: 1)", 0 },
	{ "seed", 'x', "N", 0,
		"PRNG seed; must be non-zero (default: 1)", 0 },
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
	log_info("open data file: %s", args.simul);
	const data_id fid = data_open(args.simul);
	if (fid == DATA_INVALID_FID)
		goto ex_open;

	log_info("*** Circuit: cmpsit ***");
	log_info("seed: %lu", cmd_cmpsit_dt.cp_dt.seed);
	log_info("length: %zu", cmd_cmpsit_dt.cp_dt.length);
	log_info("depth: %zu", cmd_cmpsit_dt.cp_dt.depth);
	log_info("steps: %zu", cmd_cmpsit_dt.cp_dt.steps);
	log_info("angle_det: %.16f", cmd_cmpsit_dt.cp_dt.angle_det);
	log_info("angle_rand: %.16f", cmd_cmpsit_dt.cp_dt.angle_rand);
	log_info("samples: %zu", cmd_cmpsit_dt.cp_dt.samples);
	if (cmpsit_init(&cmd_cmpsit_dt.cp, &cmd_cmpsit_dt.cp_dt, fid) < 0) {
		log_error("cmpsit: init failed");
		goto ex_close;
	}
	if (timeit(cmd_cmpsit_run, &cmd_cmpsit_dt, &cmd_cmpsit_dt.t_tot) < 0) {
		log_error("cmpsit: simulation failed after %.3f s",
			cmd_cmpsit_dt.t_tot);
		goto ex_free;
	}
	log_info("simulation finished in %.3f s", cmd_cmpsit_dt.t_tot);
	cmd_cmpsit_summary(&cmd_cmpsit_dt);

	rt = 0;
ex_free:
	cmpsit_free(&cmd_cmpsit_dt.cp);
ex_close:
	log_info("close data file: %s", args.simul);
	data_close(fid);
ex_open:
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
		/* trott2 must precede trott: strncmp uses strlen("trott")
		 * so "trott" alone would otherwise swallow "trott2". */
		cmd_parse(trott2, CMD_TROTT2);
		cmd_parse(trott, CMD_TROTT);
		cmd_parse(qdrift, CMD_QDRIFT);
		cmd_parse(cmpsit, CMD_CMPSIT);
	}

	if (world_init(nullptr, nullptr, WD_SEED) != WORLD_READY) {
		fprintf(stderr, "ph2run: world_init failed\n");
		exit(-1);
	}

	if (cmd < 0) {
		log_error("unrecognised subcommand: %s",
			args.cmd ? args.cmd : "(none)");
		world_free();
		exit(-1);
	}

	log_info("subcommand: %s", args.cmd);
	log_info("simul file: %s", args.simul);

	switch (cmd) {
	case CMD_TROTT:
		rt = cmd_trott();
		break;
	case CMD_TROTT2:
		rt = cmd_trott2();
		break;
	case CMD_QDRIFT:
		rt = cmd_qdrift();
		break;
	case CMD_CMPSIT:
		rt = cmd_cmpsit();
		break;
	}

	world_free();
	exit(rt);
}
