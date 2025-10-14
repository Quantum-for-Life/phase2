#include "c23_compat.h"

#include <argp.h>
#include <stdlib.h>

#define xstr(s) str(s)
#define str(s) #s

const char *argp_program_version = PHASE2_VERSION;
const char *argp_program_bug_address = "Marek Miller <mlm@math.ku.dk>";

#define doc "Run phase2 sumilations."
#define args_doc "CMD [CMD_OPT]"

static struct argp_option opts[] = {
	{ "verbose", 'v', 0, 0, "Print verbose output", 0 },
	     { 0 }
};

static struct args {
	bool verbose;
	char *cmd;
} args;

static error_t opts_parser(int key, char *arg, struct argp_state *state)
{
	struct args *args = state->input;

	switch (key) {
	case 'v':
		args->verbose = true;
		break;

	case ARGP_KEY_ARG:
		args->cmd = arg;
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

int main(int argc, char **argv)
{
	argp_parse(&argp, argc, argv, ARGP_IN_ORDER, nullptr, &args);

	printf("VERBOSE = %s\n", args.verbose ? "yes" : "no");
	printf("CMD = %s\n", args.cmd);

	exit(0);
}
