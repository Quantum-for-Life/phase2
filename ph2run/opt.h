#ifndef PH2RUN_ARGS_H
#define PH2RUN_ARGS_H

#define PH2RUN_DEFAULT_H5FILE "simul.h5"
#define PH2RUN_LOG_ENVVAR "PHASE2_LOG"

#define PH2RUN_SILK_DEFAULT_NUM_STEPS (8)

enum opt_cicuit {
	OPT_CICUIT_LINEN,
	OPT_CICUIT_RAYON,
	OPT_CICUIT_SILK,
};

struct args_linen {};

struct args_rayon {};

struct args_silk {
	size_t num_steps;
};

struct opt {
	enum opt_cicuit cicuit;
	const char     *dat_filename;
	union {
		struct args_linen linen;
		struct args_rayon rayon;
		struct args_silk silk;
	} circuit_args;
};

#endif // PH2RUN_ARGS_H
