#include <stdlib.h>
#include <string.h>

#ifdef DISTRIBUTED

#include "mpi.h"

#endif

#include "log/log.h"

#include "circ.h"
#include "rayon.h"
#include "linen.h"
#include "data.h"

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"

#define PHASE2_DEFAULT_H5FILE "simul.h5"

void help_page(const int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s CIRCUIT [SIMUL_FILE_H5]\n\n", argv[0]);
	fprintf(stderr, "where CIRCUIT is must be one of:\n"
			"\n"
			"    linen\n"
			"    rayon\n"
			"\n"
			"If no simulation input file (HDF5) is specified,\n"
			"the default is " PHASE2_DEFAULT_H5FILE
			" in the current directory.\n");
	fprintf(stderr, "\n"
			"To set logging level, set environment variable:\n"
			"\n    " PHASE2_LOG_ENVVAR
			"={trace, debug, info, warn, error, fatal}"
			"\n\n");
}

static void log_callback(struct log_event *ev)
{
#ifdef DISTRIBUTED
	int initialized, finalized;
	MPI_Initialized(&initialized);
	MPI_Finalized(&finalized);
	if (initialized && !finalized) {
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank > 0) {
			return;
		}
	}
#endif
	char buf[64];
	FILE *fd = ev->data;
	buf[strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", ev->time)] = '\0';
	fprintf(fd, "%s %-5s ", buf, log_level_string(ev->level));
	vfprintf(fd, ev->fmt, ev->ap);
	fprintf(fd, "\n");
	fflush(fd);
}

int run_linen(struct circ_env *env, struct data *data)
{
	(void)(data);

	return linen_simulate(env);
}

int run_rayon(struct circ_env *env, struct data *dat)
{
	int rc;

	struct rayon_data rd;
	rayon_data_init(&rd);
	if (rayon_data_from_data(&rd, dat) < 0)
		return -1;
	rc = rayon_simulate(env, &rd);
	rayon_data_write_times(&(*dat).time_series, &rd.times);
	rayon_data_destroy(&rd);

	return rc;
}

int main(const int argc, char **argv)
{
	enum log_level lvl;
	const char *lvl_str = getenv(PHASE2_LOG_ENVVAR);
	if (!lvl_str || log_level_from_lowercase(&lvl, lvl_str) < 0) {
		lvl = LOG_ERROR;
	}
	log_set_level(lvl);
	log_add_callback(log_callback, stderr, lvl);

	struct circ_env env;
	if (circ_env_init(&env) != 0) {
		log_error("Failure: %s", "initialize environment");
		exit(EXIT_FAILURE);
	}

	log_info("*** Init ***");
#ifdef DISTRIBUTED
	int rank, num_ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	log_info("MPI num_ranks: %d", num_ranks);
	log_info("This is rank no. %d", rank);
#else
	log_info("MPI mode not enabled.");
#endif

	log_debug("Parsing command line arguments");
	if (argc < 2) {
		help_page(argc, argv);
		return EXIT_FAILURE;
	}
	const char *dat_filename = PHASE2_DEFAULT_H5FILE;
	if (argc < 3)
		log_debug("No simulation input file specified; "
			  "using default: %s",
			  PHASE2_DEFAULT_H5FILE);
	else
		dat_filename = argv[2];

	log_debug("Read simulation input file: %s", dat_filename);
	data_id fid = data_file_open(dat_filename);

	struct data dat;
	data_init(&dat);
	if (data_parse(&dat, fid) != 0) {
		log_error("Failure: read data file");
		exit(EXIT_FAILURE);
	}
	data_file_close(fid);

	log_debug("State preparation:");
	log_debug("multidet, num_qubits=%zu, num_terms=%zu",
		  dat.state_prep.multidet.num_qubits,
		  dat.state_prep.multidet.num_terms);
	log_debug("Hamiltonian: num_qubits=%zu, num_terms=%zu, "
		  "norm=%f",
		  dat.pauli_hamil.num_qubits, dat.pauli_hamil.num_terms,
		  dat.pauli_hamil.norm);
	log_debug("Time series: num_steps=%zu", dat.time_series.num_steps);

	log_info("*** Circuit ***");
	int rc;
	if (strncmp(argv[1], "linen", 5) == 0) {
		log_info("Circuit: linen");
		rc = run_linen(&env, &dat);
	} else if (strncmp(argv[1], "rayon", 5) == 0) {
		log_info("Circuit: rayon");
		rc = run_rayon(&env, &dat);
	} else {
		log_error("No circ named %s", argv[1]);
		exit(EXIT_FAILURE);
	}
	if (rc < 0) {
		log_error("Failure: simulation error");
		exit(EXIT_FAILURE);
	}

	log_debug("Saving data");
	fid = data_file_open(dat_filename);
	data_time_series_write(fid, &dat.time_series);
	data_file_close(fid);

	log_info("*** Cleanup ***");
	data_destroy(&dat);

	log_info("Shut down simulation environment");
	log_info("Done");
	circ_env_destroy(&env);

	return EXIT_SUCCESS;
}
