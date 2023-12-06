#include <stdlib.h>
#include <string.h>

#ifdef DISTRIBUTED

#include "mpi.h"

#endif

#include "log.h"

#include "circ.h"
#include "rayon.h"
#include "linen.h"
#include "data.h"

#define PHASE2_LOG_ENVVAR "PHASE2_LOG"
#define PHASE2_LOG_FILE "simul.log"

#define PHASE2_DEFAULT_H5FILE "simul.h5"

void exit_failure(const char *msg)
{
	log_error("Failure: %s", msg);
	exit(EXIT_FAILURE);
}

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

void set_log_level()
{
	const char *log_level = getenv(PHASE2_LOG_ENVVAR);
	if (log_level == NULL) {
		log_set_level(LOG_ERROR);
		return;
	}
	if (strncmp(log_level, "trace", 5) == 0) {
		log_set_level(LOG_TRACE);
	}
	if (strncmp(log_level, "debug", 5) == 0) {
		log_set_level(LOG_DEBUG);
	}
	if (strncmp(log_level, "info", 4) == 0) {
		log_set_level(LOG_INFO);
	}
	if (strncmp(log_level, "warn", 4) == 0) {
		log_set_level(LOG_WARN);
	}
	if (strncmp(log_level, "error", 5) == 0) {
		log_set_level(LOG_ERROR);
	}
	if (strncmp(log_level, "fatal", 5) == 0) {
		log_set_level(LOG_FATAL);
	}
}

int main(const int argc, char **argv)
{
	data_id fid;
	struct circ_env env;
	if (circ_env_init(&env) != 0) {
		exit_failure("initialize environment");
	}

#ifdef DISTRIBUTED
	int rank, num_ranks;
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	log_info("*** Init ***");
	log_info("Initialize MPI environment");
	log_info("MPI num_ranks: %d", num_ranks);
	log_info("This is rank no. %d", rank);
#else
	log_info("*** Init ***");
	log_info("MPI mode not enabled.");
	log_info("To enable distributed mode, set "
		 "-DDISTRIBUTED "
		 "flag during compilation");
#endif

	set_log_level();
	log_info("Open log file: " PHASE2_LOG_FILE);
	FILE *log_file = fopen(PHASE2_LOG_FILE, "a");
	if (!log_file) {
		exit_failure("open log file");
	}
	log_add_fp(log_file, LOG_DEBUG);

	log_debug("Parsing command line arguments");
	if (argc < 2) {
		help_page(argc, argv);
		return EXIT_FAILURE;
	}
	const char *dat_filename = PHASE2_DEFAULT_H5FILE;
	if (argc < 3) {
		log_debug("No simulation input file specified; "
			  "using default: %s",
			  PHASE2_DEFAULT_H5FILE);
	} else {
		dat_filename = argv[2];
	}

	log_debug("Read simulation input file: %s", dat_filename);
	fid = data_file_open(dat_filename);

	struct data dat;
	data_init(&dat);
	if (data_parse(&dat, fid) != DATA_OK) {
		exit_failure("read data file");
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
	int sucess;
	if (strncmp(argv[1], "linen", 5) == 0) {
		log_info("Circuit: linen");
		sucess = linen_simulate(&env) == 0;
	} else if (strncmp(argv[1], "rayon", 5) == 0) {
		log_info("Circuit: rayon");
		struct rayon_data rd;
		rayon_data_init(&rd);
		rayon_data_from_data(&rd, &dat);
		sucess = rayon_simulate(&env, &rd) == 0;
		rayon_data_write_times(&dat.time_series, &rd.times);
		rayon_data_destroy(&rd);

	} else {
		log_error("No circ named %s", argv[1]);
		sucess = 0;
	}
	if (!sucess) {
		exit_failure("simulation error");
	}

	log_debug("Saving data");
	fid = data_file_open(dat_filename);
	data_time_series_write(fid, &dat.time_series);
	data_file_close(fid);

	log_info("*** Cleanup ***");
	data_destroy(&dat);

	log_info("Shut down simulation environment");
	circ_env_destroy(&env);

	log_info("Done");
	fclose(log_file);
	return EXIT_SUCCESS;
}
