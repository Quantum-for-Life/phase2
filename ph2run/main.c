#include <stdlib.h>
#include <string.h>

#ifdef DISTRIBUTED

#include "mpi.h"

#endif

#include "circ.h"
#include "data.h"

#include "log.h"

/* Command line arguments and env vars */
#include "argums.h"

static struct argums argums;
void print_help_page(int argc, char **argv);
int argums_parse(struct argums *a, int argc, char **argv);

/* Data */
static struct data dat;
void print_data_info(const struct data *dat);
int read_data_file(struct data *dat, const char *filename);
int save_data(const char *filename, const struct data *dat);

/* Runners */
int run_linen(struct circ_env *env, struct data *dat);
int run_rayon(struct circ_env *env, struct data *dat);

int main(const int argc, char **argv)
{
	int rc = 0;

	/* Parse command line arguments. */
	if (argums_parse(&argums, argc, argv) < 0)
		exit(EXIT_FAILURE);

	/* Initiallize logging */
	if (log_init() < 0)
		exit(EXIT_FAILURE);

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

	data_init(&dat);
	if (read_data_file(&dat, argums.dat_filename) < 0)
		goto error;
	print_data_info(&dat);

	log_info("*** Circuit ***");
	if (strncmp(argv[1], "linen", 5) == 0) {
		log_info("Circuit: linen");
		if (run_linen(&env, &dat) < 0)
			goto error;
	} else if (strncmp(argv[1], "rayon", 5) == 0) {
		log_info("Circuit: rayon");
		if (run_rayon(&env, &dat) < 0) {
			log_error("Failure: simulation error");
			goto error;
		}
	} else {
		log_error("No circ named %s", argv[1]);
		goto error;
	}

	if (save_data(argums.dat_filename, &dat) < 0)
		goto error;

	goto cleanup;
error:
	rc = -1;
cleanup:
	log_info("Shut down simulation environment");
	data_destroy(&dat);
	circ_env_destroy(&env);

	return rc;
}
