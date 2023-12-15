#include <stdlib.h>
#include <string.h>

#ifdef DISTRIBUTED

#include "mpi.h"

#endif

#include "circ.h"
#include "data.h"

#include "log.h"

/* Command line arguments and env vars */
#include "opt.h"

static struct opt opt;
void opt_help_page(int argc, char **argv);
int opt_parse(struct opt *o, int argc, char **argv);

/* Data */
static struct data dat;
void datio_print_info(const struct data *dat);
int datio_read_file(struct data *dat, const char *filename);
int datio_save_file(const char *filename, const struct data *dat);

/* Runners */
int run_linen(struct data *dat);
int run_rayon(struct data *dat);

int main(const int argc, char **argv)
{
	int rc = 0;

	/* Parse command line arguments. */
	if (opt_parse(&opt, argc, argv) < 0)
		exit(EXIT_FAILURE);

	circ_initialize();

	/* Initiallize logging */
	if (log_init() < 0)
		exit(EXIT_FAILURE);

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
	if (datio_read_file(&dat, opt.dat_filename) < 0)
		goto error;
	datio_print_info(&dat);

	log_info("*** Circuit ***");
	if (strncmp(argv[1], "linen", 5) == 0) {
		log_info("Circuit: linen");
		if (run_linen(&dat) < 0)
			goto error;
	} else if (strncmp(argv[1], "rayon", 5) == 0) {
		log_info("Circuit: rayon");
		if (run_rayon(&dat) < 0) {
			log_error("Failure: simulation error");
			goto error;
		}
	} else {
		log_error("No circ named %s", argv[1]);
		goto error;
	}

	if (datio_save_file(opt.dat_filename, &dat) < 0)
		goto error;

	goto cleanup;
error:
	rc = -1;
cleanup:
	log_info("Shut down simulation environment");
	data_destroy(&dat);
	circ_shutdown();

	return rc;
}
