#include <stdlib.h>
#include <string.h>


#include "mpi.h"

#include "circ.h"
#include "data2.h"

#include "log.h"
#include "opt.h"

static struct opt opt;

void opt_help_page(int argc, char **argv);
int  opt_parse(struct opt *o, int argc, char **argv);

/* Runners */
int run_linen(void);
int run_rayon(data2_id fid);
int run_silk(data2_id fid, size_t num_steps);

int main(const int argc, char **argv)
{
	int rc = 0;

	/* Parse command line arguments. */
	if (opt_parse(&opt, argc, argv) < 0)
		exit(EXIT_FAILURE);

	/* Initiallize logging */
	if (log_init() < 0)
		exit(EXIT_FAILURE);

	circ_initialize();

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

	data2_id fid = data2_open(opt.dat_filename);
	if (fid == DATA2_INVALID_FID)
		goto error;

	log_info("*** Circuit ***");
	switch (opt.cicuit) {
	case OPT_CICUIT_LINEN:
		log_info("Circuit: linen");
		if (run_linen() < 0)
			goto error;
		break;
	case OPT_CICUIT_RAYON:
		log_info("Circuit: rayon");
		if (run_rayon(fid) < 0) {
			log_error("Failure: simulation error");
			goto error;
		}
		break;
	case OPT_CICUIT_SILK:
		log_info("Circuit: silk");
		log_info("Num_steps: %zu", opt.circuit_args.silk.num_steps);
		if (run_silk(fid, opt.circuit_args.silk.num_steps) < 0) {
			log_error("Failure: simulation error");
			goto error;
		}
		break;
	}

	data2_close(fid);

	goto cleanup;
error:
	rc = -1;
cleanup:
	log_info("Shut down simulation environment");

	return rc;
}
