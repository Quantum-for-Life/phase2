#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "error.h"
#include "ev.h"
#include "log.h"
#include "paulis.h"
#include "qreg.h"

static int MAIN_RET = 0;

#define FATAL_ERROR(...)                                                       \
	{                                                                      \
		log_error(__VA_ARGS__);                                        \
		MAIN_RET = EXIT_FAILURE;                                       \
		goto exit;                                                     \
	}

#define NUM_RND_PAULIS (111)
static int RND_PAULI_STATE = 0;
// clang-format off
static enum pauli_op RND_PAULIS[NUM_RND_PAULIS] = {
	0 , 1 , 1 , 1 , 3 , 0 , 0 , 0 , 2 , 0 , 3 , 2 , 2 , 0 , 3 , 1 , 1 , 0 ,
	0 , 3 , 2 , 3 , 3 , 2 , 0 , 3 , 0 , 0 , 2 , 0 , 3 , 3 , 3 , 1 , 1 , 0 ,
	0 , 0 , 1 , 0 , 0 , 2 , 2 , 0 , 2 , 0 , 2 , 3 , 3 , 1 , 2 , 1 , 0 , 0 ,
	0 , 1 , 1 , 0 , 0 , 0 , 2 , 2 , 2 , 0 , 0 , 3 , 0 , 2 , 1 , 0 , 1 , 1 ,
	0 , 3 , 3 , 1 , 3 , 1 , 2 , 3 , 1 , 0 , 0 , 0 , 1 , 0 , 1 , 2 , 2 , 2 ,
	3 , 3 , 1 , 3 , 0 , 3 , 2 , 1 , 1 , 2 , 3 , 1 , 0 , 0 , 2 , 0 , 3 , 1 ,
	1 , 1 , 3 ,
};
// clang-format on

static enum pauli_op rnd_pauli(void)
{
	RND_PAULI_STATE = (RND_PAULI_STATE + 13) % NUM_RND_PAULIS;
	return RND_PAULIS[RND_PAULI_STATE];
}

static fl RND_ANGLE = 0.3;

static fl rnd_angle(void)
{
	RND_ANGLE = RND_ANGLE * -1.142 + 0.712;
	return RND_ANGLE;
}

static struct {
	unsigned long	   num_qubits;
	unsigned long long num_rots;
} opts;

void print_help_page(int argc, char **argv);
int  parse_args(int argc, char **argv);

int main(const int argc, char **argv)
{
	if (argc == 1) {
		print_help_page(argc, argv);
		exit(EXIT_SUCCESS);
	}

	struct ev ev;
	if (ev_init(&ev) < 0 || log_init() < 0) exit(EXIT_FAILURE);
	log_info("*** Init ***");

	if (parse_args(argc, argv) < 0) {
		log_error("cannot parse command line arguments");
		print_help_page(argc, argv);
		exit(EXIT_FAILURE);
	}
	log_info("num_qubits: %lu", opts.num_qubits);
	log_info("num_rots: %Lu", opts.num_rots);

	log_info("MPI num_ranks: %d", ev.num_ranks);
	log_info("This is rank no. %d", ev.rank);

	const unsigned num_ranks = ev.num_ranks;
	if ((num_ranks & (num_ranks - 1)) != 0)
		FATAL_ERROR("number of MPI ranks must be a power of two");

	struct qreg reg;
	if (qreg_init(&reg, opts.num_qubits, &ev) < 0) FATAL_ERROR("init qreg");
	qreg_zerostate(&reg);

	log_info("num_reqs=%zu, msg_count=%d", reg.num_reqs, reg.msg_count);

	log_info("*** multi pauli rotation [x1] ***");
	{
		const clock_t start = clock();
		for (size_t i = 0; i < opts.num_rots; i++) {
			log_trace("multi pauli rotation [%d]", i);

			struct paulis code_hi = paulis_new();
			for (u32 k = 0; k < reg.qb_hi; k++)
				paulis_set(&code_hi, rnd_pauli(), k);

			fl	      angles[1];
			struct paulis codes_lo[1];
			for (int j = 0; j < 1; j++) {
				codes_lo[j] = paulis_new();
				for (u32 k = 0; k < reg.qb_lo; k++)
					paulis_set(
						&codes_lo[j], rnd_pauli(), k);
				angles[j] = rnd_angle();
			}

			qreg_paulirot(&reg, code_hi, codes_lo, angles, 1);
		}
		const double elapsed =
			(double)(clock() - start) / CLOCKS_PER_SEC;
		log_info("time per pauli rotation: %f [s]",
			elapsed / opts.num_rots / 1);
	}

#define NUM_CODES_LO (16)
	log_info("*** multi pauli rotation [x%d] ***", NUM_CODES_LO);
	{
		const clock_t start = clock();
		for (size_t i = 0; i < opts.num_rots; i++) {
			log_trace("multi pauli rotation [%d]", i);

			struct paulis code_hi = paulis_new();
			for (u32 k = 0; k < reg.qb_hi; k++)
				paulis_set(&code_hi, rnd_pauli(), k);

			fl angles[NUM_CODES_LO];

			struct paulis codes_lo[NUM_CODES_LO];
			for (int j = 0; j < NUM_CODES_LO; j++) {
				codes_lo[j] = paulis_new();
				for (u32 k = 0; k < reg.qb_lo; k++)
					paulis_set(
						&codes_lo[j], rnd_pauli(), k);
				angles[j] = rnd_angle();
			}

			qreg_paulirot(
				&reg, code_hi, codes_lo, angles, NUM_CODES_LO);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		const double elapsed =
			(double)(clock() - start) / CLOCKS_PER_SEC;

		log_info("time per pauli rotation: %f [s]",
			elapsed / opts.num_rots / NUM_CODES_LO);
	}

	qreg_destroy(&reg);
exit:
	log_info("*** Shutdown ***");

	if (log_shutdown() < 0 || ev_destroy(&ev) < 0) exit(EXIT_FAILURE);

	return MAIN_RET;
}

void print_help_page(int argc, char **argv)
{
	(void)argc;
	fprintf(stderr, "usage: %s NUM_QUBITS NUM_ROTATIONS\n\n", argv[0]);
}

int parse_args(int argc, char **argv)
{
	if (argc < 3) return -EARGS;

	opts.num_qubits = strtoul(argv[1], NULL, 10);
	if (opts.num_qubits == 0) return -EARGS;

	opts.num_rots = strtoull(argv[2], NULL, 10);
	if (opts.num_rots == 0) return -EARGS;

	return OK;
}
