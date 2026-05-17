/*
 * b-algos -- end-to-end micro-benchmarks for the four
 * time-evolution algorithms under circ/.
 *
 * Each scenario builds a programmatic fixture
 * (Hamiltonian + multidet, identical across algorithms
 * so deltas are comparable), runs *_init outside the
 * timed loop, then times *_simul.  *_simul is
 * re-entrant: it calls circ_prepst at the top, which
 * zeroes the qreg, so repeated invocations measure
 * one fresh end-to-end run each.
 *
 * Where b-circ measures the per-step kernel
 * (circ_step, circ_measure) in isolation, b-algos
 * captures the whole-algorithm cost (state prep + step
 * loop + measurement) so end-to-end regressions in any
 * of those layers show up here.
 *
 * Total bench time on a quiet host: ~12 s
 * (~7 s MPI init + ~5 s measurements).
 */

#define LOG_SUBSYS "b-algos"

#include "bench.h"

#include <complex.h>
#include <stdlib.h>

#include "circ/cmpsit.h"
#include "circ/qdrift.h"
#include "circ/trott.h"
#include "circ/trott2.h"
#include "phase2.h"
#include "phase2/world.h"
#include "xoshiro256ss.h"

#define WD_SEED  UINT64_C(0x82e1d5f04ab937c6)
#define SEED     UINT64_C(0x9c4f3a261d80b75e)

#define NUM_RUNS 11
#define NQB      12
#define NTERMS   10
#define NDETS    10

static struct xoshiro256ss RNG;


/* -- fixture builders ------------------------------------------------------*/

static int build_hamil(struct circ_hamil *hm, uint32_t nqb, size_t nterms)
{
	if (circ_hamil_init(hm, nqb, nterms) < 0)
		return -1;
	for (size_t i = 0; i < nterms; i++) {
		hm->terms[i].op = paulis_new();
		for (uint32_t k = 0; k < nqb; k++)
			paulis_set(&hm->terms[i].op,
				(enum pauli_op)(xoshiro256ss_next(&RNG) % 4),
				k);
		hm->terms[i].cf = xoshiro256ss_dbl01(&RNG) - 0.5;
	}
	circ_hamil_sort_lex(hm);
	return 0;
}

static int build_muldet(struct circ_muldet *md, uint32_t nqb, size_t ndets)
{
	if (circ_muldet_init(md, ndets) < 0)
		return -1;
	const uint64_t mask = (nqb >= 64)
		? ~UINT64_C(0)
		: ((UINT64_C(1) << nqb) - 1);
	for (size_t i = 0; i < ndets; i++) {
		md->dets[i].idx = xoshiro256ss_next(&RNG) & mask;
		const double re = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
		const double im = xoshiro256ss_dbl01(&RNG) * 2.0 - 1.0;
		md->dets[i].cf = re + im * I;
	}
	return 0;
}

/* circ_init takes ownership of the carriers; clone them per
 * algorithm so each gets independent buffers. */
static int clone_hamil(struct circ_hamil *dst, const struct circ_hamil *src)
{
	if (circ_hamil_init(dst, src->qb, src->len) < 0)
		return -1;
	for (size_t i = 0; i < src->len; i++)
		dst->terms[i] = src->terms[i];
	return 0;
}

static int clone_muldet(struct circ_muldet *dst, const struct circ_muldet *src)
{
	if (circ_muldet_init(dst, src->len) < 0)
		return -1;
	for (size_t i = 0; i < src->len; i++)
		dst->dets[i] = src->dets[i];
	return 0;
}


/* -- table row + JSONL append ----------------------------------------------*/

static void record(FILE *out, const struct bench_prov *prov, int mpi_ranks,
	const char *name, const char *params_json, const char *display,
	double *samples, int num_runs, int sub_samples)
{
	const struct bench_stats st = bench_compute_stats(samples, num_runs);

	char path[256];
	bench_runs_path(prov, path, sizeof path);

	struct bench_baseline bl;
	const bool has_bl = bench_find_baseline(path, prov->hostname,
		name, params_json, sub_samples, &bl);

	bench_print_row(display, sub_samples, &st, &bl, has_bl);

	bench_append_jsonl(out, prov, BENCH_BACKEND, mpi_ranks,
		name, params_json, num_runs, sub_samples, &st);
}


/* -- per-algorithm scenarios ----------------------------------------------*/

static void measure_trott(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, const struct circ_hamil *hm_master,
	const struct circ_muldet *md_master)
{
	struct circ_hamil hm;
	struct circ_muldet md;
	if (clone_hamil(&hm, hm_master) < 0 ||
		clone_muldet(&md, md_master) < 0)
		return;

	const struct trott_data dt = { .steps = 5, .delta = 0.1 };
	struct trott tt;
	if (trott_init(&tt, &dt, hm, STPREP_MULTIDET, &md, NULL) < 0) {
		fprintf(stderr, "b-algos: trott_init failed\n");
		return;
	}

	trott_simul(&tt); /* warm-up */

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++)
		BENCH_SAMPLE_MIN(samples[r], 10, trott_simul(&tt));

	record(out, prov, mpi_ranks, "trott",
		"{\"steps\":5,\"delta\":0.1}",
		"trott steps=5 d=0.1",
		samples, NUM_RUNS, 10);

	trott_free(&tt);
}

static void measure_trott2(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, const struct circ_hamil *hm_master,
	const struct circ_muldet *md_master)
{
	struct circ_hamil hm;
	struct circ_muldet md;
	if (clone_hamil(&hm, hm_master) < 0 ||
		clone_muldet(&md, md_master) < 0)
		return;

	const struct trott2_data dt = { .steps = 5, .delta = 0.1 };
	struct trott2 t2;
	if (trott2_init(&t2, &dt, hm, STPREP_MULTIDET, &md, NULL) < 0) {
		fprintf(stderr, "b-algos: trott2_init failed\n");
		return;
	}

	trott2_simul(&t2);

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++)
		BENCH_SAMPLE_MIN(samples[r], 10, trott2_simul(&t2));

	record(out, prov, mpi_ranks, "trott2",
		"{\"steps\":5,\"delta\":0.1}",
		"trott2 steps=5 d=0.1",
		samples, NUM_RUNS, 10);

	trott2_free(&t2);
}

static void measure_qdrift(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, const struct circ_hamil *hm_master,
	const struct circ_muldet *md_master)
{
	struct circ_hamil hm;
	struct circ_muldet md;
	if (clone_hamil(&hm, hm_master) < 0 ||
		clone_muldet(&md, md_master) < 0)
		return;

	const struct qdrift_data dt = {
		.samples = 10,
		.depth = 16,
		.step_size = 0.1,
		.seed = UINT64_C(0xb40c2e1a8df37562),
	};
	struct qdrift qd;
	if (qdrift_init(&qd, &dt, hm, STPREP_MULTIDET, &md, NULL) < 0) {
		fprintf(stderr, "b-algos: qdrift_init failed\n");
		return;
	}

	qdrift_simul(&qd);

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++)
		BENCH_SAMPLE_MIN(samples[r], 5, qdrift_simul(&qd));

	record(out, prov, mpi_ranks, "qdrift",
		"{\"samples\":10,\"depth\":16,\"step_size\":0.1}",
		"qdrift n=10 d=16 s=0.1",
		samples, NUM_RUNS, 5);

	qdrift_free(&qd);
}

static void measure_cmpsit(FILE *out, const struct bench_prov *prov,
	int mpi_ranks, const struct circ_hamil *hm_master,
	const struct circ_muldet *md_master)
{
	struct circ_hamil hm;
	struct circ_muldet md;
	if (clone_hamil(&hm, hm_master) < 0 ||
		clone_muldet(&md, md_master) < 0)
		return;

	const struct cmpsit_data dt = {
		.samples = 10,
		.steps = 5,
		.length = 2,
		.depth = 16,
		.angle_det = 0.1,
		.angle_rand = 0.1,
		.seed = UINT64_C(0xa31df80c45e26b97),
	};
	struct cmpsit cp;
	if (cmpsit_init(&cp, &dt, hm, STPREP_MULTIDET, &md, NULL) < 0) {
		fprintf(stderr, "b-algos: cmpsit_init failed\n");
		return;
	}

	cmpsit_simul(&cp);

	double samples[NUM_RUNS];
	for (int r = 0; r < NUM_RUNS; r++)
		BENCH_SAMPLE_MIN(samples[r], 5, cmpsit_simul(&cp));

	record(out, prov, mpi_ranks, "cmpsit",
		"{\"samples\":10,\"steps\":5,\"length\":2,\"depth\":16}",
		"cmpsit n=10 s=5 l=2 d=16",
		samples, NUM_RUNS, 5);

	cmpsit_free(&cp);
}


/* -- main ------------------------------------------------------------------*/

int main(void)
{
	setenv("PHASE2_LOG", "warn", 0);

	bench_pin_cpu(0);

	world_init(nullptr, nullptr, WD_SEED);
	struct world_info wd;
	world_info(&wd);
	xoshiro256ss_init(&RNG, SEED);

	struct bench_prov prov;
	bench_prov_init(&prov);

	FILE *out = NULL;
	if (wd.rank == 0) {
		char path[256];
		if (bench_runs_path(&prov, path, sizeof path) < 0) {
			fprintf(stderr, "b-algos: cannot create bench/runs/\n");
			world_free();
			return 1;
		}
		out = fopen(path, "a");
		if (!out) {
			fprintf(stderr, "b-algos: cannot open %s\n", path);
			world_free();
			return 1;
		}
		bench_print_banner("b-algos (nqb=12)");
		bench_print_header();
	}

	/* Master fixture; each measure_* clones it. */
	struct circ_hamil hm_master;
	struct circ_muldet md_master;
	if (build_hamil(&hm_master, NQB, NTERMS) < 0 ||
		build_muldet(&md_master, NQB, NDETS) < 0) {
		fprintf(stderr, "b-algos: fixture build failed\n");
		if (out) fclose(out);
		world_free();
		return 1;
	}

	measure_trott (out, &prov, wd.size, &hm_master, &md_master);
	measure_trott2(out, &prov, wd.size, &hm_master, &md_master);
	measure_qdrift(out, &prov, wd.size, &hm_master, &md_master);
	measure_cmpsit(out, &prov, wd.size, &hm_master, &md_master);

	circ_hamil_free(&hm_master);
	circ_muldet_free(&md_master);

	if (out)
		fclose(out);
	world_free();
	return 0;
}
