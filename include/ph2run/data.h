#ifndef PHASE2_DATA_H
#define PHASE2_DATA_H

#include <stddef.h>
#include <stdint.h>

struct circ_hamil;

#define DATA_INVALID_FID INT64_C(-1)

/*
 * Returned by data_open() on ranks != 0.  The data subsystem
 * is rank-0-only: rank 0 owns the file handle and performs
 * all H5 calls; other ranks receive read results via MPI_Bcast
 * and short-circuit on writes.  Callers should treat both
 * DATA_INVALID_FID (failure) and DATA_FOLLOWER_FID (success on
 * follower rank) as "do not pass to H5* directly"; the data_*
 * functions handle both internally.
 */
#define DATA_FOLLOWER_FID INT64_C(-2)

/*
 * Group names and per-group scalar attribute names that the
 * circuit algorithms reference when writing per-step output.
 * Internal sub-paths (the Hamiltonian dataset names, the
 * state-prep subgroups, the values dataset under each circ
 * group) live in phase2/data.c and never leak to callers.
 */
#define DATA_CIRCTROTT "circ_trott"
#define DATA_CIRCTROTT_DELTA "delta"

#define DATA_CIRCTROTT2 "circ_trott2"
#define DATA_CIRCTROTT2_DELTA "delta"

#define DATA_CIRCQDRIFT "circ_qdrift"
#define DATA_CIRCQDRIFT_DEPTH "depth"
#define DATA_CIRCQDRIFT_NUMSAMPLES "num_samples"
#define DATA_CIRCQDRIFT_STEPSIZE "step_size"
#define DATA_CIRCQDRIFT_SEED "seed"

#define DATA_CIRCCMPSIT "circ_cmpsit"
#define DATA_CIRCCMPSIT_DEPTH "depth"
#define DATA_CIRCCMPSIT_LENGTH "length"
#define DATA_CIRCCMPSIT_ANGLEDET "angle_det"
#define DATA_CIRCCMPSIT_ANGLERAND "angle_rand"
#define DATA_CIRCCMPSIT_STEPS "steps"
#define DATA_CIRCCMPSIT_SEED "seed"

/**
 * Handle to a data file
 */
typedef int64_t data_id;

/**
 * Open data file.
 *
 * Arguments:
 *
 *	filename	Path to the data file.
 *
 * Return value:
 *
 *	A valid data_id value or DATA_INVALID_FID in case of error
 */
data_id data_open(const char *filename);

/**
 * Close data file.
 *
 * Close the file that has been open with the call to data_open().  Free the
 * resources.
 */
void data_close(data_id);

int data_grp_create(data_id fid, const char *grp_name);

int data_attr_read_dbl(data_id fid, const char *grp_name,
	const char *attr_name, double *a);

#define data_attr_read(fid, grp_name, attr_name, attr_buf)                     \
	_Generic((attr_buf),                                                   \
		double *: data_attr_read_dbl)(                                 \
		fid, grp_name, attr_name, attr_buf)

/* Create an attribute and write the value of 'a' to it. */
int data_attr_write_ul(data_id fid, const char *grp_name,
	const char *attr_name, unsigned long a);
int data_attr_write_dbl(data_id fid, const char *grp_name,
	const char *attr_name, double a);

#define data_attr_write(fid, grp_name, attr_name, attr)                        \
	_Generic((attr),                                                       \
		double: data_attr_write_dbl,                                   \
		unsigned long: data_attr_write_ul)(                            \
		fid, grp_name, attr_name, attr)

/**
 * Load /state_prep/multidet into a packed circ_muldet.
 *
 * Rank 0 reads the dimensions, the complex coefficients and
 * the per-qubit determinant bytes from disk and broadcasts to
 * followers in one collective call.  The destination dets
 * array is allocated through circ_muldet_init() and each row
 * is packed in place into { uint64_t idx; _Complex double cf }
 * (LSB = qubit 0).  Determinant bytes are validated to be 0
 * or 1; a stray value rejects the load and leaves *md zeroed.
 *
 * Release with circ_muldet_free() (declared in phase2/circ.h).
 *
 * Returns 0 on success, -1 on error.
 */
struct circ_muldet;

int data_muldet_load(data_id fid, struct circ_muldet *md);

/*
 * State-prep dispatch and coefficient-matrix loader.
 *
 * `enum stprep_kind`, `struct data_coeff_block` and
 * `struct data_coeff_matrix` are carrier types consumed by
 * the circuit layer; their definitions live in
 * phase2/circ.h so the algorithm header owns the
 * in-memory representation.  data.h only forward-declares
 * them to keep the I/O surface here while avoiding an
 * include cycle with phase2/state_prep_coeff.h.
 *
 * data_state_prep_kind() inspects the file for
 * /state_prep/multidet and /state_prep/coeff_matrix:
 *
 *  multidet | coeff_matrix | result
 *  ---------+--------------+--------
 *  absent   | absent       | -ENOENT
 *  present  | absent       | STPREP_MULTIDET
 *  absent   | present      | STPREP_COEFF_MATRIX
 *  present  | present      | -EINVAL  (ambiguous; rebuild simul.h5)
 *
 * On success *out holds the selected kind and the function
 * returns 0.  On failure *out is unchanged.
 *
 * data_coeff_matrix_load() opens /state_prep/coeff_matrix,
 * reads all scalars and arrays, validates shapes, and
 * broadcasts to followers in one call.  The matching free
 * (data_coeff_matrix_free) is a pure pointer release and
 * lives in phase2/circ.c with the other in-memory carrier
 * destructors -- see phase2/circ.h.
 *
 * Documented further in phase2/doc/simul-h5-specs.md
 * "dispatch rules".
 */
enum stprep_kind;
struct data_coeff_matrix;

int data_state_prep_kind(data_id fid, enum stprep_kind *out);
int data_coeff_matrix_load(data_id fid, struct data_coeff_matrix *cm);

/**
 * Load /pauli_hamil into a packed circ_hamil.
 *
 * Reads dimensions, normalisation, coefficients, and Pauli
 * strings from /pauli_hamil on rank 0, broadcasts them to
 * followers, then allocates the packed term array via
 * circ_hamil_init().  Each term's coefficient is multiplied
 * by the on-disk norm; its Pauli operator is packed into a
 * struct paulis via paulis_set().  Release with
 * circ_hamil_free().
 *
 * Returns 0 on success, -1 on error.
 */
int data_hamil_load(data_id fid, struct circ_hamil *hm);

/*
 * Per-step write API for /circ_{trott,trott2,qdrift,cmpsit}.
 *
 * A writer captures (fid, group, values dataset) so per-step
 * writes do not re-open the dataset each call.  The dset
 * field stores an hid_t for the open /grp/values dataset on
 * rank 0; it is unused on followers.
 *
 * data_circ_writer_init():
 *   Create the named group (idempotent) and pre-allocate the
 *   `values` dataset of shape (n_steps, 2) pre-filled with
 *   NaN, then leave it open on rank 0 with the handle cached
 *   in *w.  Passing fid == 0 zeroes the writer and returns
 *   success: a subsequent data_circ_write_step is a no-op,
 *   so callers can disable per-step output by passing fid=0.
 *   All ranks call.
 *
 * data_circ_write_step():
 *   Hyperslab-write one row of the cached dataset and
 *   H5Fflush.  All ranks call; followers and disabled
 *   writers return 0 without I/O.  A crash between steps
 *   leaves the file consistent up to the last flushed row;
 *   unwritten rows remain NaN.
 *
 * data_circ_writer_close():
 *   Close the cached dataset on rank 0 and zero the writer.
 *   Safe to call on a disabled or already-closed writer.
 *
 * The two write functions return 0 on success, -1 on error
 * (with log_error).
 */
struct data_circ_writer {
	data_id fid;
	int64_t dset;	/* H5Dopen handle on rank 0; 0 otherwise */
	int64_t mspace;	/* (1, 2) memory dataspace cached at init;
			 * reused for every per-step write. */
	size_t n_steps;
};

int data_circ_writer_init(data_id fid, const char *grp_name, size_t n_steps,
	struct data_circ_writer *w);
int data_circ_write_step(struct data_circ_writer *w, size_t step_idx,
	_Complex double z);
void data_circ_writer_close(struct data_circ_writer *w);

#endif // PHASE2_DATA_H
