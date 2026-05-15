#ifndef PHASE2_DATA_H
#define PHASE2_DATA_H

#include <stddef.h>
#include <stdint.h>

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
 * Raw multidet data read from /state_prep/multidet.
 *
 *   nqb   total qubit count of the determinants
 *   ndets number of determinants
 *   cfs   flat double array of shape (ndets, 2), real and
 *         imaginary parts of the complex coefficients
 *   dets  flat byte array of shape (ndets, nqb); each entry
 *         is the occupation of one qubit and is validated
 *         to be 0 or 1 at load time
 *
 * data_multidet_load() allocates cfs and dets and fills the
 * scalars; data_multidet_free() releases them.  The struct
 * remains owned by the caller.
 */
struct data_multidet {
	uint32_t nqb;
	size_t ndets;
	const double *cfs;
	const unsigned char *dets;
};

int data_multidet_load(data_id fid, struct data_multidet *m);
void data_multidet_free(struct data_multidet *m);

/*
 * State-prep dispatch and coefficient-matrix loader.
 *
 * `enum stprep_kind`, `struct data_coeff_block` and
 * `struct data_coeff_matrix` are carrier types consumed by the
 * circuit layer; their definitions live in phase2/circ.h so the
 * algorithm header owns the in-memory representation.  data.h
 * only forward-declares them to keep the I/O surface here while
 * avoiding an include cycle with phase2/state_prep_coeff.h.
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
 * broadcasts to followers in one call.  data_coeff_matrix_free()
 * releases every allocation and zeroes the struct.
 *
 * Documented further in phase2/doc/simul-h5-specs.md
 * "dispatch rules".
 */
enum stprep_kind;
struct data_coeff_matrix;

int data_state_prep_kind(data_id fid, enum stprep_kind *out);
int data_coeff_matrix_load(data_id fid, struct data_coeff_matrix *cm);
void data_coeff_matrix_free(struct data_coeff_matrix *cm);

/**
 * Raw Hamiltonian data read from /pauli_hamil.
 *
 *   nqb     qubit count of the Pauli strings
 *   nterms  number of terms in the Hamiltonian
 *   norm    normalisation factor; coefficients should be
 *           multiplied by this value before use
 *   cfs     flat double array of length nterms
 *   paulis  flat byte array of shape (nterms, nqb), single-
 *           qubit Pauli operators encoded as 0=I, 1=X, 2=Y, 3=Z
 *
 * data_hamil_load() allocates cfs and paulis and fills the
 * scalars; data_hamil_free() releases them.
 */
struct data_hamil {
	uint32_t nqb;
	size_t nterms;
	double norm;
	const double *cfs;
	const unsigned char *paulis;
};

int data_hamil_load(data_id fid, struct data_hamil *h);
void data_hamil_free(struct data_hamil *h);

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
	size_t n_steps;
};

int data_circ_writer_init(data_id fid, const char *grp_name, size_t n_steps,
	struct data_circ_writer *w);
int data_circ_write_step(struct data_circ_writer *w, size_t step_idx,
	_Complex double z);
void data_circ_writer_close(struct data_circ_writer *w);

#endif // PHASE2_DATA_H
