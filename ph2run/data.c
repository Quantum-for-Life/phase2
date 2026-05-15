/*
 * data -- HDF5 I/O layer between simul.h5 and the
 * simulator.
 *
 * Rank-0 ownership.  Rank 0 holds the H5 file handle and
 * performs every H5 call.  Followers receive read results
 * via MPI_Bcast and short-circuit on writes.  The
 * underlying HDF5 build is standard (serial); no
 * parallel-HDF5 driver is needed.
 *
 * Bcast invariant.  On a collective entry point, every
 * rank must reach the *same number* of MPI_Bcast calls
 * even when rank 0 takes a failure path.  The pattern is:
 * rank 0 sets a status int, all ranks BCAST(&rt, 1,
 * MPI_INT), all ranks branch on rt.  A rank-0-only goto
 * that bypasses a bcast would deadlock the followers --
 * which is why each loader has explicit `rt = -1`
 * settings before every potentially-failing call inside
 * the rank-0 block.
 *
 * Buffer ownership.  Load functions allocate through
 * local mutable pointers and assign them to the public
 * carrier struct (struct data_coeff_matrix etc., defined
 * in phase2/circ.h) at the end via const-qualified
 * fields.  The free helpers in phase2/circ.c cast back
 * through void * to release.
 *
 * Logging.  Every error path emits one log_error line
 * naming the failing H5 operation and the relevant
 * context (group, dataset, attribute, index).  Return
 * codes are normalised to {0, -1}.
 */

#define LOG_SUBSYS "data"

#include "c23_compat.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5.h"
#include "mpi.h"

#include "log.h"
#include "phase2/circ.h"
#include "ph2run/data.h"
#include "phase2/paulis.h"
#include "phase2/world.h"
#include <complex.h>

static struct world_info WD;

/*
 * Internal group/dataset/attribute names.  These are
 * implementation details of the on-disk layout; only the
 * load functions and writers in this file reference them,
 * so they live here rather than in the public header.
 */
#define DATA_STPREP "state_prep"
#define DATA_STPREP_MULTIDET "multidet"
#define DATA_STPREP_MULTIDET_COEFFS "coeffs"
#define DATA_STPREP_MULTIDET_DETS "dets"

#define DATA_STPREP_COEFFMAT "coeff_matrix"
#define DATA_STPREP_COEFFMAT_CA "C_alpha"
#define DATA_STPREP_COEFFMAT_CB "C_beta"
#define DATA_STPREP_COEFFMAT_NQB "n_qubits"
#define DATA_STPREP_COEFFMAT_NS "n_sites"
#define DATA_STPREP_COEFFMAT_NA "n_alpha"
#define DATA_STPREP_COEFFMAT_NB "n_beta"
#define DATA_STPREP_COEFFMAT_CS "closed_shell"
#define DATA_STPREP_COEFFMAT_TAP "tapered"
#define DATA_STPREP_COEFFMAT_CSF "csf"
#define DATA_STPREP_COEFFMAT_CSF_NCOMP "n_components"
#define DATA_STPREP_COEFFMAT_CSF_CF "coefficient"

#define DATA_HAMIL "pauli_hamil"
#define DATA_HAMIL_COEFFS "coeffs"
#define DATA_HAMIL_NORM "normalization"
#define DATA_HAMIL_PAULIS "paulis"

/* -- MPI broadcast helper -------------------------------------------- */

/*
 * Rank 0 originates; followers receive.  A zero count
 * short-circuits to a no-op so callers can pass empty
 * arrays without a branch.
 */
#define BCAST(ptr, count, mpi_type)                                            \
	do {                                                                   \
		const int _n = (int)(count);                                   \
		if (_n > 0)                                                    \
			MPI_Bcast((ptr), _n, (mpi_type), 0, MPI_COMM_WORLD);    \
	} while (0)

/* -- common H5 helpers (rank-0-only callers) -------------------------- */

/* Read a whole dataset by name from an already-open group into `buf`. */
static int read_dset(
	hid_t grp_id, const char *name, hid_t native_type, void *buf)
{
	const hid_t did = H5Dopen2(grp_id, name, H5P_DEFAULT);
	if (did == H5I_INVALID_HID) {
		log_error("read_dset: H5Dopen2(%s) failed", name);
		return -1;
	}
	int rt = -1;
	if (H5Dread(did, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0)
		log_error("read_dset: H5Dread(%s) failed", name);
	else
		rt = 0;
	H5Dclose(did);
	return rt;
}

/* Read one scalar attribute by name from an already-open
 * group into the caller-supplied buffer (sized for h5_type). */
static int read_attr_raw(
	hid_t grp, const char *name, hid_t h5_type, void *out)
{
	const hid_t aid = H5Aopen(grp, name, H5P_DEFAULT);
	if (aid == H5I_INVALID_HID) {
		log_error("read_attr(%s): H5Aopen failed", name);
		return -1;
	}
	int rt = -1;
	if (H5Aread(aid, h5_type, out) < 0)
		log_error("read_attr(%s): H5Aread failed", name);
	else
		rt = 0;
	H5Aclose(aid);
	return rt;
}

/* Read the two extents of a 2-D dataset by name from an already-open
 * group.  Sets *out_dims[2] on success. */
static int get_dset_dims2(
	hid_t grp_id, const char *name, hsize_t out_dims[2])
{
	const hid_t dset_id = H5Dopen2(grp_id, name, H5P_DEFAULT);
	if (dset_id == H5I_INVALID_HID) {
		log_error("get_dset_dims2: H5Dopen2(%s) failed", name);
		return -1;
	}
	int rt = -1;
	const hid_t dsp_id = H5Dget_space(dset_id);
	if (dsp_id == H5I_INVALID_HID) {
		log_error("get_dset_dims2(%s): H5Dget_space failed", name);
		goto ex_dset;
	}
	if (H5Sget_simple_extent_dims(dsp_id, out_dims, NULL) != 2) {
		log_error("get_dset_dims2(%s): not a 2-D dataset", name);
		goto ex_space;
	}
	rt = 0;
ex_space:
	H5Sclose(dsp_id);
ex_dset:
	H5Dclose(dset_id);
	return rt;
}

/* -- file open / close ------------------------------------------------ */

/*
 * Open simul.h5 RDWR.  Only rank 0 holds the underlying
 * H5 file handle; followers receive DATA_FOLLOWER_FID
 * and pass it through subsequent collective data_* calls
 * without dereferencing it.  Returns DATA_INVALID_FID on
 * open failure; callers compare against this sentinel.
 */
data_id data_open(const char *filename)
{
	if (world_info(&WD) != WORLD_READY) {
		log_error("data_open(%s): world not ready", filename);
		return DATA_INVALID_FID;
	}

	int status = 0;
	hid_t fid = H5I_INVALID_HID;
	if (WD.rank == 0) {
		const hid_t acc = H5Pcreate(H5P_FILE_ACCESS);
		fid = H5Fopen(filename, H5F_ACC_RDWR, acc);
		H5Pclose(acc);
		if (fid == H5I_INVALID_HID) {
			log_error("data_open(%s): H5Fopen failed", filename);
			status = -1;
		}
	}
	BCAST(&status, 1, MPI_INT);
	if (status < 0)
		return DATA_INVALID_FID;

	if (WD.rank == 0) {
		log_debug("data_open(%s): fid=%lld", filename, (long long)fid);
		return (data_id)fid;
	}
	return DATA_FOLLOWER_FID;
}

/*
 * Close a file id returned by data_open.  Safe to call
 * on DATA_INVALID_FID and DATA_FOLLOWER_FID -- the
 * function inspects both sentinels and skips the
 * H5Fclose on follower ranks.
 */

void data_close(const data_id fid)
{
	if (world_info(&WD) != WORLD_READY)
		return;
	if (WD.rank == 0 && fid != DATA_INVALID_FID
		&& fid != DATA_FOLLOWER_FID) {
		log_debug("data_close: fid=%lld", (long long)fid);
		H5Fclose((hid_t)fid);
	}
}

/* -- group creation --------------------------------------------------- */

/*
 * Idempotence helper: probe `name` under `fid` and
 * return
 *
 *    1 = path-clear, caller must create the group;
 *    0 = a real group already exists, treat as success;
 *   -1 = probe failed.
 *
 * If the path holds a dangling soft link (a legacy
 * fixture artifact), the link is unlinked here and the
 * function returns 1 so the caller can proceed with a
 * fresh create.
 */
static int grp_ensure_clear(hid_t fid, const char *name)
{
	const htri_t lexists = H5Lexists(fid, name, H5P_DEFAULT);
	if (lexists < 0) {
		log_error("data_grp_create(%s): H5Lexists failed", name);
		return -1;
	}
	if (lexists == 0)
		return 1;
	const htri_t oexists = H5Oexists_by_name(fid, name, H5P_DEFAULT);
	if (oexists < 0) {
		log_error("data_grp_create(%s): H5Oexists_by_name failed",
			name);
		return -1;
	}
	if (oexists > 0) {
		log_debug("data_grp_create(%s): already exists", name);
		return 0;
	}
	if (H5Ldelete(fid, name, H5P_DEFAULT) < 0) {
		log_error("data_grp_create(%s): H5Ldelete of dangling link"
			  " failed", name);
		return -1;
	}
	log_debug("data_grp_create(%s): removed dangling link", name);
	return 1;
}

/*
 * Create `name` under `fid` with the standard link-
 * creation property list (intermediate groups, UTF-8
 * link names).
 */
static int grp_make(hid_t fid, const char *name)
{
	int rt = -1;
	const hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
	if (lcpl == H5I_INVALID_HID) {
		log_error("data_grp_create(%s): H5Pcreate failed", name);
		return -1;
	}
	if (H5Pset_create_intermediate_group(lcpl, 1) < 0) {
		log_error("data_grp_create(%s):"
			  " H5Pset_create_intermediate_group failed",
			name);
		goto ex_lcpl;
	}
	if (H5Pset_char_encoding(lcpl, H5T_CSET_UTF8) < 0) {
		log_error("data_grp_create(%s): H5Pset_char_encoding failed",
			name);
		goto ex_lcpl;
	}
	const hid_t grp = H5Gcreate(fid, name, lcpl, H5P_DEFAULT, H5P_DEFAULT);
	if (grp == H5I_INVALID_HID) {
		log_error("data_grp_create(%s): H5Gcreate failed", name);
		goto ex_lcpl;
	}
	H5Gclose(grp);
	rt = 0;
ex_lcpl:
	H5Pclose(lcpl);
	return rt;
}

int data_grp_create(data_id fid, const char *grp_name)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	if (WD.rank == 0) {
		const int clear = grp_ensure_clear((hid_t)fid, grp_name);
		if (clear < 0)
			rt = -1;
		else if (clear > 0 && grp_make((hid_t)fid, grp_name) < 0)
			rt = -1;
	}
	BCAST(&rt, 1, MPI_INT);
	return rt;
}

/* -- attribute read / write ------------------------------------------- */

#define DEFINE_DATA_ATTR_READ(suff, type, h5_type, mpi_type)                   \
	int data_attr_read_##suff(data_id fid, const char *grp_name,           \
		const char *attr_name, type *a)                                \
	{                                                                      \
		if (world_info(&WD) != WORLD_READY)                            \
			return -1;                                             \
		int rt = 0;                                                    \
		type local = (type)0;                                          \
		if (WD.rank == 0) {                                            \
			rt = -1;                                               \
			const hid_t grp_id = H5Gopen(                          \
				(hid_t)fid, grp_name, H5P_DEFAULT);            \
			if (grp_id == H5I_INVALID_HID) {                       \
				log_error("data_attr_read(%s/%s):"             \
					  " H5Gopen failed",                   \
					grp_name, attr_name);                  \
			} else {                                               \
				rt = read_attr_raw(                            \
					grp_id, attr_name, h5_type, &local);   \
				H5Gclose(grp_id);                              \
			}                                                      \
		}                                                              \
		BCAST(&rt, 1, MPI_INT);                                                \
		if (rt < 0)                                                    \
			return -1;                                             \
		MPI_Bcast(&local, 1, mpi_type, 0, MPI_COMM_WORLD);             \
		*a = local;                                                    \
		return 0;                                                      \
	}

DEFINE_DATA_ATTR_READ(dbl, double, H5T_NATIVE_DOUBLE, MPI_DOUBLE);

#define DEFINE_DATA_ATTR_WRITE(suff, type, h5_type)                            \
	int data_attr_write_##suff(data_id fid, const char *grp_name,          \
		const char *attr_name, type a)                                 \
	{                                                                      \
		if (world_info(&WD) != WORLD_READY)                            \
			return -1;                                             \
		int rt = 0;                                                    \
		if (WD.rank == 0) {                                            \
			rt = -1;                                               \
			const hid_t grp_id = H5Gopen(                          \
				(hid_t)fid, grp_name, H5P_DEFAULT);            \
			if (grp_id == H5I_INVALID_HID) {                       \
				log_error("data_attr_write(%s/%s):"            \
					  " H5Gopen failed",                   \
					grp_name, attr_name);                  \
				goto ex_group;                                 \
			}                                                      \
			const hid_t acpl = H5Pcreate(H5P_ATTRIBUTE_CREATE);    \
			if (acpl == H5I_INVALID_HID) {                         \
				log_error("data_attr_write(%s/%s):"            \
					  " H5Pcreate(acpl) failed",           \
					grp_name, attr_name);                  \
				goto ex_acpl;                                  \
			}                                                      \
			if (H5Pset_char_encoding(acpl, H5T_CSET_UTF8) < 0) {   \
				log_error("data_attr_write(%s/%s):"            \
					  " H5Pset_char_encoding failed",      \
					grp_name, attr_name);                  \
				goto ex_fspace;                                \
			}                                                      \
			const hid_t fspace = H5Screate(H5S_SCALAR);            \
			if (fspace == H5I_INVALID_HID) {                       \
				log_error("data_attr_write(%s/%s):"            \
					  " H5Screate failed",                 \
					grp_name, attr_name);                  \
				goto ex_fspace;                                \
			}                                                      \
			const hid_t attr_id = H5Acreate2(grp_id, attr_name,    \
				h5_type, fspace, acpl, H5P_DEFAULT);           \
			if (attr_id == H5I_INVALID_HID) {                      \
				log_error("data_attr_write(%s/%s):"            \
					  " H5Acreate2 failed",                \
					grp_name, attr_name);                  \
				goto ex_attr;                                  \
			}                                                      \
			if (H5Awrite(attr_id, h5_type, &a) < 0) {              \
				log_error("data_attr_write(%s/%s):"            \
					  " H5Awrite failed",                  \
					grp_name, attr_name);                  \
				goto ex_write;                                 \
			}                                                      \
			rt = 0;                                                \
		ex_write:                                                      \
			H5Aclose(attr_id);                                     \
		ex_attr:                                                       \
			H5Sclose(fspace);                                      \
		ex_fspace:                                                     \
			H5Pclose(acpl);                                        \
		ex_acpl:                                                       \
			H5Gclose(grp_id);                                      \
		ex_group:;                                                     \
		}                                                              \
		BCAST(&rt, 1, MPI_INT);                                                \
		return rt;                                                     \
	}

DEFINE_DATA_ATTR_WRITE(ul, unsigned long, H5T_NATIVE_ULONG);
DEFINE_DATA_ATTR_WRITE(dbl, double, H5T_NATIVE_DOUBLE);

/* -- multidet --------------------------------------------------------- */

#define MULTIDET_PATH DATA_STPREP "/" DATA_STPREP_MULTIDET

/*
 * Load /state_prep/multidet into a packed struct
 * circ_muldet (idx, cf) view.
 *
 * Phases:
 *  1. rank 0 opens the group and reads the dimensions
 *     of the dets dataset; ndets and nqb broadcast.
 *  2. every rank malloc's the raw cfs (2*ndets doubles)
 *     and dets (ndets*nqb bytes) buffers, since the
 *     bcast must land into already-allocated memory.
 *  3. rank 0 reads the two datasets; the group handle
 *     stays open from step 1 to avoid a redundant
 *     H5Gopen.
 *  4. every rank receives the bcast, validates that
 *     each occupation byte is in {0, 1} (rejecting a
 *     malformed multidet group with -1 + log_error),
 *     then packs into the {idx, cf} form circ
 *     consumes.
 */
int data_muldet_load(data_id fid, struct circ_muldet *md)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	memset(md, 0, sizeof *md);

	int rt = 0;
	uint32_t v_nqb = 0;
	size_t v_ndets = 0;
	hid_t grp_id = H5I_INVALID_HID;
	if (WD.rank == 0) {
		rt = -1;
		grp_id = H5Gopen((hid_t)fid, MULTIDET_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_muldet_load: H5Gopen(%s) failed",
				MULTIDET_PATH);
			goto ex_dims;
		}
		hsize_t dims[2];
		if (get_dset_dims2(grp_id, DATA_STPREP_MULTIDET_DETS, dims)
			< 0)
			goto ex_dims;
		v_ndets = dims[0];
		v_nqb = (uint32_t)dims[1];
		rt = 0;
	ex_dims:;
	}
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0) {
		if (WD.rank == 0 && grp_id != H5I_INVALID_HID)
			H5Gclose(grp_id);
		return -1;
	}
	BCAST(&v_nqb, 1, MPI_UINT32_T);
	BCAST(&v_ndets, 1, MPI_UNSIGNED_LONG);

	double *cfs = malloc(sizeof *cfs * 2 * v_ndets);
	unsigned char *dets = malloc(sizeof *dets * v_ndets * v_nqb);
	if (!cfs || !dets) {
		log_error("data_muldet_load: alloc failed"
			  " (ndets=%zu, nqb=%u)", v_ndets, v_nqb);
		if (WD.rank == 0)
			H5Gclose(grp_id);
		free(cfs);
		free(dets);
		return -1;
	}

	if (WD.rank == 0) {
		if (read_dset(grp_id, DATA_STPREP_MULTIDET_COEFFS,
			    H5T_NATIVE_DOUBLE, cfs) < 0
			|| read_dset(grp_id, DATA_STPREP_MULTIDET_DETS,
				   H5T_NATIVE_UCHAR, dets) < 0)
			rt = -1;
		H5Gclose(grp_id);
	}
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0) {
		free(cfs);
		free(dets);
		return -1;
	}
	BCAST(cfs, 2 * v_ndets, MPI_DOUBLE);
	BCAST(dets, v_ndets * v_nqb, MPI_UNSIGNED_CHAR);

	/* Bit-validate the determinant occupations on every rank
	 * (cheap, ndets*nqb bytes) so a malformed multidet group
	 * fails fast instead of corrupting downstream indices. */
	for (size_t i = 0; i < v_ndets; i++) {
		for (size_t j = 0; j < v_nqb; j++) {
			const unsigned char bit = dets[i * v_nqb + j];
			if (bit > 1) {
				log_error("data_muldet_load: dets[%zu][%zu]"
					  " = %u is not 0/1; multidet group"
					  " is malformed", i, j,
					(unsigned)bit);
				free(cfs);
				free(dets);
				return -1;
			}
		}
	}

	if (circ_muldet_init(md, v_ndets) < 0) {
		log_error("data_muldet_load: circ_muldet_init failed"
			  " (ndets=%zu)", v_ndets);
		free(cfs);
		free(dets);
		return -1;
	}
	for (size_t i = 0; i < v_ndets; i++) {
		md->dets[i].cf = CMPLX(cfs[2 * i], cfs[2 * i + 1]);
		uint64_t idx = 0;
		for (uint32_t j = 0; j < v_nqb; j++)
			idx |= (uint64_t)dets[i * v_nqb + j] << j;
		md->dets[i].idx = idx;
	}
	free(cfs);
	free(dets);
	return 0;
}

/* -- pauli_hamil ------------------------------------------------------ */

/*
 * Load /pauli_hamil into a packed struct circ_hamil
 * (cf, struct paulis) view ready for circ_step.
 *
 * Phases mirror data_muldet_load: rank-0 dim read +
 * bcast, all-ranks malloc, rank-0 dataset read + bcast,
 * all-ranks packing.  Packing scales each coefficient by
 * the on-disk normalisation attribute and converts the
 * per-qubit byte array (0=I, 1=X, 2=Y, 3=Z) into the
 * struct paulis bitstring representation via
 * paulis_set().
 */
int data_hamil_load(data_id fid, struct circ_hamil *hm)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	uint32_t v_nqb = 0;
	size_t v_nterms = 0;
	double v_norm = 0.0;
	hid_t grp_id = H5I_INVALID_HID;
	if (WD.rank == 0) {
		rt = -1;
		grp_id = H5Gopen((hid_t)fid, DATA_HAMIL, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_hamil_load: H5Gopen(%s) failed",
				DATA_HAMIL);
			goto ex_dims;
		}
		hsize_t dims[2];
		if (get_dset_dims2(grp_id, DATA_HAMIL_PAULIS, dims) < 0
			|| read_attr_raw(grp_id, DATA_HAMIL_NORM,
				   H5T_NATIVE_DOUBLE, &v_norm) < 0)
			goto ex_dims;
		v_nterms = dims[0];
		v_nqb = (uint32_t)dims[1];
		rt = 0;
	ex_dims:;
	}
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0) {
		if (WD.rank == 0 && grp_id != H5I_INVALID_HID)
			H5Gclose(grp_id);
		return -1;
	}
	BCAST(&v_nqb, 1, MPI_UINT32_T);
	BCAST(&v_nterms, 1, MPI_UNSIGNED_LONG);
	BCAST(&v_norm, 1, MPI_DOUBLE);

	double *cfs = malloc(sizeof *cfs * v_nterms);
	unsigned char *paulis = malloc(sizeof *paulis * v_nterms * v_nqb);
	if (!cfs || !paulis) {
		log_error("data_hamil_load: alloc failed"
			  " (nterms=%zu, nqb=%u)", v_nterms, v_nqb);
		if (WD.rank == 0)
			H5Gclose(grp_id);
		free(cfs);
		free(paulis);
		return -1;
	}

	if (WD.rank == 0) {
		if (read_dset(grp_id, DATA_HAMIL_COEFFS,
			    H5T_NATIVE_DOUBLE, cfs) < 0
			|| read_dset(grp_id, DATA_HAMIL_PAULIS,
				   H5T_NATIVE_UCHAR, paulis) < 0)
			rt = -1;
		H5Gclose(grp_id);
	}
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0) {
		free(cfs);
		free(paulis);
		return -1;
	}
	BCAST(cfs, v_nterms, MPI_DOUBLE);
	BCAST(paulis, v_nterms * v_nqb, MPI_UNSIGNED_CHAR);

	/* Validate the per-qubit Pauli operator bytes on every
	 * rank (cheap, nterms*nqb bytes).  Without this a
	 * corrupt simul.h5 with bytes > 3 would silently
	 * truncate to one of {I, X, Y, Z} inside paulis_set
	 * and the simulator would produce nonsense results
	 * with no diagnostic. */
	for (size_t i = 0; i < v_nterms; i++) {
		for (uint32_t j = 0; j < v_nqb; j++) {
			const unsigned char b = paulis[i * v_nqb + j];
			if (b > 3) {
				log_error("data_hamil_load: paulis[%zu][%u] ="
					  " %u is not I/X/Y/Z; pauli_hamil"
					  " group is malformed",
					i, j, (unsigned)b);
				free(cfs);
				free(paulis);
				return -1;
			}
		}
	}

	if (circ_hamil_init(hm, v_nqb, v_nterms) < 0) {
		free(cfs);
		free(paulis);
		return -1;
	}
	for (size_t i = 0; i < v_nterms; i++) {
		hm->terms[i].cf = cfs[i] * v_norm;
		struct paulis op = paulis_new();
		for (uint32_t j = 0; j < v_nqb; j++)
			paulis_set(&op, paulis[i * v_nqb + j], j);
		hm->terms[i].op = op;
	}
	free(cfs);
	free(paulis);
	return 0;
}

/* -- per-step write API (/circ_{trott,trott2,qdrift,cmpsit}/values) --- */

#define CIRC_VALUES_DSET "values"

/*
 * Open an existing /grp/values dataset whose shape has
 * already been chosen by an earlier ph2run invocation
 * on the same simul.h5.  Caller owns the close.
 * Returns H5I_INVALID_HID on failure.
 */
static hid_t values_reuse(hid_t grp_id, const char *grp_name)
{
	const hid_t dset = H5Dopen2(grp_id, CIRC_VALUES_DSET, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID) {
		log_error("data_circ_writer_init(%s): H5Dopen2(values) failed",
			grp_name);
		return H5I_INVALID_HID;
	}
	log_debug("data_circ_writer_init(%s): reusing existing values"
		  " dataset", grp_name);
	return dset;
}

/*
 * Create a new /grp/values dataset of shape (n_steps,
 * 2) pre-filled with NaN.  The file is flushed before
 * return so a crash mid-simulation always observes
 * either a real row (the last write_step) or the
 * NaN-padded tail.
 */
static hid_t values_create(hid_t fid, hid_t grp_id, const char *grp_name,
	size_t n_steps)
{
	hid_t dset = H5I_INVALID_HID;
	const hid_t dspace = H5Screate_simple(2,
		(hsize_t[]){ n_steps, 2 }, NULL);
	if (dspace == H5I_INVALID_HID) {
		log_error("data_circ_writer_init(%s): H5Screate_simple"
			  " failed (n_steps=%zu)", grp_name, n_steps);
		return H5I_INVALID_HID;
	}
	const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
	if (dcpl == H5I_INVALID_HID) {
		log_error("data_circ_writer_init(%s): H5Pcreate(dcpl) failed",
			grp_name);
		goto ex_dspace;
	}
	const double nan_val = nan("");
	if (H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &nan_val) < 0) {
		log_error("data_circ_writer_init(%s): H5Pset_fill_value"
			  " failed", grp_name);
		goto ex_dcpl;
	}
	if (H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY) < 0) {
		log_error("data_circ_writer_init(%s): H5Pset_alloc_time"
			  " failed", grp_name);
		goto ex_dcpl;
	}
	if (H5Pset_fill_time(dcpl, H5D_FILL_TIME_ALLOC) < 0) {
		log_error("data_circ_writer_init(%s): H5Pset_fill_time"
			  " failed", grp_name);
		goto ex_dcpl;
	}
	dset = H5Dcreate2(grp_id, CIRC_VALUES_DSET, H5T_IEEE_F64LE, dspace,
		H5P_DEFAULT, dcpl, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID) {
		log_error("data_circ_writer_init(%s): H5Dcreate2(values)"
			  " failed", grp_name);
		goto ex_dcpl;
	}
	H5Fflush(fid, H5F_SCOPE_GLOBAL);
	log_debug("data_circ_writer_init(%s): values shape (%zu, 2),"
		  " NaN-padded", grp_name, n_steps);
ex_dcpl:
	H5Pclose(dcpl);
ex_dspace:
	H5Sclose(dspace);
	return dset;
}

/*
 * Initialise a writer at /grp/values.  The constructor
 * handles both first-run (create + NaN-pad) and re-run
 * (reuse the existing dataset) paths.  Two collective
 * BCAST calls bracket the rank-0 work so a follower
 * never returns success while rank 0 saw an error.
 *
 * fid == 0 is the documented "no per-step writes"
 * sentinel: the writer is zeroed and every subsequent
 * write_step / writer_close returns without I/O.  This
 * lets a caller drive trott_simul / qdrift_simul / etc.
 * with output disabled (in-memory smoke tests).
 */
int data_circ_writer_init(data_id fid, const char *grp_name, size_t n_steps,
	struct data_circ_writer *w)
{
	memset(w, 0, sizeof *w);
	if (fid == 0)
		return 0;
	if (world_info(&WD) != WORLD_READY)
		return -1;

	/* data_grp_create runs the bcast internally. */
	if (data_grp_create(fid, grp_name) < 0)
		return -1;

	int rt = 0;
	hid_t dset_out = H5I_INVALID_HID;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, grp_name, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_circ_writer_init(%s): H5Gopen failed",
				grp_name);
			goto ex_bcast;
		}
		const htri_t exists = H5Lexists(grp_id, CIRC_VALUES_DSET,
			H5P_DEFAULT);
		if (exists < 0)
			log_error("data_circ_writer_init(%s):"
				  " H5Lexists(values) failed", grp_name);
		else if (exists > 0)
			dset_out = values_reuse(grp_id, grp_name);
		else
			dset_out = values_create((hid_t)fid, grp_id, grp_name,
				n_steps);
		if (dset_out != H5I_INVALID_HID)
			rt = 0;
		H5Gclose(grp_id);
	ex_bcast:;
	}
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0) {
		if (WD.rank == 0 && dset_out != H5I_INVALID_HID)
			H5Dclose(dset_out);
		return -1;
	}
	hid_t mspace_out = H5I_INVALID_HID;
	if (WD.rank == 0) {
		mspace_out = H5Screate_simple(2,
			(hsize_t[]){ 1, 2 }, NULL);
		if (mspace_out == H5I_INVALID_HID) {
			log_error("data_circ_writer_init(%s):"
				  " H5Screate_simple(mspace) failed",
				grp_name);
			H5Dclose(dset_out);
			rt = -1;
		}
	}
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0)
		return -1;
	w->fid = fid;
	w->n_steps = n_steps;
	w->dset = (WD.rank == 0) ? (int64_t)dset_out : 0;
	w->mspace = (WD.rank == 0) ? (int64_t)mspace_out : 0;
	return 0;
}

/*
 * Write one row of the cached /grp/values dataset and
 * H5Fflush.  The flush is the per-step atomicity
 * guarantee: a crash between steps leaves rows
 * 0..step_idx-1 on disk and rows step_idx..n_steps-1
 * as NaN, so a SLURM timeout produces a useful partial
 * log instead of a corrupt file.  Followers and
 * disabled writers short-circuit without I/O.
 */
int data_circ_write_step(struct data_circ_writer *w, size_t step_idx,
	_Complex double z)
{
	if (w->fid == 0)
		return 0;
	if (world_info(&WD) != WORLD_READY)
		return -1;
	if (WD.rank != 0)
		return 0;

	int rt = -1;
	const hid_t dset = (hid_t)w->dset;
	const hid_t mspace = (hid_t)w->mspace;
	const hid_t fspace = H5Dget_space(dset);
	if (fspace == H5I_INVALID_HID) {
		log_error("data_circ_write_step(%zu): H5Dget_space failed",
			step_idx);
		return -1;
	}
	const hsize_t start[2] = { step_idx, 0 };
	const hsize_t count[2] = { 1, 2 };
	if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count,
		    NULL) < 0) {
		log_error("data_circ_write_step(%zu): H5Sselect_hyperslab"
			  " failed", step_idx);
		goto ex_fspace;
	}
	const double row[2] = { creal(z), cimag(z) };
	if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, row)
		< 0) {
		log_error("data_circ_write_step(%zu): H5Dwrite failed",
			step_idx);
		goto ex_fspace;
	}
	/* Atomic-on-disk: rows 0..step_idx are persisted before
	 * the simulation moves on. */
	H5Fflush((hid_t)w->fid, H5F_SCOPE_GLOBAL);

	rt = 0;
ex_fspace:
	H5Sclose(fspace);
	return rt;
}

void data_circ_writer_close(struct data_circ_writer *w)
{
	if (!w || w->fid == 0)
		return;
	if (world_info(&WD) != WORLD_READY)
		return;
	if (WD.rank == 0) {
		if (w->dset > 0)
			H5Dclose((hid_t)w->dset);
		if (w->mspace > 0)
			H5Sclose((hid_t)w->mspace);
	}
	memset(w, 0, sizeof *w);
}

/* -- state-prep dispatch ---------------------------------------------- */

/*
 * Probe simul.h5 to determine which state-prep
 * subgroup the file carries.  Both-present and
 * neither-present are rejected with a log_error line;
 * the caller (ph2run/circ) is expected to fail-fast on
 * the misbuilt file rather than guess which payload to
 * load.  See doc/simul-h5-specs.md "dispatch rules".
 */
int data_state_prep_kind(const data_id fid, enum stprep_kind *out)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	int kind = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t sp_id = H5Gopen(
			(hid_t)fid, DATA_STPREP, H5P_DEFAULT);
		if (sp_id == H5I_INVALID_HID) {
			log_error("data_state_prep_kind: %s group missing",
				DATA_STPREP);
			goto ex_bcast;
		}

		const htri_t has_md =
			H5Lexists(sp_id, DATA_STPREP_MULTIDET, H5P_DEFAULT);
		const htri_t has_cm =
			H5Lexists(sp_id, DATA_STPREP_COEFFMAT, H5P_DEFAULT);
		H5Gclose(sp_id);

		if (has_md < 0 || has_cm < 0) {
			log_error("data_state_prep_kind: H5Lexists failed"
				  " (has_md=%d has_cm=%d)",
				(int)has_md, (int)has_cm);
			goto ex_bcast;
		}
		if (has_md && has_cm) {
			log_error("simul.h5: ambiguous state prep (both"
				  " /state_prep/multidet and"
				  " /state_prep/coeff_matrix present);"
				  " rebuild simul.h5 with exactly one");
			goto ex_bcast;
		}
		if (!has_md && !has_cm) {
			log_error("simul.h5: no state-prep subgroup found"
				  " (expected /state_prep/multidet or"
				  " /state_prep/coeff_matrix)");
			goto ex_bcast;
		}

		kind = has_md ? STPREP_MULTIDET : STPREP_COEFF_MATRIX;
		log_debug("data_state_prep_kind: %s",
			(kind == STPREP_MULTIDET) ? "multidet"
						  : "coeff_matrix");
		rt = 0;
	ex_bcast:;
	}
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0)
		return -1;
	BCAST(&kind, 1, MPI_INT);
	*out = (enum stprep_kind)kind;
	return 0;
}

/* -- coeff_matrix ----------------------------------------------------- */

#define COEFFMAT_PATH DATA_STPREP "/" DATA_STPREP_COEFFMAT

static int read_coeff_attrs(hid_t grp_id, struct data_coeff_matrix *cm)
{
	uint8_t cs = 0, tap = 0;
	if (read_attr_raw(grp_id, DATA_STPREP_COEFFMAT_NQB,
		    H5T_NATIVE_UINT32, &cm->nqb) < 0
		|| read_attr_raw(grp_id, DATA_STPREP_COEFFMAT_NS,
			   H5T_NATIVE_UINT32, &cm->n_sites) < 0
		|| read_attr_raw(grp_id, DATA_STPREP_COEFFMAT_NA,
			   H5T_NATIVE_UINT32, &cm->n_alpha) < 0
		|| read_attr_raw(grp_id, DATA_STPREP_COEFFMAT_NB,
			   H5T_NATIVE_UINT32, &cm->n_beta) < 0
		|| read_attr_raw(grp_id, DATA_STPREP_COEFFMAT_CS,
			   H5T_NATIVE_UINT8, &cs) < 0
		|| read_attr_raw(grp_id, DATA_STPREP_COEFFMAT_TAP,
			   H5T_NATIVE_UINT8, &tap) < 0)
		return -1;
	cm->closed_shell = cs ? 1 : 0;
	cm->tapered = tap ? 1 : 0;
	return 0;
}

static int read_C_dset(hid_t grp_id, const char *name, uint32_t n_sites,
	uint32_t n_occ, double *buf)
{
	const hid_t did = H5Dopen2(grp_id, name, H5P_DEFAULT);
	if (did == H5I_INVALID_HID) {
		log_error("read_C_dset: H5Dopen2(%s) failed", name);
		return -1;
	}
	int rt = -1;

	const hid_t sid = H5Dget_space(did);
	if (sid == H5I_INVALID_HID) {
		log_error("read_C_dset(%s): H5Dget_space failed", name);
		goto ex_dset;
	}
	hsize_t dims[2] = { 0, 0 };
	if (H5Sget_simple_extent_dims(sid, dims, NULL) != 2) {
		log_error("read_C_dset(%s): dataset is not 2-D", name);
		goto ex_space;
	}
	if (dims[0] != n_sites || dims[1] != n_occ) {
		log_error("read_C_dset(%s): shape (%llu,%llu) does not match"
			  " (%u,%u)",
			name, (unsigned long long)dims[0],
			(unsigned long long)dims[1], n_sites, n_occ);
		goto ex_space;
	}
	const hid_t tid = H5Dget_type(did);
	if (tid == H5I_INVALID_HID) {
		log_error("read_C_dset(%s): H5Dget_type failed", name);
		goto ex_space;
	}
	const H5T_class_t cls = H5Tget_class(tid);
	const size_t sz = H5Tget_size(tid);
	H5Tclose(tid);
	if (cls != H5T_FLOAT || sz != sizeof(double)) {
		log_error("read_C_dset(%s): expected double, got class=%d"
			  " size=%zu",
			name, (int)cls, sz);
		goto ex_space;
	}

	if (H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf)
		< 0) {
		log_error("read_C_dset: H5Dread(%s) failed", name);
		goto ex_space;
	}
	rt = 0;
ex_space:
	H5Sclose(sid);
ex_dset:
	H5Dclose(did);
	return rt;
}

/*
 * rank-0 metadata read: opens /state_prep/coeff_matrix,
 * fills the scalar fields on cm, and reports whether a
 * csf/ subgroup is present plus its declared component
 * count.  An explicit csf/ with n_components == 0 is
 * non-physical (a single-block file must omit the
 * subgroup) and is rejected here.
 */
static int read_coeff_meta(hid_t fid, struct data_coeff_matrix *cm,
	int *has_csf, uint32_t *ncomp)
{
	*has_csf = 0;
	*ncomp = 0;
	const hid_t grp_id = H5Gopen(fid, COEFFMAT_PATH, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		log_error("data_coeff_matrix_load: H5Gopen(%s) failed",
			COEFFMAT_PATH);
		return -1;
	}
	int rt = -1;
	if (read_coeff_attrs(grp_id, cm) < 0)
		goto ex_grp;
	const htri_t has = H5Lexists(
		grp_id, DATA_STPREP_COEFFMAT_CSF, H5P_DEFAULT);
	if (has < 0) {
		log_error("data_coeff_matrix_load: H5Lexists(csf) failed");
		goto ex_grp;
	}
	if (has == 0) {
		rt = 0;
		goto ex_grp;
	}
	const hid_t cg = H5Gopen(
		grp_id, DATA_STPREP_COEFFMAT_CSF, H5P_DEFAULT);
	if (cg == H5I_INVALID_HID) {
		log_error("data_coeff_matrix_load: H5Gopen(csf) failed");
		goto ex_grp;
	}
	const int rt_n = read_attr_raw(cg, DATA_STPREP_COEFFMAT_CSF_NCOMP,
		H5T_NATIVE_UINT32, ncomp);
	H5Gclose(cg);
	if (rt_n < 0)
		goto ex_grp;
	if (*ncomp == 0) {
		log_error("simul.h5: /state_prep/coeff_matrix/csf present"
			  " with n_components=0; remove the csf subgroup"
			  " or list at least one component");
		goto ex_grp;
	}
	*has_csf = 1;
	rt = 0;
ex_grp:
	H5Gclose(grp_id);
	return rt;
}

/*
 * All-ranks bcast of the metadata read on rank 0.
 * Followers receive zeroed values until this completes.
 */
static void bcast_coeff_meta(struct data_coeff_matrix *cm,
	int *has_csf, uint32_t *ncomp)
{
	BCAST(&cm->nqb, 1, MPI_UINT32_T);
	BCAST(&cm->n_sites, 1, MPI_UINT32_T);
	BCAST(&cm->n_alpha, 1, MPI_UINT32_T);
	BCAST(&cm->n_beta, 1, MPI_UINT32_T);
	BCAST(&cm->closed_shell, 1, MPI_INT);
	BCAST(&cm->tapered, 1, MPI_INT);
	BCAST(has_csf, 1, MPI_INT);
	BCAST(ncomp, 1, MPI_UINT32_T);
}

/*
 * Release helper for the error path.  Mirrors the
 * allocator below: frees the top-level Ca/Cb of a
 * single-block layout OR the per-block buffers + the
 * blocks array of a CSF layout, then zeros the struct
 * so the caller is never left with dangling pointers.
 */
static void free_coeff_buffers(uint32_t ncomp, double *Ca, double *Cb,
	struct data_coeff_block *blocks)
{
	free(Ca);
	free(Cb);
	if (blocks) {
		for (size_t k = 0; k < ncomp; k++) {
			free((void *)blocks[k].C_alpha);
			free((void *)blocks[k].C_beta);
		}
		free(blocks);
	}
}

/*
 * All-ranks allocator.  On success returns 0 with the
 * out-pointers populated.  On failure frees any partial
 * allocation and returns -1.
 *
 * The single-block layout uses *Ca / *Cb; the CSF layout
 * uses *blocks with per-block C_alpha / C_beta.
 */
static int alloc_coeff_buffers(const struct data_coeff_matrix *cm,
	int has_csf, uint32_t ncomp,
	double **Ca, double **Cb, struct data_coeff_block **blocks)
{
	const size_t sz_a = (size_t)cm->n_sites * cm->n_alpha;
	const size_t sz_b = (size_t)cm->n_sites * cm->n_beta;
	*Ca = NULL;
	*Cb = NULL;
	*blocks = NULL;
	if (!has_csf) {
		*Ca = malloc(sizeof(double) * (sz_a ? sz_a : 1));
		if (!*Ca) {
			log_error("data_coeff_matrix_load: alloc C_alpha"
				  " failed");
			return -1;
		}
		if (!cm->closed_shell) {
			*Cb = malloc(sizeof(double) * (sz_b ? sz_b : 1));
			if (!*Cb) {
				log_error("data_coeff_matrix_load: alloc"
					  " C_beta failed");
				free_coeff_buffers(0, *Ca, NULL, NULL);
				*Ca = NULL;
				return -1;
			}
		}
		return 0;
	}
	*blocks = calloc(ncomp, sizeof **blocks);
	if (!*blocks) {
		log_error("data_coeff_matrix_load: alloc blocks (n=%u)"
			  " failed", ncomp);
		return -1;
	}
	for (size_t k = 0; k < ncomp; k++) {
		double *bca = malloc(sizeof(double) * (sz_a ? sz_a : 1));
		if (!bca) {
			log_error("data_coeff_matrix_load: alloc C_alpha[%zu]"
				  " failed", k);
			free_coeff_buffers(ncomp, NULL, NULL, *blocks);
			*blocks = NULL;
			return -1;
		}
		(*blocks)[k].C_alpha = bca;
		if (!cm->closed_shell) {
			double *bcb = malloc(sizeof(double)
				* (sz_b ? sz_b : 1));
			if (!bcb) {
				log_error("data_coeff_matrix_load: alloc"
					  " C_beta[%zu] failed", k);
				free_coeff_buffers(ncomp, NULL, NULL, *blocks);
				*blocks = NULL;
				return -1;
			}
			(*blocks)[k].C_beta = bcb;
		}
	}
	return 0;
}

/*
 * rank-0 dataset read into the pre-allocated buffers.
 * In the CSF case opens csf/<k>/ once per component and
 * reads the per-block weight + coefficient datasets.
 */
static int read_coeff_dsets(hid_t fid, const struct data_coeff_matrix *cm,
	int has_csf, uint32_t ncomp,
	double *Ca, double *Cb, struct data_coeff_block *blocks)
{
	const hid_t grp_id = H5Gopen(fid, COEFFMAT_PATH, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		log_error("data_coeff_matrix_load: H5Gopen(%s) failed (dsets)",
			COEFFMAT_PATH);
		return -1;
	}
	int rt = -1;
	if (!has_csf) {
		if (read_C_dset(grp_id, DATA_STPREP_COEFFMAT_CA,
			    cm->n_sites, cm->n_alpha, Ca) < 0)
			goto ex_grp;
		if (!cm->closed_shell
			&& read_C_dset(grp_id, DATA_STPREP_COEFFMAT_CB,
				   cm->n_sites, cm->n_beta, Cb) < 0)
			goto ex_grp;
		rt = 0;
		goto ex_grp;
	}
	for (size_t k = 0; k < ncomp; k++) {
		char sub[32];
		snprintf(sub, sizeof sub, "%s/%zu",
			DATA_STPREP_COEFFMAT_CSF, k);
		const hid_t cg = H5Gopen(grp_id, sub, H5P_DEFAULT);
		if (cg == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_load: H5Gopen(%s/%s)"
				  " failed", COEFFMAT_PATH, sub);
			goto ex_grp;
		}
		const int ok = read_attr_raw(cg, DATA_STPREP_COEFFMAT_CSF_CF,
				    H5T_NATIVE_DOUBLE, &blocks[k].cf) >= 0
			&& read_C_dset(cg, DATA_STPREP_COEFFMAT_CA,
				   cm->n_sites, cm->n_alpha,
				   (double *)blocks[k].C_alpha) >= 0
			&& (cm->closed_shell
				|| read_C_dset(cg, DATA_STPREP_COEFFMAT_CB,
					   cm->n_sites, cm->n_beta,
					   (double *)blocks[k].C_beta)
					>= 0);
		H5Gclose(cg);
		if (!ok)
			goto ex_grp;
	}
	rt = 0;
ex_grp:
	H5Gclose(grp_id);
	return rt;
}

/*
 * All-ranks bcast of the coefficient datasets read on
 * rank 0.  Mirrors the layout chosen by
 * alloc_coeff_buffers.
 */
static void bcast_coeff_dsets(const struct data_coeff_matrix *cm,
	int has_csf, uint32_t ncomp,
	double *Ca, double *Cb, struct data_coeff_block *blocks)
{
	const size_t sz_a = (size_t)cm->n_sites * cm->n_alpha;
	const size_t sz_b = (size_t)cm->n_sites * cm->n_beta;
	if (!has_csf) {
		BCAST(Ca, sz_a, MPI_DOUBLE);
		if (!cm->closed_shell)
			BCAST(Cb, sz_b, MPI_DOUBLE);
		return;
	}
	for (size_t k = 0; k < ncomp; k++) {
		BCAST(&blocks[k].cf, 1, MPI_DOUBLE);
		BCAST((double *)blocks[k].C_alpha, sz_a, MPI_DOUBLE);
		if (!cm->closed_shell)
			BCAST((double *)blocks[k].C_beta, sz_b, MPI_DOUBLE);
	}
}

/*
 * Load /state_prep/coeff_matrix into a struct
 * data_coeff_matrix.
 *
 * The five-phase orchestration:
 *
 *   1. read_coeff_meta   -- rank-0 attribute reads
 *      determine nqb / n_sites / n_alpha / n_beta /
 *      closed_shell / tapered, and probe whether a
 *      CSF subgroup is present plus its arity.
 *   2. bcast_coeff_meta  -- followers receive those
 *      scalars so they can size the same allocations
 *      rank 0 will fill.
 *   3. alloc_coeff_buffers -- all ranks allocate
 *      either a single-block (Ca, Cb) layout or a
 *      CSF blocks[] array.
 *   4. read_coeff_dsets  -- rank 0 reads the C arrays
 *      into the pre-allocated buffers.
 *   5. bcast_coeff_dsets -- followers receive the C
 *      arrays.
 *
 * Each phase has its own collective failure point.  If
 * any rank fails, every rank frees what it allocated
 * and returns the struct to its zero-initialised state.
 */
int data_coeff_matrix_load(const data_id fid, struct data_coeff_matrix *cm)
{
	memset(cm, 0, sizeof *cm);
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	int has_csf = 0;
	uint32_t ncomp = 0;
	if (WD.rank == 0
		&& read_coeff_meta((hid_t)fid, cm, &has_csf, &ncomp) < 0)
		rt = -1;
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0) {
		memset(cm, 0, sizeof *cm);
		return -1;
	}
	bcast_coeff_meta(cm, &has_csf, &ncomp);

	double *Ca = NULL, *Cb = NULL;
	struct data_coeff_block *blocks = NULL;
	if (alloc_coeff_buffers(cm, has_csf, ncomp, &Ca, &Cb, &blocks) < 0) {
		memset(cm, 0, sizeof *cm);
		return -1;
	}

	rt = 0;
	if (WD.rank == 0
		&& read_coeff_dsets((hid_t)fid, cm, has_csf, ncomp,
			   Ca, Cb, blocks) < 0)
		rt = -1;
	BCAST(&rt, 1, MPI_INT);
	if (rt < 0) {
		free_coeff_buffers(ncomp, Ca, Cb, blocks);
		memset(cm, 0, sizeof *cm);
		return -1;
	}
	bcast_coeff_dsets(cm, has_csf, ncomp, Ca, Cb, blocks);

	/* Cast through const: the struct fields are
	 * const-qualified to the public consumer, but data.c
	 * owns the lifetime via data_coeff_matrix_free in
	 * phase2/circ.c (which casts back to void * to call
	 * free). */
	cm->C_alpha = Ca;
	cm->C_beta = Cb;
	cm->blocks = blocks;
	cm->n_components = has_csf ? ncomp : 0;
	return 0;
}

