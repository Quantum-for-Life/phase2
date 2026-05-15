/*
 * data -- HDF5 I/O layer between simul.h5 and the simulator.
 *
 * Rank 0 owns the file handle and performs every H5 call.
 * Other ranks receive read results via MPI_Bcast and short-
 * circuit on writes.  The underlying HDF5 build is standard
 * (serial); no parallel-HDF5 driver is needed.
 *
 * Every error path emits a log_error line naming the failing
 * H5 operation and the relevant context (group, dataset,
 * attribute, index).  Return codes are normalised to {0, -1}.
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
#include "phase2/data.h"
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

/* -- MPI broadcast helpers -------------------------------------------- */

static inline void bcast_int(int *v)
{
	MPI_Bcast(v, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

static inline void bcast_u32(uint32_t *v)
{
	MPI_Bcast(v, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
}

static inline void bcast_sz(size_t *v)
{
	MPI_Bcast(v, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
}

static inline void bcast_doubles(double *buf, size_t n)
{
	if (n > 0)
		MPI_Bcast(buf, (int)n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

static inline void bcast_uchars(unsigned char *buf, size_t n)
{
	if (n > 0)
		MPI_Bcast(buf, (int)n, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
}

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
	bcast_int(&status);
	if (status < 0)
		return DATA_INVALID_FID;

	if (WD.rank == 0) {
		log_debug("data_open(%s): fid=%lld", filename, (long long)fid);
		return (data_id)fid;
	}
	return DATA_FOLLOWER_FID;
}

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

int data_grp_create(data_id fid, const char *grp_name)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	if (WD.rank == 0) {
		/* Idempotent: if a real group already exists at this
		 * path, the caller (e.g. re-running ph2run on the
		 * same simul.h5) gets success without a fresh create.
		 * Stale dangling soft links (left over from old fixture
		 * builds) are unlinked so the create can proceed. */
		const htri_t lexists = H5Lexists(
			(hid_t)fid, grp_name, H5P_DEFAULT);
		if (lexists < 0) {
			log_error("data_grp_create(%s): H5Lexists failed",
				grp_name);
			rt = -1;
			goto ex_done;
		}
		if (lexists > 0) {
			const htri_t oexists = H5Oexists_by_name(
				(hid_t)fid, grp_name, H5P_DEFAULT);
			if (oexists > 0) {
				log_debug("data_grp_create(%s): already exists",
					grp_name);
				goto ex_done;
			}
			if (oexists < 0) {
				log_error("data_grp_create(%s):"
					  " H5Oexists_by_name failed",
					grp_name);
				rt = -1;
				goto ex_done;
			}
			/* Dangling link -- unlink and proceed. */
			if (H5Ldelete((hid_t)fid, grp_name, H5P_DEFAULT) < 0) {
				log_error("data_grp_create(%s): H5Ldelete of"
					  " dangling link failed", grp_name);
				rt = -1;
				goto ex_done;
			}
			log_debug("data_grp_create(%s): removed dangling"
				  " link", grp_name);
		}

		rt = -1;
		const hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
		if (lcpl == H5I_INVALID_HID) {
			log_error("data_grp_create(%s): H5Pcreate failed",
				grp_name);
			goto ex_lcpl;
		}
		if (H5Pset_create_intermediate_group(lcpl, 1) < 0) {
			log_error("data_grp_create(%s):"
				  " H5Pset_create_intermediate_group failed",
				grp_name);
			goto ex_prop;
		}
		if (H5Pset_char_encoding(lcpl, H5T_CSET_UTF8) < 0) {
			log_error("data_grp_create(%s):"
				  " H5Pset_char_encoding failed",
				grp_name);
			goto ex_prop;
		}
		const hid_t group = H5Gcreate((hid_t)fid, grp_name, lcpl,
			H5P_DEFAULT, H5P_DEFAULT);
		if (group == H5I_INVALID_HID) {
			log_error("data_grp_create(%s): H5Gcreate failed",
				grp_name);
			goto ex_group;
		}

		rt = 0;
		H5Gclose(group);
	ex_group:
	ex_prop:
		H5Pclose(lcpl);
	ex_lcpl:
	ex_done:;
	}
	bcast_int(&rt);
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
		bcast_int(&rt);                                                \
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
		bcast_int(&rt);                                                \
		return rt;                                                     \
	}

DEFINE_DATA_ATTR_WRITE(ul, unsigned long, H5T_NATIVE_ULONG);
DEFINE_DATA_ATTR_WRITE(dbl, double, H5T_NATIVE_DOUBLE);

/* -- multidet --------------------------------------------------------- */

#define MULTIDET_PATH DATA_STPREP "/" DATA_STPREP_MULTIDET

void data_multidet_free(struct data_multidet *m)
{
	if (!m)
		return;
	free((void *)m->cfs);
	free((void *)m->dets);
	m->cfs = NULL;
	m->dets = NULL;
}

int data_multidet_load(data_id fid, struct data_multidet *m)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	memset(m, 0, sizeof *m);

	int rt = 0;
	uint32_t v_nqb = 0;
	size_t v_ndets = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, MULTIDET_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_multidet_load: H5Gopen(%s) failed",
				MULTIDET_PATH);
			goto ex_dims;
		}
		hsize_t dims[2];
		if (get_dset_dims2(grp_id, DATA_STPREP_MULTIDET_DETS, dims)
			< 0) {
			H5Gclose(grp_id);
			goto ex_dims;
		}
		H5Gclose(grp_id);
		v_ndets = dims[0];
		v_nqb = (uint32_t)dims[1];
		rt = 0;
	ex_dims:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_u32(&v_nqb);
	bcast_sz(&v_ndets);

	double *cfs = malloc(sizeof *cfs * 2 * v_ndets);
	unsigned char *dets = malloc(sizeof *dets * v_ndets * v_nqb);
	if (!cfs || !dets) {
		log_error("data_multidet_load: alloc failed"
			  " (ndets=%zu, nqb=%u)", v_ndets, v_nqb);
		free(cfs);
		free(dets);
		return -1;
	}

	if (WD.rank == 0) {
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, MULTIDET_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_multidet_load: H5Gopen(%s) failed",
				MULTIDET_PATH);
			rt = -1;
		} else {
			if (read_dset(grp_id, DATA_STPREP_MULTIDET_COEFFS,
				    H5T_NATIVE_DOUBLE, cfs) < 0
				|| read_dset(grp_id, DATA_STPREP_MULTIDET_DETS,
					   H5T_NATIVE_UCHAR, dets) < 0)
				rt = -1;
			H5Gclose(grp_id);
		}
	}
	bcast_int(&rt);
	if (rt < 0) {
		free(cfs);
		free(dets);
		return -1;
	}
	bcast_doubles(cfs, 2 * v_ndets);
	bcast_uchars(dets, v_ndets * v_nqb);

	/* Bit-validate the determinant occupations on every rank
	 * (cheap, ndets*nqb bytes) so a malformed multidet group
	 * fails fast instead of corrupting downstream indices. */
	for (size_t i = 0; i < v_ndets; i++) {
		for (size_t j = 0; j < v_nqb; j++) {
			const unsigned char bit = dets[i * v_nqb + j];
			if (bit > 1) {
				log_error("data_multidet_load: dets[%zu][%zu]"
					  " = %u is not 0/1; multidet group"
					  " is malformed", i, j,
					(unsigned)bit);
				free(cfs);
				free(dets);
				return -1;
			}
		}
	}

	m->nqb = v_nqb;
	m->ndets = v_ndets;
	m->cfs = cfs;
	m->dets = dets;
	return 0;
}

/* -- pauli_hamil ------------------------------------------------------ */

void data_hamil_free(struct data_hamil *h)
{
	if (!h)
		return;
	free((void *)h->cfs);
	free((void *)h->paulis);
	h->cfs = NULL;
	h->paulis = NULL;
}

int data_hamil_load(data_id fid, struct data_hamil *h)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	memset(h, 0, sizeof *h);

	int rt = 0;
	uint32_t v_nqb = 0;
	size_t v_nterms = 0;
	double v_norm = 0.0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, DATA_HAMIL, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_hamil_load: H5Gopen(%s) failed",
				DATA_HAMIL);
			goto ex_dims;
		}
		hsize_t dims[2];
		if (get_dset_dims2(grp_id, DATA_HAMIL_PAULIS, dims) < 0
			|| read_attr_raw(grp_id, DATA_HAMIL_NORM,
				   H5T_NATIVE_DOUBLE, &v_norm) < 0) {
			H5Gclose(grp_id);
			goto ex_dims;
		}
		H5Gclose(grp_id);
		v_nterms = dims[0];
		v_nqb = (uint32_t)dims[1];
		rt = 0;
	ex_dims:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_u32(&v_nqb);
	bcast_sz(&v_nterms);
	MPI_Bcast(&v_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double *cfs = malloc(sizeof *cfs * v_nterms);
	unsigned char *paulis = malloc(sizeof *paulis * v_nterms * v_nqb);
	if (!cfs || !paulis) {
		log_error("data_hamil_load: alloc failed"
			  " (nterms=%zu, nqb=%u)", v_nterms, v_nqb);
		free(cfs);
		free(paulis);
		return -1;
	}

	if (WD.rank == 0) {
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, DATA_HAMIL, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_hamil_load: H5Gopen(%s) failed",
				DATA_HAMIL);
			rt = -1;
		} else {
			if (read_dset(grp_id, DATA_HAMIL_COEFFS,
				    H5T_NATIVE_DOUBLE, cfs) < 0
				|| read_dset(grp_id, DATA_HAMIL_PAULIS,
					   H5T_NATIVE_UCHAR, paulis) < 0)
				rt = -1;
			H5Gclose(grp_id);
		}
	}
	bcast_int(&rt);
	if (rt < 0) {
		free(cfs);
		free(paulis);
		return -1;
	}
	bcast_doubles(cfs, v_nterms);
	bcast_uchars(paulis, v_nterms * v_nqb);

	h->nqb = v_nqb;
	h->nterms = v_nterms;
	h->norm = v_norm;
	h->cfs = cfs;
	h->paulis = paulis;
	return 0;
}

/* -- per-step write API (/circ_{trott,trott2,qdrift,cmpsit}/values) --- */

#define CIRC_VALUES_DSET "values"

int data_circ_writer_init(data_id fid, const char *grp_name, size_t n_steps,
	struct data_circ_writer *w)
{
	memset(w, 0, sizeof *w);
	if (fid == 0)
		return 0;
	if (world_info(&WD) != WORLD_READY)
		return -1;

	/* Group create is collective (rank-0 does H5; all bcast). */
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

		/* Idempotence: reuse an existing values dataset. */
		const htri_t exists = H5Lexists(
			grp_id, CIRC_VALUES_DSET, H5P_DEFAULT);
		if (exists < 0) {
			log_error("data_circ_writer_init(%s): H5Lexists(values)"
				  " failed", grp_name);
			goto ex_grp;
		}
		if (exists > 0) {
			dset_out = H5Dopen2(
				grp_id, CIRC_VALUES_DSET, H5P_DEFAULT);
			if (dset_out == H5I_INVALID_HID) {
				log_error("data_circ_writer_init(%s):"
					  " H5Dopen2(values) failed", grp_name);
				goto ex_grp;
			}
			log_debug("data_circ_writer_init(%s): reusing existing"
				  " values dataset", grp_name);
			rt = 0;
			goto ex_grp;
		}

		const hid_t dspace = H5Screate_simple(
			2, (hsize_t[]){ n_steps, 2 }, NULL);
		if (dspace == H5I_INVALID_HID) {
			log_error("data_circ_writer_init(%s): H5Screate_simple"
				  " failed (n_steps=%zu)",
				grp_name, n_steps);
			goto ex_grp;
		}

		const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
		if (dcpl == H5I_INVALID_HID) {
			log_error("data_circ_writer_init(%s): H5Pcreate(dcpl)"
				  " failed", grp_name);
			goto ex_dspace;
		}
		const double nan_val = nan("");
		if (H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &nan_val) < 0) {
			log_error("data_circ_writer_init(%s):"
				  " H5Pset_fill_value failed", grp_name);
			goto ex_dcpl;
		}
		if (H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY) < 0) {
			log_error("data_circ_writer_init(%s):"
				  " H5Pset_alloc_time failed", grp_name);
			goto ex_dcpl;
		}
		if (H5Pset_fill_time(dcpl, H5D_FILL_TIME_ALLOC) < 0) {
			log_error("data_circ_writer_init(%s):"
				  " H5Pset_fill_time failed", grp_name);
			goto ex_dcpl;
		}

		dset_out = H5Dcreate2(grp_id, CIRC_VALUES_DSET, H5T_IEEE_F64LE,
			dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);
		if (dset_out == H5I_INVALID_HID) {
			log_error("data_circ_writer_init(%s):"
				  " H5Dcreate2(values) failed", grp_name);
			goto ex_dcpl;
		}

		/* Flush so the NaN-padded dataset is on disk before
		 * any per-step write.  Cheap for an empty dataset. */
		H5Fflush((hid_t)fid, H5F_SCOPE_GLOBAL);

		log_debug("data_circ_writer_init(%s): values shape (%zu, 2),"
			  " NaN-padded", grp_name, n_steps);
		rt = 0;
	ex_dcpl:
		H5Pclose(dcpl);
	ex_dspace:
		H5Sclose(dspace);
	ex_grp:
		H5Gclose(grp_id);
	ex_bcast:;
	}
	bcast_int(&rt);
	if (rt < 0) {
		if (WD.rank == 0 && dset_out != H5I_INVALID_HID)
			H5Dclose(dset_out);
		return -1;
	}
	w->fid = fid;
	w->n_steps = n_steps;
	w->dset = (WD.rank == 0) ? (int64_t)dset_out : 0;
	return 0;
}

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
	const hid_t mspace = H5Screate_simple(2, count, NULL);
	if (mspace == H5I_INVALID_HID) {
		log_error("data_circ_write_step(%zu): H5Screate_simple failed",
			step_idx);
		goto ex_fspace;
	}
	const double row[2] = { creal(z), cimag(z) };
	if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, row)
		< 0) {
		log_error("data_circ_write_step(%zu): H5Dwrite failed",
			step_idx);
		goto ex_mspace;
	}
	/* Atomic-on-disk: rows 0..step_idx are persisted before
	 * the simulation moves on. */
	H5Fflush((hid_t)w->fid, H5F_SCOPE_GLOBAL);

	rt = 0;
ex_mspace:
	H5Sclose(mspace);
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
	if (WD.rank == 0 && w->dset > 0)
		H5Dclose((hid_t)w->dset);
	memset(w, 0, sizeof *w);
}

/* -- state-prep dispatch ---------------------------------------------- */

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
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_int(&kind);
	*out = (enum stprep_kind)kind;
	return 0;
}

/* -- coeff_matrix ----------------------------------------------------- */

#define COEFFMAT_PATH DATA_STPREP "/" DATA_STPREP_COEFFMAT

static int read_u32_attr(hid_t grp_id, const char *name, uint32_t *out)
{
	return read_attr_raw(grp_id, name, H5T_NATIVE_UINT32, out);
}

static int read_u8_attr(hid_t grp_id, const char *name, int *out)
{
	uint8_t v = 0;
	if (read_attr_raw(grp_id, name, H5T_NATIVE_UINT8, &v) < 0)
		return -1;
	*out = v ? 1 : 0;
	return 0;
}

static int read_double_attr(hid_t grp_id, const char *name, double *out)
{
	return read_attr_raw(grp_id, name, H5T_NATIVE_DOUBLE, out);
}

static int read_coeff_attrs(hid_t grp_id, struct data_coeff_matrix *cm)
{
	if (read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NQB, &cm->nqb) < 0
		|| read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NS,
			   &cm->n_sites) < 0
		|| read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NA,
			   &cm->n_alpha) < 0
		|| read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NB,
			   &cm->n_beta) < 0
		|| read_u8_attr(grp_id, DATA_STPREP_COEFFMAT_CS,
			   &cm->closed_shell) < 0
		|| read_u8_attr(grp_id, DATA_STPREP_COEFFMAT_TAP,
			   &cm->tapered) < 0)
		return -1;
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

int data_coeff_matrix_load(const data_id fid, struct data_coeff_matrix *cm)
{
	memset(cm, 0, sizeof *cm);
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	uint32_t v_ncomp = 0;
	int v_has_csf = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, COEFFMAT_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_load: H5Gopen(%s) failed",
				COEFFMAT_PATH);
			goto ex_attrs;
		}
		if (read_coeff_attrs(grp_id, cm) < 0)
			goto ex_attrs_grp;
		const htri_t has = H5Lexists(
			grp_id, DATA_STPREP_COEFFMAT_CSF, H5P_DEFAULT);
		if (has < 0) {
			log_error("data_coeff_matrix_load: H5Lexists(csf)"
				  " failed");
			goto ex_attrs_grp;
		}
		if (has > 0) {
			v_has_csf = 1;
			const hid_t cg = H5Gopen(grp_id,
				DATA_STPREP_COEFFMAT_CSF, H5P_DEFAULT);
			if (cg == H5I_INVALID_HID) {
				log_error("data_coeff_matrix_load:"
					  " H5Gopen(csf) failed");
				goto ex_attrs_grp;
			}
			const int rt_n = read_u32_attr(
				cg, DATA_STPREP_COEFFMAT_CSF_NCOMP, &v_ncomp);
			H5Gclose(cg);
			if (rt_n < 0)
				goto ex_attrs_grp;
			/* An explicit csf/ subgroup with n_components==0
			 * is non-physical: a single-block file must omit
			 * the subgroup, not declare an empty
			 * superposition.  Reject early so a misbuilt
			 * simul.h5 never reaches the simulator. */
			if (v_ncomp == 0) {
				log_error("simul.h5: /state_prep/coeff_matrix"
					  "/csf present with n_components=0;"
					  " remove the csf subgroup or list"
					  " at least one component");
				goto ex_attrs_grp;
			}
		}
		rt = 0;
	ex_attrs_grp:
		H5Gclose(grp_id);
	ex_attrs:;
	}
	bcast_int(&rt);
	if (rt < 0) {
		memset(cm, 0, sizeof *cm);
		return -1;
	}
	bcast_u32(&cm->nqb);
	bcast_u32(&cm->n_sites);
	bcast_u32(&cm->n_alpha);
	bcast_u32(&cm->n_beta);
	bcast_int(&cm->closed_shell);
	bcast_int(&cm->tapered);
	bcast_int(&v_has_csf);
	bcast_u32(&v_ncomp);

	const size_t sz_a = (size_t)cm->n_sites * cm->n_alpha;
	const size_t sz_b = (size_t)cm->n_sites * cm->n_beta;

	/* Allocate through local mutable pointers; the struct
	 * fields are const-qualified and assigned only after
	 * loads + bcasts complete.  On failure we free locally
	 * and zero the struct. */
	double *Ca = NULL, *Cb = NULL;
	struct data_coeff_block *blocks = NULL;
	if (v_has_csf) {
		blocks = calloc(v_ncomp, sizeof *blocks);
		if (!blocks) {
			log_error("data_coeff_matrix_load: alloc blocks"
				  " (n=%u) failed", v_ncomp);
			goto err_alloc;
		}
		for (size_t k = 0; k < v_ncomp; k++) {
			double *bca = malloc(sizeof(double) * (sz_a ? sz_a : 1));
			if (!bca) {
				log_error("data_coeff_matrix_load: alloc"
					  " C_alpha[%zu] failed", k);
				goto err_alloc;
			}
			blocks[k].C_alpha = bca;
			if (!cm->closed_shell) {
				double *bcb = malloc(
					sizeof(double) * (sz_b ? sz_b : 1));
				if (!bcb) {
					log_error("data_coeff_matrix_load:"
						  " alloc C_beta[%zu] failed",
						k);
					goto err_alloc;
				}
				blocks[k].C_beta = bcb;
			}
		}
	} else {
		Ca = malloc(sizeof(double) * (sz_a ? sz_a : 1));
		if (!Ca) {
			log_error("data_coeff_matrix_load: alloc C_alpha"
				  " failed");
			goto err_alloc;
		}
		if (!cm->closed_shell) {
			Cb = malloc(sizeof(double) * (sz_b ? sz_b : 1));
			if (!Cb) {
				log_error("data_coeff_matrix_load: alloc"
					  " C_beta failed");
				goto err_alloc;
			}
		}
	}

	rt = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, COEFFMAT_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_load: H5Gopen(%s) failed"
				  " (dsets)", COEFFMAT_PATH);
			goto ex_dsets;
		}
		if (v_has_csf) {
			int ok = 1;
			for (size_t k = 0; k < v_ncomp && ok; k++) {
				char sub[32];
				snprintf(sub, sizeof sub, "%s/%zu",
					DATA_STPREP_COEFFMAT_CSF, k);
				const hid_t cg = H5Gopen(
					grp_id, sub, H5P_DEFAULT);
				if (cg == H5I_INVALID_HID) {
					log_error("data_coeff_matrix_load:"
						  " H5Gopen(%s/%s) failed",
						COEFFMAT_PATH, sub);
					ok = 0;
					break;
				}
				if (read_double_attr(cg,
					    DATA_STPREP_COEFFMAT_CSF_CF,
					    &blocks[k].cf) < 0
					|| read_C_dset(cg,
						   DATA_STPREP_COEFFMAT_CA,
						   cm->n_sites, cm->n_alpha,
						   (double *)blocks[k].C_alpha)
						< 0
					|| (!cm->closed_shell
						&& read_C_dset(cg,
							   DATA_STPREP_COEFFMAT_CB,
							   cm->n_sites,
							   cm->n_beta,
							   (double *)blocks[k]
								   .C_beta)
							< 0))
					ok = 0;
				H5Gclose(cg);
			}
			if (ok)
				rt = 0;
		} else {
			if (read_C_dset(grp_id, DATA_STPREP_COEFFMAT_CA,
				    cm->n_sites, cm->n_alpha, Ca) >= 0
				&& (cm->closed_shell
					|| read_C_dset(grp_id,
						   DATA_STPREP_COEFFMAT_CB,
						   cm->n_sites, cm->n_beta,
						   Cb) >= 0))
				rt = 0;
		}
		H5Gclose(grp_id);
	ex_dsets:;
	}
	bcast_int(&rt);
	if (rt < 0)
		goto err_alloc;

	if (v_has_csf) {
		for (size_t k = 0; k < v_ncomp; k++) {
			MPI_Bcast(&blocks[k].cf, 1, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
			bcast_doubles((double *)blocks[k].C_alpha, sz_a);
			if (!cm->closed_shell)
				bcast_doubles((double *)blocks[k].C_beta, sz_b);
		}
		cm->n_components = v_ncomp;
		cm->blocks = blocks;
	} else {
		bcast_doubles(Ca, sz_a);
		if (!cm->closed_shell)
			bcast_doubles(Cb, sz_b);
		cm->C_alpha = Ca;
		cm->C_beta = Cb;
	}
	return 0;

err_alloc:
	free(Ca);
	free(Cb);
	if (blocks) {
		for (size_t k = 0; k < v_ncomp; k++) {
			free((void *)blocks[k].C_alpha);
			free((void *)blocks[k].C_beta);
		}
		free(blocks);
	}
	memset(cm, 0, sizeof *cm);
	return -1;
}

void data_coeff_matrix_free(struct data_coeff_matrix *cm)
{
	if (!cm)
		return;
	if (cm->blocks) {
		for (size_t k = 0; k < cm->n_components; k++) {
			free((void *)cm->blocks[k].C_alpha);
			free((void *)cm->blocks[k].C_beta);
		}
		free(cm->blocks);
	}
	free((void *)cm->C_alpha);
	free((void *)cm->C_beta);
	memset(cm, 0, sizeof *cm);
}
