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

#include "hdf5.h"
#include "mpi.h"

#include "log.h"
#include "phase2/data.h"
#include "phase2/world.h"
#include <complex.h>

static struct world_info WD;

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

int data_multidet_getnums(data_id fid, uint32_t *nqb, size_t *ndets)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	uint32_t v_nqb = 0;
	size_t v_ndets = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, MULTIDET_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_multidet_getnums: H5Gopen(%s) failed",
				MULTIDET_PATH);
			goto ex_open;
		}
		hsize_t dims[2];
		if (get_dset_dims2(grp_id, DATA_STPREP_MULTIDET_DETS, dims) < 0)
			goto ex_grp;
		v_ndets = dims[0];
		v_nqb = (uint32_t)dims[1];
		rt = 0;
	ex_grp:
		H5Gclose(grp_id);
	ex_open:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_u32(&v_nqb);
	bcast_sz(&v_ndets);
	*nqb = v_nqb;
	*ndets = v_ndets;
	return 0;
}

static int multidet_read_data(
	hid_t fid, double *cfs, unsigned char *dets)
{
	const hid_t grp_id = H5Gopen(fid, MULTIDET_PATH, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		log_error("multidet_read_data: H5Gopen(%s) failed",
			MULTIDET_PATH);
		return -1;
	}
	int rt = -1;
	if (read_dset(grp_id, DATA_STPREP_MULTIDET_COEFFS, H5T_NATIVE_DOUBLE,
		    cfs) < 0)
		goto ex;
	if (read_dset(grp_id, DATA_STPREP_MULTIDET_DETS, H5T_NATIVE_UCHAR,
		    dets) < 0)
		goto ex;
	rt = 0;
ex:
	H5Gclose(grp_id);
	return rt;
}

int data_multidet_foreach(data_id fid,
	int (*op)(_Complex double cf, uint64_t idx, void *), void *op_data)
{
	int rt = -1, rc = 0;
	uint32_t nqb;
	size_t ndets;

	if (data_multidet_getnums(fid, &nqb, &ndets) < 0)
		return -1;

	double *cfs = malloc(sizeof *cfs * 2 * ndets);
	if (!cfs) {
		log_error("data_multidet_foreach: alloc cfs (%zu bytes)",
			sizeof *cfs * 2 * ndets);
		goto ex_alloc_coeffs;
	}
	unsigned char *dets = malloc(sizeof *dets * ndets * nqb);
	if (!dets) {
		log_error("data_multidet_foreach: alloc dets (%zu bytes)",
			sizeof *dets * ndets * nqb);
		goto ex_alloc_dets;
	}

	int status = 0;
	if (WD.rank == 0) {
		if (multidet_read_data((hid_t)fid, cfs, dets) < 0)
			status = -1;
	}
	bcast_int(&status);
	if (status < 0)
		goto ex_data_read;
	bcast_doubles(cfs, 2 * ndets);
	bcast_uchars(dets, ndets * nqb);

	for (size_t i = 0; i < ndets; i++) {
		uint64_t idx = 0;
		for (size_t j = 0; j < nqb; j++) {
			const unsigned char bit = dets[i * nqb + j];
			if (bit > 1) {
				log_error("data_multidet_foreach: dets[%zu][%zu]"
					  " = %u is not 0/1; multidet group"
					  " is malformed",
					i, j, (unsigned)bit);
				rt = -1;
				goto ex_data_read;
			}
			idx += (uint64_t)bit << j;
		}
		_Complex double cf = CMPLX(cfs[2 * i], cfs[2 * i + 1]);
		rc = op(cf, idx, op_data);
		/* Caller short-circuits iteration with a non-zero
		 * return; that isn't an error. */
		if (rc != 0)
			break;
	}
	rt = rc;

ex_data_read:
	free(dets);
ex_alloc_dets:
	free(cfs);
ex_alloc_coeffs:
	return rt;
}

/* -- pauli_hamil ------------------------------------------------------ */

int data_hamil_getnums(data_id fid, uint32_t *nqb, size_t *nterms)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	uint32_t v_nqb = 0;
	size_t v_nterms = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, DATA_HAMIL, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_hamil_getnums: H5Gopen(%s) failed",
				DATA_HAMIL);
			goto ex_grp;
		}
		hsize_t dims[2];
		if (get_dset_dims2(grp_id, DATA_HAMIL_PAULIS, dims) < 0)
			goto ex_dims;
		v_nterms = dims[0];
		v_nqb = (uint32_t)dims[1];
		rt = 0;
	ex_dims:
		H5Gclose(grp_id);
	ex_grp:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_u32(&v_nqb);
	bcast_sz(&v_nterms);
	*nqb = v_nqb;
	*nterms = v_nterms;
	return 0;
}

int data_hamil_getnorm(data_id fid, double *norm)
{
	return data_attr_read_dbl(fid, DATA_HAMIL, DATA_HAMIL_NORM, norm);
}

static int hamil_read_data(hid_t fid, double *cfs, unsigned char *paulis)
{
	const hid_t grp_id = H5Gopen(fid, DATA_HAMIL, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		log_error("hamil_read_data: H5Gopen(%s) failed", DATA_HAMIL);
		return -1;
	}
	int rt = -1;
	if (read_dset(grp_id, DATA_HAMIL_COEFFS, H5T_NATIVE_DOUBLE, cfs) < 0)
		goto ex;
	if (read_dset(grp_id, DATA_HAMIL_PAULIS, H5T_NATIVE_UCHAR, paulis) < 0)
		goto ex;
	rt = 0;
ex:
	H5Gclose(grp_id);
	return rt;
}

int data_hamil_foreach(const data_id fid,
	int (*op)(double, unsigned char *, void *), void *op_data)
{
	int rt = -1, rc = 0;
	uint32_t nqb;
	size_t nterms;

	if (data_hamil_getnums(fid, &nqb, &nterms) < 0)
		return -1;

	double *cfs = malloc(sizeof *cfs * nterms);
	if (!cfs) {
		log_error("data_hamil_foreach: alloc cfs (%zu bytes)",
			sizeof *cfs * nterms);
		goto ex_coeffs_alloc;
	}
	unsigned char *paulis = malloc(sizeof *paulis * nqb * nterms);
	if (!paulis) {
		log_error("data_hamil_foreach: alloc paulis (%zu bytes)",
			sizeof *paulis * nqb * nterms);
		goto ex_paulis_alloc;
	}

	int status = 0;
	if (WD.rank == 0) {
		if (hamil_read_data((hid_t)fid, cfs, paulis) < 0)
			status = -1;
	}
	bcast_int(&status);
	if (status < 0)
		goto ex_hamil_read;
	bcast_doubles(cfs, nterms);
	bcast_uchars(paulis, nqb * nterms);

	unsigned char *paustr = malloc(sizeof *paustr * nqb);
	if (!paustr) {
		log_error("data_hamil_foreach: alloc paustr (%zu bytes)",
			sizeof *paustr * nqb);
		goto ex_paustr_alloc;
	}

	for (size_t i = 0; i < nterms; i++) {
		for (size_t j = 0; j < nqb; j++)
			paustr[j] = paulis[i * nqb + j];
		rc = op(cfs[i], paustr, op_data);
		if (rc != 0)
			break;
	}

	rt = rc;
	free(paustr);
ex_paustr_alloc:
ex_hamil_read:
	free(paulis);
ex_paulis_alloc:
	free(cfs);
ex_coeffs_alloc:
	return rt;
}

/* -- per-step write API (/circ_{trott,trott2,qdrift,cmpsit}/values) --- */

#define CIRC_VALUES_DSET "values"

int data_circ_init(data_id fid, const char *grp_name, size_t n_steps)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	/* Group create is collective (rank-0 does H5; all bcast). */
	if (data_grp_create(fid, grp_name) < 0)
		return -1;

	int rt = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, grp_name, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_circ_init(%s): H5Gopen failed",
				grp_name);
			goto ex_bcast;
		}

		/* Idempotence: skip if the dataset already exists. */
		const htri_t exists = H5Lexists(
			grp_id, CIRC_VALUES_DSET, H5P_DEFAULT);
		if (exists < 0) {
			log_error("data_circ_init(%s): H5Lexists(values)"
				  " failed", grp_name);
			goto ex_grp;
		}
		if (exists > 0) {
			log_debug("data_circ_init(%s): values dataset"
				  " already exists", grp_name);
			rt = 0;
			goto ex_grp;
		}

		const hid_t dspace = H5Screate_simple(
			2, (hsize_t[]){ n_steps, 2 }, NULL);
		if (dspace == H5I_INVALID_HID) {
			log_error("data_circ_init(%s): H5Screate_simple"
				  " failed (n_steps=%zu)",
				grp_name, n_steps);
			goto ex_grp;
		}

		const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
		if (dcpl == H5I_INVALID_HID) {
			log_error("data_circ_init(%s): H5Pcreate(dcpl)"
				  " failed", grp_name);
			goto ex_dspace;
		}
		const double nan_val = nan("");
		if (H5Pset_fill_value(dcpl, H5T_NATIVE_DOUBLE, &nan_val) < 0) {
			log_error("data_circ_init(%s): H5Pset_fill_value"
				  " failed", grp_name);
			goto ex_dcpl;
		}
		if (H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_EARLY) < 0) {
			log_error("data_circ_init(%s): H5Pset_alloc_time"
				  " failed", grp_name);
			goto ex_dcpl;
		}
		if (H5Pset_fill_time(dcpl, H5D_FILL_TIME_ALLOC) < 0) {
			log_error("data_circ_init(%s): H5Pset_fill_time"
				  " failed", grp_name);
			goto ex_dcpl;
		}

		const hid_t dset = H5Dcreate2(grp_id, CIRC_VALUES_DSET,
			H5T_IEEE_F64LE, dspace, H5P_DEFAULT, dcpl,
			H5P_DEFAULT);
		if (dset == H5I_INVALID_HID) {
			log_error("data_circ_init(%s): H5Dcreate2(values)"
				  " failed", grp_name);
			goto ex_dcpl;
		}
		H5Dclose(dset);

		/* Flush so the NaN-padded dataset is on disk before
		 * any per-step write.  Cheap for an empty dataset. */
		H5Fflush((hid_t)fid, H5F_SCOPE_GLOBAL);

		log_debug("data_circ_init(%s): values shape (%zu, 2),"
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
	return rt;
}

int data_circ_write_step(data_id fid, const char *grp_name, size_t step_idx,
	_Complex double z)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;
	if (WD.rank != 0)
		return 0;

	int rt = -1;
	const hid_t grp_id = H5Gopen((hid_t)fid, grp_name, H5P_DEFAULT);
	if (grp_id == H5I_INVALID_HID) {
		log_error("data_circ_write_step(%s, %zu): H5Gopen failed",
			grp_name, step_idx);
		return -1;
	}
	const hid_t dset = H5Dopen2(grp_id, CIRC_VALUES_DSET, H5P_DEFAULT);
	if (dset == H5I_INVALID_HID) {
		log_error("data_circ_write_step(%s, %zu): H5Dopen2(values)"
			  " failed", grp_name, step_idx);
		goto ex_grp;
	}
	const hid_t fspace = H5Dget_space(dset);
	if (fspace == H5I_INVALID_HID) {
		log_error("data_circ_write_step(%s, %zu): H5Dget_space"
			  " failed", grp_name, step_idx);
		goto ex_dset;
	}
	const hsize_t start[2] = { step_idx, 0 };
	const hsize_t count[2] = { 1, 2 };
	if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count,
		    NULL) < 0) {
		log_error("data_circ_write_step(%s, %zu): H5Sselect_hyperslab"
			  " failed", grp_name, step_idx);
		goto ex_fspace;
	}
	const hid_t mspace = H5Screate_simple(2, count, NULL);
	if (mspace == H5I_INVALID_HID) {
		log_error("data_circ_write_step(%s, %zu): H5Screate_simple"
			  " failed", grp_name, step_idx);
		goto ex_fspace;
	}
	const double row[2] = { creal(z), cimag(z) };
	if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, row)
		< 0) {
		log_error("data_circ_write_step(%s, %zu): H5Dwrite failed",
			grp_name, step_idx);
		goto ex_mspace;
	}
	/* Atomic-on-disk: the next H5Fflush call ensures rows
	 * 0..step_idx are persisted before the simulation moves
	 * on. */
	H5Fflush((hid_t)fid, H5F_SCOPE_GLOBAL);

	rt = 0;
ex_mspace:
	H5Sclose(mspace);
ex_fspace:
	H5Sclose(fspace);
ex_dset:
	H5Dclose(dset);
ex_grp:
	H5Gclose(grp_id);
	return rt;
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

int data_coeff_matrix_getnums(const data_id fid, uint32_t *nqb,
	uint32_t *n_sites, uint32_t *n_alpha, uint32_t *n_beta,
	int *closed_shell, int *tapered)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	uint32_t v_nqb = 0, v_ns = 0, v_na = 0, v_nb = 0;
	int v_cs = 0, v_tap = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, COEFFMAT_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_getnums:"
				  " H5Gopen(%s) failed", COEFFMAT_PATH);
			goto ex_open;
		}

		if (read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NQB, &v_nqb) < 0
			|| read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NS,
				   &v_ns) < 0
			|| read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NA,
				   &v_na) < 0
			|| read_u32_attr(grp_id, DATA_STPREP_COEFFMAT_NB,
				   &v_nb) < 0
			|| read_u8_attr(grp_id, DATA_STPREP_COEFFMAT_CS,
				   &v_cs) < 0
			|| read_u8_attr(grp_id, DATA_STPREP_COEFFMAT_TAP,
				   &v_tap) < 0)
			goto ex_grp;
		rt = 0;
	ex_grp:
		H5Gclose(grp_id);
	ex_open:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_u32(&v_nqb);
	bcast_u32(&v_ns);
	bcast_u32(&v_na);
	bcast_u32(&v_nb);
	bcast_int(&v_cs);
	bcast_int(&v_tap);
	*nqb = v_nqb;
	*n_sites = v_ns;
	*n_alpha = v_na;
	*n_beta = v_nb;
	*closed_shell = v_cs;
	*tapered = v_tap;
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

int data_coeff_matrix_read(
	const data_id fid, double *C_alpha, double *C_beta)
{
	uint32_t nqb, n_sites, n_alpha, n_beta;
	int closed_shell, tapered;
	if (data_coeff_matrix_getnums(fid, &nqb, &n_sites, &n_alpha, &n_beta,
		    &closed_shell, &tapered) < 0)
		return -1;

	if (closed_shell && C_beta != NULL) {
		log_error("data_coeff_matrix_read: closed_shell=1 but"
			  " C_beta != NULL");
		return -1;
	}
	if (!closed_shell && C_beta == NULL) {
		log_error("data_coeff_matrix_read: closed_shell=0 but"
			  " C_beta == NULL");
		return -1;
	}

	int rt = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, COEFFMAT_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_read: H5Gopen(%s) failed",
				COEFFMAT_PATH);
			goto ex_open;
		}
		if (read_C_dset(grp_id, DATA_STPREP_COEFFMAT_CA, n_sites,
			    n_alpha, C_alpha) < 0)
			goto ex_grp;
		if (!closed_shell
			&& read_C_dset(grp_id, DATA_STPREP_COEFFMAT_CB, n_sites,
				   n_beta, C_beta) < 0)
			goto ex_grp;
		rt = 0;
	ex_grp:
		H5Gclose(grp_id);
	ex_open:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_doubles(C_alpha, (size_t)n_sites * n_alpha);
	if (!closed_shell)
		bcast_doubles(C_beta, (size_t)n_sites * n_beta);
	return 0;
}

int data_coeff_matrix_csf_count(const data_id fid, size_t *n)
{
	if (world_info(&WD) != WORLD_READY)
		return -1;

	int rt = 0;
	uint32_t ncomp = 0;
	int present = 0;
	if (WD.rank == 0) {
		rt = -1;
		const hid_t grp_id = H5Gopen(
			(hid_t)fid, COEFFMAT_PATH, H5P_DEFAULT);
		if (grp_id == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_csf_count:"
				  " H5Gopen(%s) failed", COEFFMAT_PATH);
			goto ex_open;
		}

		const htri_t has_csf = H5Lexists(
			grp_id, DATA_STPREP_COEFFMAT_CSF, H5P_DEFAULT);
		if (has_csf < 0) {
			log_error("data_coeff_matrix_csf_count: H5Lexists"
				  " failed");
			goto ex_grp;
		}
		if (!has_csf) {
			present = 0;
			ncomp = 0;
			rt = 0;
			goto ex_grp;
		}
		present = 1;

		const hid_t cg = H5Gopen(grp_id, DATA_STPREP_COEFFMAT_CSF,
			H5P_DEFAULT);
		if (cg == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_csf_count: H5Gopen(%s)"
				  " failed", DATA_STPREP_COEFFMAT_CSF);
			goto ex_grp;
		}
		if (read_u32_attr(
			    cg, DATA_STPREP_COEFFMAT_CSF_NCOMP, &ncomp) < 0) {
			H5Gclose(cg);
			goto ex_grp;
		}
		H5Gclose(cg);

		/* An explicit csf/ subgroup with n_components==0 is
		 * non-physical: a single-block file must omit the
		 * subgroup, not declare an empty superposition.
		 * Reject early so a misbuilt simul.h5 never reaches
		 * the simulator. */
		if (ncomp == 0) {
			log_error("simul.h5: /state_prep/coeff_matrix/csf"
				  " present with n_components=0; remove the"
				  " csf subgroup or list at least one"
				  " component");
			goto ex_grp;
		}

		rt = 0;
	ex_grp:
		H5Gclose(grp_id);
	ex_open:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	bcast_int(&present);
	bcast_u32(&ncomp);
	*n = present ? (size_t)ncomp : 0;
	return 0;
}

int data_coeff_matrix_csf_read(const data_id fid, const size_t k,
	double *coefficient, double *C_alpha, double *C_beta)
{
	uint32_t nqb, n_sites, n_alpha, n_beta;
	int closed_shell, tapered;
	if (data_coeff_matrix_getnums(fid, &nqb, &n_sites, &n_alpha, &n_beta,
		    &closed_shell, &tapered) < 0)
		return -1;

	if (closed_shell && C_beta != NULL) {
		log_error("data_coeff_matrix_csf_read[%zu]: closed_shell=1"
			  " but C_beta != NULL", k);
		return -1;
	}
	if (!closed_shell && C_beta == NULL) {
		log_error("data_coeff_matrix_csf_read[%zu]: closed_shell=0"
			  " but C_beta == NULL", k);
		return -1;
	}

	int rt = 0;
	if (WD.rank == 0) {
		rt = -1;
		char path[64];
		snprintf(path, sizeof path, "%s/%s/%zu", COEFFMAT_PATH,
			DATA_STPREP_COEFFMAT_CSF, k);
		const hid_t cg_k = H5Gopen((hid_t)fid, path, H5P_DEFAULT);
		if (cg_k == H5I_INVALID_HID) {
			log_error("data_coeff_matrix_csf_read[%zu]:"
				  " H5Gopen(%s) failed", k, path);
			goto ex_open;
		}

		if (read_double_attr(cg_k, DATA_STPREP_COEFFMAT_CSF_CF,
			    coefficient) < 0)
			goto ex_kgrp;
		if (read_C_dset(cg_k, DATA_STPREP_COEFFMAT_CA, n_sites,
			    n_alpha, C_alpha) < 0)
			goto ex_kgrp;
		if (!closed_shell
			&& read_C_dset(cg_k, DATA_STPREP_COEFFMAT_CB, n_sites,
				   n_beta, C_beta) < 0)
			goto ex_kgrp;
		rt = 0;
	ex_kgrp:
		H5Gclose(cg_k);
	ex_open:;
	}
	bcast_int(&rt);
	if (rt < 0)
		return -1;
	MPI_Bcast(coefficient, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	bcast_doubles(C_alpha, (size_t)n_sites * n_alpha);
	if (!closed_shell)
		bcast_doubles(C_beta, (size_t)n_sites * n_beta);
	return 0;
}
