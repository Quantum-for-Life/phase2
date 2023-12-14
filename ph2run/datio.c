#include "data.h"

#include "log.h"

int datio_read_file(struct data *dat, const char *filename)
{
	int rc = 0;

	log_debug("Read simulation input file: %s", filename);
	data_id fid = data_file_open(filename);
	if (fid == DATA_INVALID_FID) {
		log_error("Open data file");
		return -1;
	}
	if (data_parse(dat, fid) < 0) {
		log_error("Failure: read data file");
		rc = -1;
	}
	data_file_close(fid);

	return rc;
}

void datio_print_info(const struct data *dat)
{
	log_debug("State preparation:");
	log_debug("multidet, num_qubits=%zu, num_terms=%zu",
		  dat->state_prep.multidet.num_qubits,
		  dat->state_prep.multidet.num_terms);
	log_debug("Hamiltonian: num_qubits=%zu, num_terms=%zu, "
		  "norm=%f",
		  dat->pauli_hamil.num_qubits, dat->pauli_hamil.num_terms,
		  dat->pauli_hamil.norm);
	log_debug("Time series: num_steps=%zu", dat->time_series.num_steps);
}

int datio_save_file(const char *filename, const struct data *dat)
{
	int rc = 0;

	log_debug("Saving data");
	data_id fid = data_file_open(filename);
	if (fid == DATA_INVALID_FID) {
		log_error("Open data file");
		return -1;
	}
	if (data_time_series_write(fid, &dat->time_series) < 0) {
		log_error("writing data");
		rc = -1;
	}
	data_file_close(fid);

	return rc;
}
