# -----------------------------------------------------------------------------
# qDRIFT randomised product formula: prepare simul.h5, run ph2run,
# estimate the ground-state energy from the averaged sample overlap
# (energy_mc.py).  One run at fixed step size and depth; the analysis
# averages over SAMPLES independent draws.
#
# The estimate is accurate while the effective time
# T = depth * asin(step_size) is order 1 (the averaged overlap stays
# coherent).  Small step_size keeps the per-gate Trotter error low;
# raise SAMPLES to shrink the 1/sqrt(N) statistical error.
#
# Simulation parameters (override on the command line):
STEP_SIZE	?= 0.0625
DEPTH		?= 16
SAMPLES		?= 1024
SEED		?= 3924873492
MPI_RANKS	?= 2
# -----------------------------------------------------------------------------

ROOT_DIR	= ../../..
PH2RUN		?= $(ROOT_DIR)/build/ph2run/ph2run
SCRIPTS		= $(ROOT_DIR)/examples/scripts
DATA		= $(ROOT_DIR)/examples/data

HAMIL_H5	= $(DATA)/hamil.h5
FCIDUMP		= $(DATA)/FCIDUMP
INPUTST		= $(DATA)/INPUTST
SIMUL_H5	:= simul.h5
SIMUL_PROC	:= simul.h5.proc

PYTHON		?= python3
MPIRUN		?= mpirun
PHASE2_LOG	?= info

PARSE_FCIDUMP	= $(PYTHON) $(SCRIPTS)/parse_fcidump.py
PARSE_INPUTST	= $(PYTHON) $(SCRIPTS)/parse_inputst.py
ENERGY_MC	= $(PYTHON) $(SCRIPTS)/energy_mc.py --group circ_qdrift

PH2RUN_FLAGS	= --simul=$(SIMUL_H5) qdrift			\
			--step-size=$(STEP_SIZE) --depth=$(DEPTH)	\
			--samples=$(SAMPLES) --seed=$(SEED)

.DEFAULT_GOAL := all
.PHONY: all run analyze regen clean

all: $(SIMUL_PROC)

$(SIMUL_H5): $(HAMIL_H5)
	cp $(HAMIL_H5) $@
	$(MPIRUN) -n $(MPI_RANKS) -x PHASE2_LOG=$(PHASE2_LOG)	\
		$(PH2RUN) $(PH2RUN_FLAGS)

# Ground-state energy (CSV: E0, E0 + offset, in Hartree).
$(SIMUL_PROC): $(SIMUL_H5)
	$(ENERGY_MC) $(SIMUL_H5) | tee $@

run: $(SIMUL_H5)
analyze: $(SIMUL_PROC)

regen:
	$(PARSE_FCIDUMP) --sort-terms -o $(HAMIL_H5) $(FCIDUMP)
	$(PARSE_INPUTST) -o $(HAMIL_H5) $(INPUTST)

clean:
	$(RM) $(SIMUL_H5) $(SIMUL_PROC)
