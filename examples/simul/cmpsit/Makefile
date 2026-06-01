# -----------------------------------------------------------------------------
# Composite (partially randomised) 2nd-order Trotter: prepare
# simul.h5, run ph2run, estimate the ground-state energy.
#
# The energy estimate is accurate when ANGLE_RAND ~ ANGLE_DET (the
# randomised tail then shares the deterministic time scale) and
# STEPS * ANGLE_DET is small enough to keep the phase in (-pi, pi).
#
# Simulation parameters (override on the command line):
LENGTH		?= 200
DEPTH		?= 4
ANGLE_DET	?= 0.0625
ANGLE_RAND	?= 0.0625
STEPS		?= 8
SAMPLES		?= 256
SEED		?= 9393992002
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
ENERGY_MC	= $(PYTHON) $(SCRIPTS)/energy_mc.py --group circ_cmpsit

PH2RUN_FLAGS	= --simul=$(SIMUL_H5) cmpsit		\
			--length=$(LENGTH) --depth=$(DEPTH)	\
			--angle-det=$(ANGLE_DET)		\
			--angle-rand=$(ANGLE_RAND)		\
			--steps=$(STEPS) --samples=$(SAMPLES)	\
			--seed=$(SEED)

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
