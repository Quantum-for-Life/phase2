# -----------------------------------------------------------------------------
# qDRIFT randomised product formula with robust phase estimation.
#
# RPE drives a sweep of J+1 runs at a fixed step size x = 2^-J and
# doubling depth (2^J, 2^(J+1), ..., 2^(2J)), where
# J = ceil(log2(DELTA / EPSILON)).  qdrift_rpe.py unwraps the
# accumulated phase across the sweep into a ground-state energy.
#
# Simulation parameters (override on the command line):
DELTA		?= 1.0
EPSILON		?= 0.125
SAMPLES		?= 64
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
QDRIFT_RPE	= $(PYTHON) $(SCRIPTS)/qdrift_rpe.py

# Sweep schedule derived from DELTA / EPSILON (GNU make + coreutils).
J		:= $(shell $(PYTHON) -c "import math;print(math.ceil(math.log2($(DELTA)/$(EPSILON))))")
STEP_SIZE	:= $(shell $(PYTHON) -c "print(2.0**-$(J))")
LEVELS		:= $(shell seq 0 $(J))

.DEFAULT_GOAL := all
.PHONY: all sweep analyze regen clean

all: $(SIMUL_PROC)

# One run per RPE level i: depth = 2^(i+J), constant step size.
# Each level lands in its own simul.h5-<i> file.
sweep:
	@for i in $(LEVELS); do						\
		depth=$$($(PYTHON) -c "print(int(2**($$i+$(J))))");	\
		echo "qdrift RPE level $$i: depth=$$depth"		\
			"step_size=$(STEP_SIZE)";			\
		cp $(HAMIL_H5) $(SIMUL_H5)-$$i;				\
		$(MPIRUN) -n $(MPI_RANKS) -x PHASE2_LOG=$(PHASE2_LOG)	\
			$(PH2RUN) --simul=$(SIMUL_H5)-$$i qdrift	\
			--step-size=$(STEP_SIZE) --depth=$$depth	\
			--samples=$(SAMPLES) --seed=$(SEED);		\
	done

# The base simul.h5 carries /pauli_hamil (normalization, offset);
# qdrift_rpe.py reads the per-level simul.h5-<i> circuit results.
$(SIMUL_PROC): sweep
	cp $(HAMIL_H5) $(SIMUL_H5)
	$(QDRIFT_RPE) --delta=$(DELTA) --epsilon=$(EPSILON)		\
		--num-samples=$(SAMPLES) $(SIMUL_H5) | tee $@

analyze: $(SIMUL_PROC)

regen:
	$(PARSE_FCIDUMP) --sort-terms -o $(HAMIL_H5) $(FCIDUMP)
	$(PARSE_INPUTST) -o $(HAMIL_H5) $(INPUTST)

clean:
	$(RM) $(SIMUL_H5) $(SIMUL_H5)-* $(SIMUL_PROC)
