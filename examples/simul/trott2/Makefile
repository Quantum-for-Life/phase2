# -----------------------------------------------------------------------------
# 2nd-order symmetric (Strang) Trotter: prepare simul.h5, run ph2run,
# estimate the ground-state energy by robust phase estimation.
#
# Simulation parameters (override on the command line):
DELTA		?= 0.1
TROTT_STEPS	?= 1024
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
TROTT_RPE	= $(PYTHON) $(SCRIPTS)/trott_rpe.py --group circ_trott2

PH2RUN_FLAGS	= --simul=$(SIMUL_H5) trott2		\
			--delta=$(DELTA) --steps=$(TROTT_STEPS)

.DEFAULT_GOAL := all
.PHONY: all run analyze regen clean

all: $(SIMUL_PROC)

$(SIMUL_H5): $(HAMIL_H5)
	cp $(HAMIL_H5) $@
	$(MPIRUN) -n $(MPI_RANKS) -x PHASE2_LOG=$(PHASE2_LOG)	\
		$(PH2RUN) $(PH2RUN_FLAGS)

$(SIMUL_PROC): $(SIMUL_H5)
	$(TROTT_RPE) $(SIMUL_H5) | tee $@

run: $(SIMUL_H5)
analyze: $(SIMUL_PROC)

regen:
	$(PARSE_FCIDUMP) --sort-terms -o $(HAMIL_H5) $(FCIDUMP)
	$(PARSE_INPUTST) -o $(HAMIL_H5) $(INPUTST)

clean:
	$(RM) $(SIMUL_H5) $(SIMUL_PROC)
