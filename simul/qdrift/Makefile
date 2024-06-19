# -----------------------------------------------------------------------------
# Simulation parameters:
DELTA		= 1.0
EPSILON		= 0.01
NUM_SAMPLES 	= 256
MPI_RANKS	= 4
# -----------------------------------------------------------------------------

# System
ROOT_DIR	= ../..
SCRIPTS		= $(ROOT_DIR)/scripts
BUILD		= $(ROOT_DIR)/ph2run

FCIDUMP		?= FCIDUMP
INPUTST		?= INPUTST
SIMUL_H5	:= simul.h5
SIMUL_ORIG	:= $(SIMUL_H5).orig

PYTHON		?= python3
MPIRUN		?= mpirun
PH2RUN		?= $(BUILD)/ph2run-qdrift

PARSE_FCIDUMP	= $(PYTHON) $(SCRIPTS)/parse_fcidump.py
PARSE_INPUTST	= $(PYTHON) $(SCRIPTS)/parse_inputst.py
QDRIFT_SAMPLE	= $(PYTHON) $(SCRIPTS)/qdrift_sample.py
QDRIFT_RPE	= $(PYTHON) $(SCRIPTS)/qdrift_rpe.py

# Application flags
PARSE_FCIDUMP_FLAGS	=
PARSE_INPUTST_FLAGS	=
QDRIFT_SAMPLE_FLAGS	= --delta=$(DELTA) --epsilon=$(EPSILON) \
				--num-samples=$(NUM_SAMPLES)
QDRIFT_RPE_FLAGS	= --delta=$(DELTA) --epsilon=$(EPSILON) \
				--num-samples=$(NUM_SAMPLES)
MPIRUN_FLAGS		= -n $(MPI_RANKS)
PH2RUN_FLAGS		=
PR2RUN_LOGLVL		= info

.DEFAULT_GOAL := all
.PHONY: all clean simulate

all: simulate

$(SIMUL_ORIG): $(FCIDUMP) $(INPUTST)
	$(PARSE_FCIDUMP) $(PARSE_FCIDUMP_FLAGS) -o $@ $(FCIDUMP)
	$(PARSE_INPUTST) $(PARSE_INPUTST_FLAGS) -o $@ $(INPUTST)

$(SIMUL_H5): $(SIMUL_ORIG)
	cp $(SIMUL_ORIG) $(SIMUL_H5)
	$(QDRIFT_SAMPLE) $(QDRIFT_SAMPLE_FLAGS) $(SIMUL_H5)

simulate: $(SIMUL_H5)
	find . -iname "$(SIMUL_H5)-*" \
		-exec $(MPIRUN) $(MPIRUN_FLAGS) \
			-x PHASE2_LOG=$(PR2RUN_LOGLVL) \
			$(PH2RUN) {} $(PH2RUN_FLAGS) \;

rpe:
	$(QDRIFT_RPE) $(QDRIFT_RPE_FLAGS) $(SIMUL_H5)

clean:
	$(RM) $(SIMUL_ORIG) $(SIMUL_H5)*
