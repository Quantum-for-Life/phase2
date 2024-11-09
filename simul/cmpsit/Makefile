# -----------------------------------------------------------------------------
# Simulation parameters:
DEPTH		:= 16
LENGTH		:= 4
SAMPLES 	:= 32
STEP_SIZE	:= 0.0625
STEPS		:= 8
MPI_RANKS	?= 2
# -----------------------------------------------------------------------------

# System
ROOT_DIR	= ../..
SCRIPTS		= $(ROOT_DIR)/scripts
BUILD		= $(ROOT_DIR)/ph2run

FCIDUMP		?= FCIDUMP
INPUTST		?= INPUTST
SIMUL_H5	:= simul.h5
SIMUL_ORIG	:= $(SIMUL_H5).orig
SIMUL_PROC	:= $(SIMUL_H5).proc

PYTHON		?= python3
MPIRUN		?= mpirun
PH2RUN		?= $(BUILD)/ph2run-cmpsit

PARSE_FCIDUMP	= $(PYTHON) $(SCRIPTS)/parse_fcidump.py
PARSE_INPUTST	= $(PYTHON) $(SCRIPTS)/parse_inputst.py

# Application flags
PARSE_FCIDUMP_FLAGS	= --sort-terms
PARSE_INPUTST_FLAGS	=
MPIRUN_FLAGS		= -n $(MPI_RANKS)
PH2RUN_FLAGS		= --depth=$(DEPTH) 			\
				--length=$(LENGTH)		\
				--samples=$(SAMPLES)		\
				--step-size=$(STEP_SIZE)	\
				--steps=$(STEPS)
PHASE2_LOG		?= info

.DEFAULT_GOAL := all
.PHONY: all clean rpe

all: $(SIMUL_H5)

$(SIMUL_ORIG): $(FCIDUMP) $(INPUTST)
	$(PARSE_FCIDUMP) $(PARSE_FCIDUMP_FLAGS) -o $@ $(FCIDUMP)
	$(PARSE_INPUTST) $(PARSE_INPUTST_FLAGS) -o $@ $(INPUTST)

$(SIMUL_H5): $(SIMUL_ORIG)
	cp $< $@
	$(MPIRUN) $(MPIRUN_FLAGS) -x PHASE2_LOG=$(PHASE2_LOG)	\
		$(PH2RUN) $(PH2RUN_FLAGS) $@

clean:
	$(RM) $(SIMUL_ORIG) $(SIMUL_H5)*
