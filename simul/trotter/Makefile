# -----------------------------------------------------------------------------
# Simulation parameters:
DELTA		= 0.01
TROTT_STEPS	= 1024
MPI_RANKS	= 2
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
PH2RUN		?= $(BUILD)/ph2run-trott

PARSE_FCIDUMP	= $(PYTHON) $(SCRIPTS)/parse_fcidump.py
PARSE_INPUTST	= $(PYTHON) $(SCRIPTS)/parse_inputst.py
TROTT_PREP	= $(PYTHON) $(SCRIPTS)/trott_prep.py
TROTT_RPE	= $(PYTHON) $(SCRIPTS)/trott_rpe.py

# Application flags
PARSE_FCIDUMP_FLAGS	= --sort-terms
PARSE_INPUTST_FLAGS	=
TROTT_PREP_FLAGS	= --time-factor=$(TIME_FACT)
TROTT_RPE_FLAGS		=
MPIRUN_FLAGS		= -n $(MPI_RANKS)
PH2RUN_FLAGS		= --steps=$(TROTT_STEPS) --delta=$(DELTA)
PHASE2_LOG		?= info


.DEFAULT_GOAL: all
.PHONY: all clean

all: $(SIMUL_PROC)

$(SIMUL_ORIG): $(FCIDUMP) $(INPUTST)
	$(PARSE_FCIDUMP) $(PARSE_FCIDUMP_FLAGS) -o $@ $(FCIDUMP)
	$(PARSE_INPUTST) $(PARSE_INPUTST_FLAGS) -o $@ $(INPUTST)

$(SIMUL_H5): $(SIMUL_ORIG)
	cp $< $@
	$(MPIRUN) $(MPIRUN_FLAGS) -x PHASE2_LOG=$(PHASE2_LOG) \
		$(PH2RUN) $(PH2RUN_FLAGS) $(SIMUL_H5) 

$(SIMUL_PROC): $(SIMUL_H5)
	@echo -n "$(TIME_FACT),$(TROTT_STEPS)," >> $@
	$(TROTT_RPE) $(TROTT_RPE_FLAGS) $< | tee -a $@

clean:
	$(RM) $(SIMUL_ORIG) $(SIMUL_H5) $(SIMUL_PROC)
