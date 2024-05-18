# ############################################################################
# Simulation parameters:
DELTA		= 1.0
EPSILON		= 0.001
QDIRFT_NUM_SAMPLES = 1
MPI_RANKS	= 4
# #############################################################################

# System
ROOT_DIR	= ../
SCRIPTS		= $(ROOT_DIR)/scripts
BUILD		= $(ROOT_DIR)/build

FCIDUMP		= FCIDUMP
INPUTST		= INPUTST
SIMUL_H5	= simul.h5
SIMUL_ORIG	= $(SIMUL_H5).orig
SIMUL_PROC	= $(SIMUL_H5).proc

PYTHON			= python3
PARSE_FCIDP		= $(PYTHON) $(SCRIPTS)/parse_fcidump.py
PARSE_INPUTST	= $(PYTHON) $(SCRIPTS)/parse_inputst.py
QDRIFT_SAMPLE	= $(PYTHON) $(SCRIPTS)/qdrift_sample.py
QDRIFT_RPE		= $(PYTHON) $(SCRIPTS)/qdrift_rpe.py
TROTT_VALID		= $(PYTHON) $(SCRIPTS)/trotter_validate.py
MPIRUN			= mpirun
PH2RUN			= $(BUILD)/ph2run
CP				= cp -f
RM				= rm -fv

# Application flags
PARSE_FCIDP_FLAGS	=
PARSE_INPUTST_FLAGS	=
QDRIFT_SAMPLE_FLAGS	= --delta=$(DELTA) --epsilon=$(EPSILON) --num-samples=$(QDIRFT_NUM_SAMPLES)
QDRIFT_RPE_FLAGS	= --delta=$(DELTA) --epsilon=$(EPSILON) --num-samples=$(QDIRFT_NUM_SAMPLES)
TROTT_VALID_FLAGS	=
MPIRUN_FLAGS		= -n $(MPI_RANKS)
PH2RUN_FLAGS		= 1
PR2RUN_LOGLVL		= info


.DEFAULT_GOAL := all
.PHONY: all clean simulate

all: simulate

$(SIMUL_ORIG): $(FCIDUMP) $(INPUTST)
	$(PARSE_FCIDP) $(PARSE_FCIDP_FLAGS) -o $@ $(FCIDUMP)
	$(PARSE_INPUTST) $(PARSE_INPUTST_FLAGS) -o $@ $(INPUTST)
	$(CP) $@ $(SIMUL_H5)
	$(QDRIFT_SAMPLE) $(QDRIFT_SAMPLE_FLAGS) $(SIMUL_H5)

simulate: $(SIMUL_ORIG)
	find . -iname "$(SIMUL_H5)-*" -exec $(MPIRUN) $(MPIRUN_FLAGS) -x PHASE2_LOG=$(PR2RUN_LOGLVL) $(PH2RUN) {} $(PH2RUN_FLAGS) \;

rpe:
	$(QDRIFT_RPE) $(QDRIFT_RPE_FLAGS) $(SIMUL_H5)
clean:
	$(RM) $(SIMUL_ORIG)* $(SIMUL_H5)* $(SIMUL_PROC)*
