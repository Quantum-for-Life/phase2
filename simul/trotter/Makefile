# ############################################################################
# Simulation parameters:
SIMUL_TAG	= 0001
TIME_FACT	= 0.1
TROTT_STEPS	= 1024
MPI_RANKS	= 4
# #############################################################################

# System
ROOT_DIR	= ../..
SCRIPTS		= $(ROOT_DIR)/scripts
BUILD		= $(ROOT_DIR)/build

FCIDUMP		= FCIDUMP
INPUTST		= INPUTST
SIMUL_H5	= simul.h5
SIMUL_ORIG	= $(SIMUL_H5).orig
SIMUL_PROC	= $(SIMUL_H5).proc
SIMUL_ALL	= $(SIMUL_H5).all

PYTHON		= python3
PARSE_FCIDP	= $(PYTHON) $(SCRIPTS)/parse_fcidump.py
PARSE_INPUTST	= $(PYTHON) $(SCRIPTS)/parse_inputst.py
TROTT_PREP	= $(PYTHON) $(SCRIPTS)/trott_prep.py
TROTT_VALID	= $(PYTHON) $(SCRIPTS)/trott_validate.py
MPIRUN		= mpirun
PH2RUN		= $(BUILD)/ph2run-trott
CP			= cp -f
RM			= rm -fv

# Application flags
PARSE_FCIDP_FLAGS	= --sort-terms
PARSE_INPUTST_FLAGS =
TROTT_PREP_FLAGS	= --time-factor=$(TIME_FACT)
TROTT_VALID_FLAGS	=
MPIRUN_FLAGS		= -n $(MPI_RANKS)
PH2RUN_FLAGS		= $(TROTT_STEPS)
PR2RUN_LOGLVL		= info


.DEFAULT_GOAL: all
.PHONY: all clean dist-clean

all: $(SIMUL_PROC)

$(SIMUL_ORIG): $(FCIDUMP) $(INPUTST)
	$(PARSE_FCIDP) $(PARSE_FCIDP_FLAGS) -o $@ $(FCIDUMP)
	$(PARSE_INPUTST) $(PARSE_INPUTST_FLAGS) -o $@ $(INPUTST)
	$(TROTT_PREP) $(TROTT_PREP_FLAGS) -o $@

$(SIMUL_H5): $(SIMUL_ORIG)
	$(CP) $< $@
	$(MPIRUN) $(MPIRUN_FLAGS) -x PHASE2_LOG=$(PR2RUN_LOGLVL) \
 		$(PH2RUN) $(SIMUL_H5) $(PH2RUN_FLAGS)

$(SIMUL_PROC): $(SIMUL_H5)
	@echo -n "$(SIMUL_TAG),$(TIME_FACT),$(TROTT_STEPS)," > $@
	$(TROTT_VALID) $(TROTT_VALID_FLAGS) $< >> $@
	@cat $@ | tee -a $(SIMUL_ALL)

clean:
	$(RM) $(SIMUL_ORIG) $(SIMUL_H5) $(SIMUL_PROC)

dist-clean: clean
	$(RM) $(SIMUL_ALL)