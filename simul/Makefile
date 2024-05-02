# ############################################################################
# Simulation parameters:
TIME_FACT	= 1.0
TROTT_STEPS	= 128
MPI_RANKS	= 4
# #############################################################################

# System
ROOT_DIR	= ..
SCRIPTS		= $(ROOT_DIR)/scripts
BUILD		= $(ROOT_DIR)/build

FCIDUMP		= FCIDUMP
SIMUL_H5	= simul.h5
SIMUL_ORIG	= $(SIMUL_H5).orig
SIMUL_PROC	= $(SIMUL_H5).proc

PYTHON		= python3
PARSE_FCIDP	= $(PYTHON) $(SCRIPTS)/parse_fcidump.py
TROTT_VALID	= $(PYTHON) $(SCRIPTS)/trotter_validate.py
MPIRUN		= mpirun
PH2RUN		= $(BUILD)/ph2run
CP		= cp -f
RM		= rm -fv

# Application flags
PARSE_FCIDP_FLAGS	= --sort-terms --time-factor=$(TIME_FACT)
TROTT_VALID_FLAGS	=
MPIRUN_FLAGS		= -n $(MPI_RANKS)
PH2RUN_FLAGS		= $(TROTT_STEPS)
PR2RUN_LOGLVL		= info


.DEFAULT_GOAL: all
.PHONY: all clean

all: $(SIMUL_PROC)

$(SIMUL_ORIG): FCIDUMP
	$(PARSE_FCIDP) $(PARSE_FCIDP_FLAGS) -o $@ $<

$(SIMUL_H5): $(SIMUL_ORIG)
	$(CP) $< $@
	$(MPIRUN) $(MPIRUN_FLAGS) -x PHASE2_LOG=$(PR2RUN_LOGLVL) \
 		$(PH2RUN) $(SIMUL_H5) $(PH2RUN_FLAGS)

$(SIMUL_PROC): $(SIMUL_H5)
	$(TROTT_VALID) $(TROTT_VALID_FLAGS) $< > $@
	@cat $@

clean:
	$(RM) $(SIMUL_ORIG) $(SIMUL_H5) $(SIMUL_PROC)