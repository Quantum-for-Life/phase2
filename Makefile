# --------------------------------------------------------------------------- #
# Configuration

CC	?= gcc
CFLAGS	+= -Wall -Wextra -O3 -march=native
INCLUDE	+= -I./include
LDFLAGS += 
LDLIBS	+= -lm
LIB64	:= /usr/lib/x86_64-linux-gnu

# If you're unsure where to find the compiled MPI libraries 
# or headers, but have OpenMPI installed in your system,
# you can query:
#
#	mpicc --showme
#
MPI_INCLUDE	= -I$(LIB64)/openmpi/include
MPI_LDFLAGS	= -L$(LIB64)/openmpi/lib
MPI_LDLIBS	= -lmpi

# Make sure all paths point to the _parallel_ version of 
# HDF5.  You can the correct paths for your system by querying:
#
#	h5pcc --show
#
HDF5_INCLUDE	= -I/usr/include/hdf5/openmpi
HDF5_LDFLAGS	= -Wl,-rpath -Wl,$(LIB64)/hdf5/openmpi
HDF5_LDLIBS	= -lcrypto -lcurl -lsz -lz -ldl -lm
HDF5_LIBS	= $(LIB64)/hdf5/openmpi/libhdf5_hl.a \
		 	$(LIB64)/hdf5/openmpi/libhdf5.a 
# Update flags
INCLUDE	+= $(MPI_INCLUDE) $(HDF5_INCLUDE)
LDFLAGS	+= $(MPI_LDFLAGS) $(HDF5_LDFLAGS)
LDLIBS	+= $(MPI_LDLIBS) $(HDF5_LIBS) $(HDF5_LDLIBS)
CFLAGS 	+= $(INCLUDE)

# --------------------------------------------------------------------------- #
# Build dependencies
PH2RUNDIR 	= ./ph2run
PH2RUNOBJS	= \
	$(PH2RUNDIR)/circ.o 		\
	$(PH2RUNDIR)/circ_qdrift.o 	\
	$(PH2RUNDIR)/circ_trott.o 	\
	$(PH2RUNDIR)/data.o 		\
	$(PH2RUNDIR)/log.o		\
	$(PH2RUNDIR)/qreg.o 		\
	$(PH2RUNDIR)/xoshiro256ss.o

$(PH2RUNDIR)/ph2run-trott: $(PH2RUNOBJS)
$(PH2RUNDIR)/ph2run-qdrift: $(PH2RUNOBJS)

# --------------------------------------------------------------------------- #
# Targets

PROGS	= $(PH2RUNDIR)/ph2run-qdrift \
		$(PH2RUNDIR)/ph2run-trott

.DEFAULT_GOAL := all
.PHONY: all			\
	build build-test 	\
	clean			\
	debug			\
	test

all: build

debug:	build build-test
debug:	CFLAGS	+= -g -Og -DDEBUG

build: $(PROGS)

# --------------------------------------------------------------------------- #
# Testing

TESTDIR	= ./test

CFLAGS	+= -I$(TESTDIR) -DPH2_TESTDIR=\"$(TESTDIR)\"

$(TESTDIR)/test-data: $(TESTDIR)/test-data_hamil.o	\
		$(TESTDIR)/test-data_multidet.o		\
		$(TESTDIR)/test-data_open.o		\
		$(TESTDIR)/test-data_trott_steps.o	\
			$(PH2RUNOBJS)

$(TESTDIR)/test-trott_caserand: $(PH2RUNOBJS)

TESTS	= $(TESTDIR)/test-data \
		$(TESTDIR)/test-trott_caserand

build-test: $(TESTS)

test: build-test
	@for tt in $(TESTS); do \
		echo \--- $$tt ; \
		$$tt && echo OK || ( echo FAIL; exit 1 ) ; \
	done 

clean:
	$(RM) $(PH2RUNDIR)/*.o $(PH2RUNDIR)/*.d
	$(RM) $(PROGS)
	$(RM) $(TESTDIR)/*.o $(TESTDIR)/*.d
	$(RM) $(TESTS)

