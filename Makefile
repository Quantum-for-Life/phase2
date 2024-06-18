# --------------------------------------------------------------------------- #
# Configuration                                                               #
#                                                                             #
# Specify compile flags and the path to both MPI and HDF5 dynamic libraries   #
# and headers.                                                                #
# ----------------------------------------------------------------------------#

CC	?= gcc
INCLUDE	:= ./include
CFLAGS	+= -std=c17 -Wall -Wextra -O2 -march=native -I$(INCLUDE)
LDFLAGS +=
LDLIBS	+= -lm
LIB64	:= /usr/lib/x86_64-linux-gnu

# If you're unsure where to find the compiled MPI libraries or headers,
# but have OpenMPI installed in your system, you can query:
#
# $ mpicc -showme
#
MPI_CFLAGS	= -I$(LIB64)/openmpi/include
MPI_LDFLAGS	= -L$(LIB64)/openmpi/lib
MPI_LDLIBS	= -lmpi

# Make sure you use the _parallel_ version of HDF5.
# You can find the correct paths for your system by querying:
#
# $ h5pcc -shlib -show
#
HDF5_CFLAGS	= -I/usr/include/hdf5/openmpi
HDF5_LDFLAGS	= -L$(LIB64)/hdf5/openmpi -Wl,-rpath -Wl,$(LIB64)/hdf5/openmpi
HDF5_LDLIBS	= -lhdf5 -lhdf5_hl -lcrypto -lcurl -lsz -lz -ldl -lm

# Update flags
CFLAGS	+= $(MPI_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS	+= $(MPI_LDFLAGS) $(HDF5_LDFLAGS)
LDLIBS	+= $(MPI_LDLIBS) $(HDF5_LDLIBS)

# --------------------------------------------------------------------------- #
# Build dependencies                                                          #
# --------------------------------------------------------------------------- #

PH2RUNDIR	=	./ph2run
PH2RUNOBJS	=	$(PH2RUNDIR)/circ.o	\
			$(PH2RUNDIR)/data.o	\
			$(PH2RUNDIR)/log.o	\
			$(PH2RUNDIR)/qreg.o	\
			$(PH2RUNDIR)/xoshiro256ss.o

$(PH2RUNDIR)/circ.o:		$(INCLUDE)/circ.h
$(PH2RUNDIR)/data.o:		$(INCLUDE)/data.h
$(PH2RUNDIR)/log.o:		$(INCLUDE)/log.h
$(PH2RUNDIR)/qreg.o:		$(INCLUDE)/qreg.h
$(PH2RUNDIR)/xoshiro256ss.o:	$(INCLUDE)/xoshiro256ss.h

$(PH2RUNDIR)/ph2run-trott:	$(PH2RUNOBJS)	\
					$(PH2RUNDIR)/circ_trott.o
$(PH2RUNDIR)/ph2run-qdrift:	$(PH2RUNOBJS)	\
					$(PH2RUNDIR)/circ_qdrift.o

# --------------------------------------------------------------------------- #
# Targets                                                                     #
# --------------------------------------------------------------------------- #

PROGS	= 	$(PH2RUNDIR)/ph2run-qdrift \
		$(PH2RUNDIR)/ph2run-trott

.DEFAULT_GOAL := all
.PHONY: all			\
	build build-test 	\
	clean			\
	debug			\
	test

all: build build-test

debug:	build build-test
debug:	CFLAGS	+= -g -Og -DDEBUG

build: $(PROGS)

clean:
	$(RM) $(PH2RUNDIR)/*.o $(PH2RUNDIR)/*.d
	$(RM) $(PROGS)
	$(RM) $(TESTDIR)/*.o $(TESTDIR)/*.d
	$(RM) $(TESTS)

# --------------------------------------------------------------------------- #
# Testing                                                                     #
# --------------------------------------------------------------------------- #

TESTDIR	:= ./test
CFLAGS	+= -I$(TESTDIR) -DPH2_TESTDIR=\"$(TESTDIR)\"

$(TESTDIR)/test-data:	$(TESTDIR)/test-data_hamil.o		\
			$(TESTDIR)/test-data_multidet.o		\
			$(TESTDIR)/test-data_open.o		\
			$(TESTDIR)/test-data_trott_steps.o	\
				$(PH2RUNOBJS)

$(TESTDIR)/test-trott_caserand: $(PH2RUNOBJS) $(PH2RUNDIR)/circ_trott.o

TESTS	= 	$(TESTDIR)/test-data \
		$(TESTDIR)/test-trott_caserand

build-test: $(TESTS)

test: build-test
	@for tt in $(TESTS); do \
		$$tt && echo "$$tt: OK" || ( echo "$$tt: FAIL"; exit 1 ) ; \
	done

