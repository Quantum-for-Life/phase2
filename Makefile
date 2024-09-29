# --------------------------------------------------------------------------- #
# Configuration                                                               #
#                                                                             #
# Specify compile flags and the path to both MPI and HDF5 dynamic libraries   #
# and headers.                                                                #
# ----------------------------------------------------------------------------#
AS	:= nasm
ASFLAGS	+= -felf64 -w+all -w-reloc-rel-dword -Ox
CC	?= gcc
CFLAGS	+= -std=c17 -MP -MMD -Wall -Wextra -O2 -march=native -mavx2
INCLUDE	:= -I./include -I./ph2run
DEPS	:= $(wildcard *.d)
LDFLAGS +=
LDLIBS	+= -lm
LIB64	:= /usr/lib/x86_64-linux-gnu

# If you're unsure where to find the compiled MPI libraries or headers,
# but have OpenMPI installed in your system, you can query:
#
# $ mpicc -showme
#
MPI_CFLAGS	:= -I$(LIB64)/openmpi/include
MPI_LDFLAGS	:= -L$(LIB64)/openmpi/lib
MPI_LDLIBS	:= -lmpi

# Make sure you use the _parallel_ version of HDF5.
# You can find the correct paths for your system by querying:
#
# $ h5pcc -shlib -show
#
HDF5_CFLAGS	:= -I/usr/include/hdf5/openmpi
HDF5_LDFLAGS	:= -L$(LIB64)/hdf5/openmpi -Wl,-rpath -Wl,$(LIB64)/hdf5/openmpi
HDF5_LDLIBS	:= -lhdf5 -lhdf5_hl -lcrypto -lcurl -lsz -lz -ldl -lm

# Update flags
CFLAGS	+= $(INCLUDE) $(MPI_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS	+= $(MPI_LDFLAGS) $(HDF5_LDFLAGS)
LDLIBS	+= $(MPI_LDLIBS) $(HDF5_LDLIBS)

# --------------------------------------------------------------------------- #
# Build dependencies                                                          #
# --------------------------------------------------------------------------- #
PH2RUNDIR	=	./ph2run
PH2RUNOBJS	=	$(PH2RUNDIR)/circ.o	\
			$(PH2RUNDIR)/data.o	\
			$(PH2RUNDIR)/world.o	\
			$(PH2RUNDIR)/paulis.o	\
			$(PH2RUNDIR)/qreg.o	\
			$(PH2RUNDIR)/xoshiro256ss.o

$(PH2RUNDIR)/ph2run-trott:	$(PH2RUNOBJS)	\
					$(PH2RUNDIR)/circ_trott.o
$(PH2RUNDIR)/ph2run-qdrift:	$(PH2RUNOBJS)	\
					$(PH2RUNDIR)/circ_qdrift.o

# --------------------------------------------------------------------------- #
# Targets                                                                     #
# --------------------------------------------------------------------------- #
PROGS	:= 	$(PH2RUNDIR)/ph2run-qdrift \
		$(PH2RUNDIR)/ph2run-trott

ifneq ($(DEPS),)
include $(DEPS)
endif

.DEFAULT_GOAL := all
.PHONY: all			\
	build build-test 	\
	clean			\
	debug			\
	check

all: build build-test

debug:	build build-test
debug:	ASFLAGS	+= -DDEBUG -Og -Fdwarf
debug:	CFLAGS	+= -DDEBUG -g -Og

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

TESTS	:= 	$(TESTDIR)/test-data \
		$(TESTDIR)/test-trott_caserand

build-test: $(TESTS)

check: build-test
	@for tt in $(TESTS); do						\
		./$$tt && echo "$$tt: OK" || ( echo "$$tt: FAIL"; exit 1 ) ; \
	done

