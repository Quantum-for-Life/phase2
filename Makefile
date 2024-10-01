# --------------------------------------------------------------------------- #
# Configuration                                                               #
#                                                                             #
# Specify compile flags and the path to both MPI and HDF5 dynamic libraries   #
# and headers.                                                                #
# ----------------------------------------------------------------------------#
AS		:= nasm
ASFLAGS		+= -felf64 -w+all -w-reloc-rel-dword -Ox
CC		?= gcc
CFLAGS		+= -std=c17 -Wall -Wextra -O2 -march=native -mavx2
INCLUDE		:= ./include
LDFLAGS 	+=
LDLIBS		+= -lm
LIB64		:= /usr/lib/x86_64-linux-gnu

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

# --------------------------------------------------------------------------- #
# Build dependencies                                                          #
# --------------------------------------------------------------------------- #
PHASE2DIR	:= ./phase2
PH2RUNDIR	:= ./ph2run
LIBDIR		:= ./lib

# APIs
$(PHASE2DIR)/circ.o:	$(INCLUDE)/phase2/circ.h
$(PHASE2DIR)/data.o:	$(INCLUDE)/phase2/data.h
$(PHASE2DIR)/paulis.o:	$(INCLUDE)/phase2/paulis.h
$(PHASE2DIR)/qreg.o:	$(INCLUDE)/phase2/qreg.h
$(PHASE2DIR)/world.o:	$(INCLUDE)/phase2/world.h

$(LIBDIR)/xoshiro256ss.o:	$(INCLUDE)/xoshiro256ss.h

# Object files
PHASE2OBJS	:= $(PHASE2DIR)/circ.o					\
			$(PHASE2DIR)/circ_trott.o			\
			$(PHASE2DIR)/circ_qdrift.o			\
			$(PHASE2DIR)/data.o				\
			$(PHASE2DIR)/paulis.o				\
			$(PHASE2DIR)/qreg.o				\
			$(PHASE2DIR)/world.o
				
UTILSOBJS	:= $(LIBDIR)/xoshiro256ss.o

# Applications
PROGS		:= $(PH2RUNDIR)/ph2run-qdrift				\
			$(PH2RUNDIR)/ph2run-trott
$(PROGS):	$(PHASE2OBJS) $(UTILSOBJS)

# Update flags
CFLAGS		+= -I$(INCLUDE) $(MPI_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS		+= $(MPI_LDFLAGS) $(HDF5_LDFLAGS)
LDLIBS		+= $(MPI_LDLIBS) $(HDF5_LDLIBS)
# --------------------------------------------------------------------------- #
# Targets                                                                     #
# --------------------------------------------------------------------------- #
.DEFAULT_GOAL	:= all
.PHONY: all			\
	build build-test 	\
	clean			\
	debug			\
	check			\
	check-mpi

all: build build-test

debug: build build-test
debug: ASFLAGS	+= -DDEBUG -Og -Fdwarf
debug: CFLAGS	+= -DDEBUG -g -Og

build: $(PROGS)

clean:
	$(RM) $(PHASE2DIR)/*.o $(PHASE2DIR)/*.d
	$(RM) $(PH2RUNDIR)/*.o $(PH2RUNDIR)/*.d
	$(RM) $(LIBDIR)/*.o $(LIBDIR)/*.d
	$(RM) $(TESTDIR)/*.o $(TESTDIR)/*.d
	$(RM) $(PROGS)
	$(RM) $(TESTS)

# --------------------------------------------------------------------------- #
# Testing                                                                     #
# --------------------------------------------------------------------------- #
TESTDIR		:= ./test
CFLAGS		+= -I$(TESTDIR) -DPH2_TESTDIR=\"$(TESTDIR)\"

TESTS		:= $(TESTDIR)/t-circ_trott				\
			$(TESTDIR)/t-data_hamil				\
			$(TESTDIR)/t-data_multidet			\
			$(TESTDIR)/t-data_open				\
			$(TESTDIR)/t-data_trott_steps			\
			$(TESTDIR)/t-paulis				\
			$(TESTDIR)/t-qreg				\
			$(TESTDIR)/t-world

$(TESTS):	$(TESTDIR)/test.h					\
		$(TESTDIR)/t-data.h					\
		$(PHASE2OBJS) $(UTILSOBJS)

build-test: $(TESTS)

check: build-test
	@for tt in $(TESTS); do						\
		./$$tt && echo "$$tt: OK" || ( echo "$$tt: FAIL"; exit 1 ) ; \
	done

MPIRUN		:= mpirun
MPIFLAGS	:=
MPIRANKS	?= 2

check-mpi: build-test
	@for tt in $(TESTS); do						\
		$(MPIRUN) -n $(MPIRANKS) $(MPIFLAGS) ./$$tt && 		\
		echo "$$tt: OK" || ( echo "$$tt: FAIL"; exit 1 ) ; 	\
	done

