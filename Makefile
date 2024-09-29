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
PH2DIR	=	./ph2run

# APIs
$(PH2DIR)/circ.o:		$(INCLUDE)/circ.h
$(PH2DIR)/data.o:		$(INCLUDE)/data.h
$(PH2DIR)/paulis.o:		$(INCLUDE)/paulis.h
$(PH2DIR)/qreg.o:		$(INCLUDE)/qreg.h
$(PH2DIR)/world.o:		$(INCLUDE)/world.h
$(PH2DIR)/xoshiro256ss.o:	$(INCLUDE)/xoshiro256ss.h

# Object files
PH2OBJS	=			$(PH2DIR)/circ.o			\
				$(PH2DIR)/data.o			\
				$(PH2DIR)/world.o			\
				$(PH2DIR)/paulis.o			\
				$(PH2DIR)/qreg.o			\
				$(PH2DIR)/xoshiro256ss.o

# Applications
$(PH2DIR)/ph2run-trott:		$(PH2OBJS) $(PH2DIR)/circ_trott.o
$(PH2DIR)/ph2run-qdrift:	$(PH2OBJS) $(PH2DIR)/circ_qdrift.o
PROGS		:= $(PH2DIR)/ph2run-qdrift				\
			$(PH2DIR)/ph2run-trott

# Update flags
CFLAGS		+= -I$(INCLUDE) -I$(PH2DIR)				\
			$(MPI_CFLAGS) $(HDF5_CFLAGS)
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
	check

all: build build-test

debug:	build build-test
debug:	ASFLAGS	+= -DDEBUG -Og -Fdwarf
debug:	CFLAGS	+= -DDEBUG -g -Og

build: $(PROGS)

clean:
	$(RM) $(PH2DIR)/*.o $(PH2DIR)/*.d
	$(RM) $(PROGS)
	$(RM) $(TESTDIR)/*.o $(TESTDIR)/*.d
	$(RM) $(TESTS)

# --------------------------------------------------------------------------- #
# Testing                                                                     #
# --------------------------------------------------------------------------- #
TESTDIR		:= ./test
CFLAGS		+= -I$(TESTDIR) -DPH2_TESTDIR=\"$(TESTDIR)\"

$(TESTDIR)/t-data_hamil:	$(PH2OBJS)				\
				$(TESTDIR)/test.h			\
				$(TESTDIR)/test-data.h
$(TESTDIR)/t-data_multidet:	$(PH2OBJS)				\
				$(TESTDIR)/test.h			\
				$(TESTDIR)/test-data.h
$(TESTDIR)/t-data_open:		$(PH2OBJS)				\
				$(TESTDIR)/test.h			\
				$(TESTDIR)/test-data.h
$(TESTDIR)/t-data_trott_steps:	$(PH2OBJS)				\
				$(TESTDIR)/test.h			\
				$(TESTDIR)/test-data.h
$(TESTDIR)/t-paulis:		$(TESTDIR)/test.h			\
				$(PH2DIR)/paulis.o			\
				$(PH2DIR)/xoshiro256ss.o
$(TESTDIR)/t-trott_caserand:	$(PH2OBJS)				\
				$(PH2DIR)/circ_trott.o			\
				$(TESTDIR)/test.h

TESTS			:=	$(TESTDIR)/t-data_hamil			\
				$(TESTDIR)/t-data_multidet		\
				$(TESTDIR)/t-data_open			\
				$(TESTDIR)/t-data_trott_steps		\
				$(TESTDIR)/t-paulis			\
				$(TESTDIR)/t-trott_caserand

build-test: $(TESTS)

check: build-test
	@for tt in $(TESTS); do						\
		./$$tt && echo "$$tt: OK" || ( echo "$$tt: FAIL"; exit 1 ) ; \
	done

