# --------------------------------------------------------------------------- #
# Configuration                                                               #
#                                                                             #
# Specify compile flags and the path to both MPI and HDF5 dynamic libraries   #
# and headers.                                                                #
# ----------------------------------------------------------------------------#
AS		:= nasm
ASFLAGS		+= -felf64 -w+all -w-reloc-rel-dword -Ox
CC		?= gcc
CFLAGS		+= -std=c11 -Wall -Wextra -O3 -march=native -mavx2
NVCC		?= nvcc
NVCCFLAGS	+= -O3
INCLUDE		:= ./include
LDFLAGS 	+=
LDLIBS		+= -lm
LIB64		:= /usr/lib/x86_64-linux-gnu

MPIRUN		:= mpirun
MPIFLAGS	:=
MPIRANKS	?= 2

RM		= rm -fv
MKDIR		= mkdir -p
# --------------------------------------------------------------------------- #
# Build dependencies                                                          #
# --------------------------------------------------------------------------- #
PHASE2DIR	:= ./phase2
PH2RUNDIR	:= ./ph2run
LIBDIR		:= ./lib

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
HDF5_LDLIBS	:= -lhdf5 -lhdf5_hl -lcurl -lsz -lz -ldl -lm

# Backends
#
# You can specify which backend the qreg API should use.  Available options
# are (case-sensitive):
#
# 	* qreg      - native engine
# 	* quest	    - QuEST library
#	* cuda      - NVIDIA GPU driver and runtime
#
# See below for how to specify the dependencies.
BACKEND		:= qreg
#BACKEND	:= quest
#BACKEND	:= cuda

BACKEND_OBJS	:=
BACKEND_CFLAGS	:=
BACKEND_LDFLAGS	:=
BACKEND_LDLIBS	:=

ifeq ($(BACKEND),qreg)
BACKEND_N	:= 0
BACKEND_OBJS	+= $(PHASE2DIR)/qreg.o
BACKEND_CFLAGS	+=
BACKEND_LDFLAGS	+=
BACKEND_LDLIBS	+=
endif

ifeq ($(BACKEND),quest)
# Specify QUEST_PREFIX, if you have QuEST installed in a nonstandard location.
QUEST_PREFIX	:=
QUEST_INCLUDE	:= $(QUEST_PREFIX)/usr/include/QuEST
QUEST_LIBDIR	:= $(QUEST_PREFIX)/usr/lib
BACKEND_N	:= 1
BACKEND_OBJS	+= $(PHASE2DIR)/qreg_quest.o				\
			$(PHASE2DIR)/world_quest.o
BACKEND_CFLAGS	+= -I$(QUEST_INCLUDE)
BACKEND_LDFLAGS	+= -L$(QUEST_LIBDIR) -Wl,-rpath -Wl,$(QUEST_LIBDIR)
BACKEND_LDLIBS	+= -lQuEST
$(BACKEND_OBJS): $(PHASE2DIR)/world_quest.h
endif

ifeq ($(BACKEND),cuda)
CUDA_PREFIX	:= /usr/local/cuda
CUDA_INCLUDE	:= $(CUDA_PREFIX)/include
CUDA_LIBDIR	:= $(CUDA_PREFIX)/lib64
BACKEND_N	:= 2
BACKEND_OBJS	+= $(PHASE2DIR)/qreg_cuda.o				\
		   	$(PHASE2DIR)/qreg_cuda_dlink.o			\
		   	$(PHASE2DIR)/qreg_cuda_target.o			\
			$(PHASE2DIR)/world_cuda.o
BACKEND_CFLAGS	+= -I$(CUDA_INCLUDE)
BACKEND_LDFLAGS	+= -L$(CUDA_LIBDIR) -Wl,-rpath -Wl,$(CUDA_LIBDIR)
BACKEND_LDLIBS	+= -lcudart -lstdc++

$(BACKEND_OBJS): $(PHASE2DIR)/qreg_cuda.h				\
       			$(PHASE2DIR)/world_cuda.h

$(PHASE2DIR)/qreg_cuda_target.o: $(PHASE2DIR)/qreg_cuda.cu
	$(NVCC) $(NVCCFLAGS) $(MPI_CFLAGS) $(BACKEND_CFLAGS) 		\
	       -I$(INCLUDE) -c $< -o $@

$(PHASE2DIR)/qreg_cuda_dlink.o: $(PHASE2DIR)/qreg_cuda_target.o
	$(NVCC) $(NVCCFLAGS) $< -o $@ -dlink

endif

BACKEND_CFLAGS	+= -DPHASE2_BACKEND=$(BACKEND_N)


# APIs
$(PHASE2DIR)/circ.o:	$(INCLUDE)/phase2/circ.h
$(PHASE2DIR)/data.o:	$(INCLUDE)/phase2/data.h
$(PHASE2DIR)/paulis.o:	$(INCLUDE)/phase2/paulis.h
$(PHASE2DIR)/qreg.o:	$(INCLUDE)/phase2/qreg.h
$(PHASE2DIR)/world.o:	$(INCLUDE)/phase2/world.h
$(PHASE2DIR)/world_log.o: $(INCLUDE)/phase2/world.h

$(LIBDIR)/xoshiro256ss.o:	$(INCLUDE)/xoshiro256ss.h

# Object files
PHASE2OBJS	:= $(PHASE2DIR)/circ.o					\
			$(PHASE2DIR)/circ_trott.o			\
			$(PHASE2DIR)/circ_qdrift.o			\
			$(PHASE2DIR)/data.o				\
			$(PHASE2DIR)/paulis.o				\
			$(BACKEND_OBJS)					\
			$(PHASE2DIR)/world.o				\
			$(PHASE2DIR)/world_log.o

UTILSOBJS	:= $(LIBDIR)/xoshiro256ss.o

# Applications
PROGS		:= $(PH2RUNDIR)/ph2run-qdrift				\
			$(PH2RUNDIR)/ph2run-trott
$(PROGS):	$(PHASE2OBJS) $(UTILSOBJS)

# Update flags
CFLAGS		+= -I$(INCLUDE)						\
			$(MPI_CFLAGS)					\
			$(HDF5_CFLAGS)					\
			$(BACKEND_CFLAGS)
LDFLAGS		+= $(MPI_LDFLAGS)					\
			$(HDF5_LDFLAGS)					\
			$(BACKEND_LDFLAGS)
LDLIBS		+= $(MPI_LDLIBS)					\
			$(HDF5_LDLIBS)					\
			$(BACKEND_LDLIBS)
# --------------------------------------------------------------------------- #
# Targets                                                                     #
# --------------------------------------------------------------------------- #
.DEFAULT_GOAL	:= all
.PHONY: all			\
	bench			\
	build			\
	bulid-bench		\
	build-test 		\
	clean			\
	debug			\
	check			\
	check-mpi		\
	format

all: build build-bench build-test

debug: build build-bench build-test
debug: ASFLAGS	+= -DDEBUG -Og -Fdwarf
debug: CFLAGS	+= -DDEBUG -g -Og

build: $(PROGS)

clean:
	$(RM) $(PHASE2DIR)/*.o $(PHASE2DIR)/*.d
	$(RM) $(PH2RUNDIR)/*.o $(PH2RUNDIR)/*.d
	$(RM) $(LIBDIR)/*.o $(LIBDIR)/*.d
	$(RM) $(PROGS)
	$(RM) $(BENCHDIR)/*.o $(BENCHDIR)/*.d
	$(RM) $(BENCHES)
	$(RM) $(TESTDIR)/*.o $(TESTDIR)/*.d
	$(RM) $(TESTS)

format:
	@find ./ -name "*.c" 						\
		-or -name "*.h"						\
		-or -name "*.cpp"					\
		-or -name "*.cu" | 					\
		while read f ; do					\
			clang-format --style=file -i $$f ;		\
		done
# --------------------------------------------------------------------------- #
# Benchmarks                                                                  #
# --------------------------------------------------------------------------- #
BENCHDIR	:= ./bench
CFLAGS		+= -I$(BENCHDIR)

BENCHES		:= $(BENCHDIR)/b-paulis					\
			$(BENCHDIR)/b-qreg

$(BENCHES):	$(BENCHDIR)/bench.h					\
		$(BENCHDIR)/bench.o					\
		$(PHASE2OBJS) $(UTILSOBJS)

build-bench: $(BENCHES)

bench: build-bench
	@for bb in $(BENCHES); do					\
		./$$bb &&						\
			echo "$$bb: OK" ||				\
 			( echo "$$bb: FAIL"; exit 1 );			\
	done

bench-mpi: build-bench
	@for bb in $(BENCHES); do					\
		$(MPIRUN) -n $(MPIRANKS) $(MPIFLAGS) ./$$bb && 		\
			echo "$$bb: OK" ||				\
			( echo "$$bb: FAIL"; exit 1 );	 		\
	done

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
		./$$tt &&						\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

check-mpi: build-test
	@for tt in $(TESTS); do						\
		$(MPIRUN) -n $(MPIRANKS) $(MPIFLAGS) ./$$tt && 		\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

