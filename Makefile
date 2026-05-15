# --------------------------------------------------------------------------- #
# Configuration                                                               #
#                                                                             #
# Specify compile flags and the path to both MPI and HDF5 dynamic libraries   #
# and headers.                                                                #
# ----------------------------------------------------------------------------#
VERSION_MAJOR	:= 1
VERSION_MINOR	:= 1
VERSION_PATCH	:= 0

AS		:= nasm
ASFLAGS		+= -felf64 -w+all -w-reloc-rel-dword -Ox
CC		?= gcc
CFLAGS		+= -std=c11 -Wall -Wextra -O3 -march=native -mavx2
# EXTRA_CFLAGS / EXTRA_LDFLAGS allow command-line overrides
# (e.g. for sanitizer builds) without clobbering the rest of
# the build flag pipeline.
EXTRA_CFLAGS	?=
EXTRA_LDFLAGS	?=
CFLAGS		+= $(EXTRA_CFLAGS)
INCLUDE		:= ./include
LDFLAGS 	+= $(EXTRA_LDFLAGS)
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
CIRCDIR		:= ./circ
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

# Standard (serial) HDF5.  Find the correct paths via:
#
# $ h5cc -shlib -show
#
HDF5_CFLAGS	:= -I/usr/include/hdf5/serial
HDF5_LDFLAGS	:= -L$(LIB64)/hdf5/serial -Wl,-rpath -Wl,$(LIB64)/hdf5/serial
HDF5_LDLIBS	:= -lhdf5 -lhdf5_hl -lcurl -lsz -lz -ldl -lm

# Backends
#
# You can specify which backend the qreg API should use.  Available options
# are (case-sensitive):
#
# 	* qreg      - native engine
#	* cuda      - NVIDIA GPU driver and runtime
#
# See below for how to specify the dependencies.
BACKEND		:= qreg
#BACKEND	:= cuda

BACKEND_OBJS	:=
BACKEND_CFLAGS	:=
BACKEND_LDFLAGS	:=
BACKEND_LDLIBS	:=

ifeq ($(BACKEND),qreg)
BACKEND_N	:= 0
BACKEND_OBJS	+= $(PHASE2DIR)/qreg_qreg.o
BACKEND_CFLAGS	+=
BACKEND_LDFLAGS	+=
BACKEND_LDLIBS	+=
endif

ifeq ($(BACKEND),cuda)
NVCC		?= nvcc
NVCCFLAGS	+= -O3 -dopt=on -arch=native
## Compile for NVIDIA H100
## https://docs.nvidia.com/cuda/hopper-compatibility-guide/index.html
#NVCCFLAGS      += -O3 -dopt=on -arch=sm_90a
CUDA_PREFIX	:=/usr/local/cuda
CUDA_INCLUDE	:=$(CUDA_PREFIX)/include
CUDA_LIBDIR	:=$(CUDA_PREFIX)/lib64
BACKEND_N	:= 2
BACKEND_OBJS	+= $(PHASE2DIR)/qreg_cuda.o				\
		   	$(PHASE2DIR)/qreg_cuda_lo.o			\
			$(PHASE2DIR)/qreg_cuda_lo_dlink.o		\
			$(PHASE2DIR)/world_cuda.o
BACKEND_CFLAGS	+= -I$(CUDA_INCLUDE)
BACKEND_LDFLAGS	+= -L$(CUDA_LIBDIR) -Wl,-rpath -Wl,$(CUDA_LIBDIR)
BACKEND_LDLIBS	+= -lcudart -lstdc++

$(BACKEND_OBJS): $(PHASE2DIR)/qreg_cuda.h				\
       			$(PHASE2DIR)/world_cuda.h

NVCCFLAGS	+= $(MPI_CFLAGS) $(HDF5_CFLAGS)
$(PHASE2DIR)/qreg_cuda_lo.o: $(PHASE2DIR)/qreg_cuda_lo.cu
	$(NVCC) $(NVCCFLAGS) $(BACKEND_CFLAGS) 				\
	       -I$(INCLUDE) -c $< -o $@

$(PHASE2DIR)/qreg_cuda_lo_dlink.o: $(PHASE2DIR)/qreg_cuda_lo.o
	$(NVCC) $(NVCCFLAGS) -dlink $< -o $@

endif

$(BACKEND_OBJS): $(PHASE2DIR)/qreg.h

BACKEND_CFLAGS	+= -DPHASE2_BACKEND=$(BACKEND_N)


# phase2 public API
$(PHASE2DIR)/circ.o:	$(INCLUDE)/phase2/circ.h			\
			$(INCLUDE)/phase2/state_prep_coeff.h
$(PH2RUNDIR)/data.o:	$(INCLUDE)/phase2/circ.h $(INCLUDE)/ph2run/data.h	\
			$(INCLUDE)/phase2/paulis.h
$(PHASE2DIR)/paulis.o:	$(INCLUDE)/phase2/paulis.h
$(PHASE2DIR)/prob.o:	$(INCLUDE)/phase2/prob.h
$(PHASE2DIR)/qreg.o:	$(INCLUDE)/phase2/qreg.h $(PHASE2DIR)/qreg.h
$(PHASE2DIR)/state_prep_coeff.o:	$(INCLUDE)/phase2/state_prep_coeff.h	\
			$(INCLUDE)/combinations.h			\
			$(INCLUDE)/det_small.h				\
			$(INCLUDE)/phase2/circ.h			\
			$(INCLUDE)/phase2/qreg.h			\
			$(PHASE2DIR)/qreg.h
$(PHASE2DIR)/world.o:	$(INCLUDE)/phase2/world.h

# internal API
$(PHASE2DIR)/circ_cache.o:	$(PHASE2DIR)/circ_cache.h
$(PHASE2DIR)/phase2_run.o:	$(INCLUDE)/phase2/phase2_run.h		\
				$(PHASE2DIR)/circ_cache.h

PHASE2OBJS	:= $(PHASE2DIR)/circ.o					\
			$(PHASE2DIR)/circ_cache.o			\
			$(PHASE2DIR)/paulis.o				\
			$(PHASE2DIR)/prob.o				\
			$(PHASE2DIR)/qreg.o				\
			$(PHASE2DIR)/state_prep_coeff.o			\
			$(PHASE2DIR)/world.o				\
			$(BACKEND_OBJS)

$(PHASE2OBJS):	$(INCLUDE)/phase2.h

PH2RUN_DATA_OBJS := $(PH2RUNDIR)/data.o


# Circuits
$(CIRCDIR)/cmpsit.o: $(INCLUDE)/circ/cmpsit.h
$(CIRCDIR)/qdrift.o: $(INCLUDE)/circ/qdrift.h
$(CIRCDIR)/trott.o: $(INCLUDE)/circ/trott.h
$(CIRCDIR)/trott2.o: $(INCLUDE)/circ/trott2.h

CIRCOBJS	:= $(CIRCDIR)/cmpsit.o					\
			$(CIRCDIR)/qdrift.o				\
			$(CIRCDIR)/trott.o				\
			$(CIRCDIR)/trott2.o


# Library / utilities
$(LIBDIR)/combinations.o: $(INCLUDE)/combinations.h
$(LIBDIR)/det_small.o:	$(INCLUDE)/det_small.h
$(LIBDIR)/log.o:	$(INCLUDE)/log.h
$(LIBDIR)/xoshiro256ss.o: $(INCLUDE)/xoshiro256ss.h

LIBOBJS		:= $(LIBDIR)/combinations.o				\
			$(LIBDIR)/det_small.o				\
			$(LIBDIR)/log.o					\
			$(LIBDIR)/xoshiro256ss.o


# Applications
PROGS		:=  $(PH2RUNDIR)/ph2run

$(PH2RUNDIR)/ph2run: $(CIRCDIR)/trott.o					\
			$(CIRCDIR)/trott2.o				\
			$(CIRCDIR)/qdrift.o				\
			$(CIRCDIR)/cmpsit.o

$(PROGS):	$(PHASE2OBJS)						\
			$(LIBOBJS)					\
			$(PH2RUN_DATA_OBJS)

# Update flags
VERSION		:= $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)
CFLAGS		+= -I$(INCLUDE)						\
			$(MPI_CFLAGS)					\
			$(HDF5_CFLAGS)					\
			$(BACKEND_CFLAGS)				\
			-DPHASE2_VERSION=\"$(VERSION)\"
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
	check			\
	check-mpi		\
	debug			\
	distclean		\
	format			\
	shared

all: build build-bench build-test

debug: build build-bench build-test
debug: ASFLAGS	+= -DDEBUG -Og -Fdwarf
debug: CFLAGS	+= -DDEBUG -g -Og

build: $(PROGS)

# --------------------------------------------------------------------------- #
# Shared library (Python interface)                                           #
# --------------------------------------------------------------------------- #
# libphase2.so carries the pure compute surface only; HDF5 stays
# on the ph2run side, so the shared object links without HDF5
# and a Python caller can drive phase2 over ctypes without ever
# touching an HDF5 file.
SHARED_LDFLAGS	:= $(MPI_LDFLAGS) $(BACKEND_LDFLAGS)
SHARED_LDLIBS	:= $(MPI_LDLIBS) $(BACKEND_LDLIBS)

shared: CFLAGS += -fPIC
shared: $(PHASE2DIR)/phase2_run.o $(PHASE2OBJS) $(LIBOBJS)
	$(CC) -shared -o libphase2.so $^ $(SHARED_LDFLAGS) $(SHARED_LDLIBS)

# --------------------------------------------------------------------------- #
# Benchmarks                                                                  #
# --------------------------------------------------------------------------- #
BENCHDIR	:= ./bench
CFLAGS		+= -I$(BENCHDIR)

BENCHES		:= $(BENCHDIR)/b-paulis					\
			$(BENCHDIR)/b-qreg

$(BENCHES):	$(BENCHDIR)/bench.h					\
		$(BENCHDIR)/bench.o					\
		$(PHASE2OBJS) $(LIBOBJS) $(PH2RUN_DATA_OBJS)

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
CFLAGS		+= -I$(TESTDIR) -I$(PHASE2DIR) -DPH2_TESTDIR=\"$(TESTDIR)\"

TESTS		:= $(TESTDIR)/t-bitstring_index			\
			$(TESTDIR)/t-circ_cache				\
			$(TESTDIR)/t-circ_prepst_coeff			\
			$(TESTDIR)/t-circ_trott				\
			$(TESTDIR)/t-circ_trott2			\
			$(TESTDIR)/t-circ_trott2_coeff			\
			$(TESTDIR)/t-circ_trott_coeff			\
			$(TESTDIR)/t-circ				\
			$(TESTDIR)/t-combinations			\
			$(TESTDIR)/t-data_attr				\
			$(TESTDIR)/t-data_coeff_matrix			\
			$(TESTDIR)/t-data_hamil				\
			$(TESTDIR)/t-data_hamil_validate		\
			$(TESTDIR)/t-data_multidet			\
			$(TESTDIR)/t-data_open				\
			$(TESTDIR)/t-data_dets_validate			\
			$(TESTDIR)/t-data_idempotence			\
			$(TESTDIR)/t-data_mpi				\
			$(TESTDIR)/t-data_trott_steps			\
			$(TESTDIR)/t-det_small				\
			$(TESTDIR)/t-log				\
			$(TESTDIR)/t-log_release			\
			$(TESTDIR)/t-paulis				\
			$(TESTDIR)/t-prob				\
			$(TESTDIR)/t-qreg				\
			$(TESTDIR)/t-ref-bendazzoli			\
			$(TESTDIR)/t-state_prep_coeff_csf		\
			$(TESTDIR)/t-state_prep_coeff_expand		\
			$(TESTDIR)/t-world

TESTS_SLOW	:= $(TESTDIR)/t-state_prep_coeff_large

$(TESTS) $(TESTS_SLOW):	$(TESTDIR)/test.h				\
			$(TESTDIR)/t-data.h				\
			$(PHASE2OBJS) $(LIBOBJS) $(PH2RUN_DATA_OBJS)

$(TESTDIR)/t-circ_cache: $(CIRCDIR)/trott.o
$(TESTDIR)/t-circ_trott: $(CIRCDIR)/trott.o
$(TESTDIR)/t-circ_trott2: $(CIRCDIR)/trott2.o
$(TESTDIR)/t-circ_trott_coeff: $(CIRCDIR)/trott.o
$(TESTDIR)/t-circ_trott2_coeff: $(CIRCDIR)/trott2.o

build-test: check-tests-coverage $(TESTS)

# Guard against a t-*.c file being added to test/ but
# forgotten in TESTS / TESTS_SLOW above -- a silent
# omission would otherwise produce a partial suite that
# CI cannot tell from a full one.
.PHONY: check-tests-coverage
check-tests-coverage:
	@tmp=$$(mktemp -d);					\
	ls $(TESTDIR)/t-*.c 2>/dev/null				\
		| sed 's|\.c$$||' | sort > $$tmp/expected;	\
	for t in $(TESTS) $(TESTS_SLOW); do echo $$t; done	\
		| sort > $$tmp/declared;			\
	missing=$$(comm -23 $$tmp/expected $$tmp/declared);	\
	rc=0;							\
	if [ -n "$$missing" ]; then				\
		echo "test/: t-*.c files not in TESTS:";	\
		echo "$$missing";				\
		rc=1;						\
	fi;							\
	rm -rf $$tmp;						\
	exit $$rc

CHECKS	:= $(TESTS:$(TESTDIR)/%=check/%)

.PHONY: $(CHECKS)
# Tests are always built with -DDEBUG so log_trace / log_debug
# in include/log.h expand to real emits.  t-log_release
# below cancels this so the release-build strip of
# trace/debug macros can be verified end-to-end -- the
# override must come after the blanket -DDEBUG so the
# trailing -UDEBUG wins on the gcc command line.
$(TESTS): CFLAGS += -DDEBUG -g -Og
$(CHECKS): CFLAGS += -DDEBUG -g -Og
$(TESTDIR)/t-log_release: CFLAGS += -UDEBUG
check/t-log_release: CFLAGS += -UDEBUG
$(CHECKS): check/%: $(TESTDIR)/%
	@./$< && echo "$< OK" || (echo "$<: FAIL"; exit 1)

.PHONY: check
check: $(CHECKS) check-python

# Python harness cross-validates the C expansion path against
# an independent in-tree reference oracle.  Depends on
# build-test (for the t-state_prep_coeff_expand binary).
.PHONY: check-python
check-python: build-test
	@python3 $(TESTDIR)/t-ref-coeff_matrix.py	\
		&& echo "$(TESTDIR)/t-ref-coeff_matrix.py OK"	\
		|| (echo "$(TESTDIR)/t-ref-coeff_matrix.py: FAIL"; exit 1)

# Slow tests: built but not part of the default check target.
.PHONY: build-test-slow
build-test-slow: $(TESTS_SLOW)

CHECKS_SLOW	:= $(TESTS_SLOW:$(TESTDIR)/%=check/%)
.PHONY: $(CHECKS_SLOW)
$(CHECKS_SLOW): CFLAGS += -DDEBUG -g -Og
$(CHECKS_SLOW): check/%: $(TESTDIR)/%
	@./$< && echo "$< OK" || (echo "$<: FAIL"; exit 1)
.PHONY: check-slow
check-slow: $(CHECKS_SLOW)
test-slow: check-slow

# Sanitizer / valgrind targets.  These rebuild from scratch
# under the requested instrumentation and re-run the full
# test suite.
.PHONY: test-asan test-valgrind test-mpi-asan
ASAN_FLAGS	:= -fsanitize=address,undefined -fno-omit-frame-pointer -g
# OpenMPI's startup path triggers internal "leaks" via
# libevent / libopen-pal that we can't fix from here.  Disable
# leak detection but keep all the other ASan + UBSan
# diagnostics active.
ASAN_OPTIONS_VAL := detect_leaks=0:halt_on_error=1

test-asan:
	$(MAKE) distclean
	ASAN_OPTIONS=$(ASAN_OPTIONS_VAL)			\
	$(MAKE) check						\
		EXTRA_CFLAGS="$(ASAN_FLAGS)"			\
		EXTRA_LDFLAGS="$(ASAN_FLAGS)"

test-valgrind: build-test
	@for tt in $(TESTS); do						\
		valgrind --quiet --error-exitcode=1			\
			--leak-check=full				\
			--errors-for-leak-kinds=all			\
			--track-origins=yes ./$$tt &&			\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

test-mpi-asan:
	$(MAKE) distclean
	$(MAKE) build-test					\
		EXTRA_CFLAGS="$(ASAN_FLAGS)"			\
		EXTRA_LDFLAGS="$(ASAN_FLAGS)"
	@for tt in $(TESTDIR)/t-circ_trott_coeff; do			\
		$(MPIRUN) -n 4 $(MPIFLAGS)				\
			-x ASAN_OPTIONS=$(ASAN_OPTIONS_VAL) ./$$tt &&	\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

#check: build-test
#	@for tt in $(TESTS); do						\
#		./$$tt &&						\
#			echo "$$tt: OK" ||				\
#			( echo "$$tt: FAIL"; exit 1 );			\
#	done

check-mpi: build-test
	@for tt in $(TESTS); do						\
		$(MPIRUN) -n $(MPIRANKS) $(MPIFLAGS) ./$$tt && 		\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done


# --------------------------------------------------------------------------- #
# Clean up.                                                                   #
# --------------------------------------------------------------------------- #

clean:
	@$(RM) $(CIRCDIR)/*.o $(CIRCDIR)/*.d
	@$(RM) $(BENCHDIR)/*.o $(BENCHDIR)/*.d
	@$(RM) $(LIBDIR)/*.o $(LIBDIR)/*.d
	@$(RM) $(PH2RUNDIR)/*.o $(PH2RUNDIR)/*.d
	@$(RM) $(PHASE2DIR)/*.o $(PHASE2DIR)/*.d
	@$(RM) $(TESTDIR)/*.o $(TESTDIR)/*.d

distclean: clean
	@$(RM) $(BENCHES)
	@$(RM) $(TESTS)
	@$(RM) $(PROGS)
	@$(RM) libphase2.so

format:
	@find ./ -name "*.c" 						\
		-or -name "*.h"						\
		-or -name "*.cpp"					\
		-or -name "*.cu" | 					\
		while read f ; do					\
			clang-format --style=file -i $$f ;		\
		done

