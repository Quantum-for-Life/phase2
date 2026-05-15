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
CFLAGS		+= -std=c11 -Wall -Wextra -O3 -march=native -mavx2 -MMD -MP
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

# Out-of-tree build root.  Every .o and .d lands
# under $(BUILDDIR)/<srcdir>/, mirroring the source
# layout.  Final binaries still live next to their
# sources (ph2run/ph2run, test/t-paulis, ...).
BUILDDIR	:= ./build

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
BACKEND_OBJS	+= $(BUILDDIR)/phase2/qreg_qreg.o
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
BACKEND_OBJS	+= $(BUILDDIR)/phase2/qreg_cuda.o			\
		   	$(BUILDDIR)/phase2/qreg_cuda_lo.o		\
			$(BUILDDIR)/phase2/qreg_cuda_lo_dlink.o		\
			$(BUILDDIR)/phase2/world_cuda.o
BACKEND_CFLAGS	+= -I$(CUDA_INCLUDE)
BACKEND_LDFLAGS	+= -L$(CUDA_LIBDIR) -Wl,-rpath -Wl,$(CUDA_LIBDIR)
BACKEND_LDLIBS	+= -lcudart -lstdc++

NVCCFLAGS	+= $(MPI_CFLAGS) $(HDF5_CFLAGS)
$(BUILDDIR)/phase2/qreg_cuda_lo.o: $(PHASE2DIR)/qreg_cuda_lo.cu
	@$(MKDIR) $(@D)
	$(NVCC) $(NVCCFLAGS) $(BACKEND_CFLAGS) 				\
	       -I$(INCLUDE) -c $< -o $@

$(BUILDDIR)/phase2/qreg_cuda_lo_dlink.o: $(BUILDDIR)/phase2/qreg_cuda_lo.o
	$(NVCC) $(NVCCFLAGS) -dlink $< -o $@

endif

BACKEND_CFLAGS	+= -DPHASE2_BACKEND=$(BACKEND_N)


# Header dependencies for every compiled .c are tracked
# automatically via -MMD/-MP (see CFLAGS below); the
# emitted .d files are pulled in at the foot of this
# file.  The Makefile only declares structural deps:
# object-set unions, inter-target link prereqs, and
# generated-header prereqs (none today).

PHASE2OBJS	:= $(BUILDDIR)/phase2/circ.o				\
			$(BUILDDIR)/phase2/circ_cache.o			\
			$(BUILDDIR)/phase2/paulis.o			\
			$(BUILDDIR)/phase2/prob.o			\
			$(BUILDDIR)/phase2/qreg.o			\
			$(BUILDDIR)/phase2/state_prep_coeff.o		\
			$(BUILDDIR)/phase2/world.o			\
			$(BACKEND_OBJS)

PH2RUN_DATA_OBJS := $(BUILDDIR)/ph2run/data.o

CIRCOBJS	:= $(BUILDDIR)/circ/cmpsit.o				\
			$(BUILDDIR)/circ/qdrift.o			\
			$(BUILDDIR)/circ/trott.o			\
			$(BUILDDIR)/circ/trott2.o

LIBOBJS		:= $(BUILDDIR)/lib/combinations.o			\
			$(BUILDDIR)/lib/det_small.o			\
			$(BUILDDIR)/lib/log.o				\
			$(BUILDDIR)/lib/xoshiro256ss.o

# Manifest of every .c / .cu source that the build
# system knows how to compile.  Used by the
# `check-srcs-coverage` drift guard to fail the
# build when a new source file is added to a
# subsystem directory but not wired into an OBJS
# variable.
DECLARED_SRCS	:= $(patsubst $(BUILDDIR)/%.o,%.c,			\
			$(PHASE2OBJS) $(CIRCOBJS)			\
			$(LIBOBJS) $(PH2RUN_DATA_OBJS))			\
		   phase2/phase2_run.c					\
		   ph2run/ph2run.c					\
		   bench/bench.c					\
		   bench/b-paulis.c bench/b-qreg.c			\
		   phase2/qreg_cuda.c					\
		   phase2/qreg_cuda_lo.cu				\
		   phase2/world_cuda.c

# Generic compile rules.  Two variants:
#  - $(BUILDDIR)/%.o:  regular objects.
#  - $(BUILDDIR)/shared/%.o:  same sources rebuilt
#    with -fPIC on a disjoint path, so the shared
#    library does not race the regular build's
#    objects.
$(BUILDDIR)/%.o: %.c
	@$(MKDIR) $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/shared/%.o: %.c
	@$(MKDIR) $(@D)
	$(CC) $(CFLAGS) -fPIC -c $< -o $@


# Applications
PROGS		:=  $(PH2RUNDIR)/ph2run

$(PH2RUNDIR)/ph2run: $(BUILDDIR)/circ/trott.o				\
			$(BUILDDIR)/circ/trott2.o			\
			$(BUILDDIR)/circ/qdrift.o			\
			$(BUILDDIR)/circ/cmpsit.o

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
	build-bench		\
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

build: check-srcs-coverage $(PROGS)

# --------------------------------------------------------------------------- #
# Shared library (Python interface)                                           #
# --------------------------------------------------------------------------- #
# libphase2.so carries the pure compute surface only; HDF5 stays
# on the ph2run side, so the shared object links without HDF5
# and a Python caller can drive phase2 over ctypes without ever
# touching an HDF5 file.
SHARED_LDFLAGS	:= $(MPI_LDFLAGS) $(BACKEND_LDFLAGS)
SHARED_LDLIBS	:= $(MPI_LDLIBS) $(BACKEND_LDLIBS)

SHARED_OBJS	:= $(BUILDDIR)/shared/phase2/phase2_run.o		\
		   $(patsubst $(BUILDDIR)/%, $(BUILDDIR)/shared/%,	\
			$(PHASE2OBJS) $(LIBOBJS))

shared: $(SHARED_OBJS)
	$(CC) -shared -o libphase2.so $^ $(SHARED_LDFLAGS) $(SHARED_LDLIBS)

# --------------------------------------------------------------------------- #
# Benchmarks                                                                  #
# --------------------------------------------------------------------------- #
BENCHDIR	:= ./bench
CFLAGS		+= -I$(BENCHDIR)

BENCHES		:= $(BENCHDIR)/b-paulis					\
			$(BENCHDIR)/b-qreg

$(BENCHES):	$(BENCHDIR)/bench.h					\
		$(BUILDDIR)/bench/bench.o				\
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
RUNFIXDIR	:= $(TESTDIR)/run-fixtures
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
			$(TESTDIR)/t-run				\
			$(TESTDIR)/t-state_prep_coeff_csf		\
			$(TESTDIR)/t-state_prep_coeff_expand		\
			$(TESTDIR)/t-world

# Synthetic fixtures used by t-run to drive the runner
# against controlled pass / fail / signal / banner
# outcomes.  Each fixture is a tiny self-contained C
# program with no phase2 / MPI / HDF5 deps.
RUNFIX		:= $(RUNFIXDIR)/pass					\
			$(RUNFIXDIR)/fail				\
			$(RUNFIXDIR)/sleep				\
			$(RUNFIXDIR)/abort				\
			$(RUNFIXDIR)/banner

$(RUNFIXDIR)/%: $(RUNFIXDIR)/%.c
	$(CC) -std=c11 -Wall -Wextra -O2 -o $@ $<

# t-run is a meta-test: it shells out to ./test/run, so
# both the runner and its fixtures must exist before
# t-run is invoked.  Use order-only prereqs (after `|`)
# so the implicit %: %.c link line does not try to feed
# the runner / fixture binaries to ld.
$(TESTDIR)/t-run: | $(TESTDIR)/run $(RUNFIX)

TESTS_SLOW	:= $(TESTDIR)/t-state_prep_coeff_large

$(TESTS) $(TESTS_SLOW):	$(TESTDIR)/test.h				\
			$(TESTDIR)/t-data.h				\
			$(PHASE2OBJS) $(LIBOBJS) $(PH2RUN_DATA_OBJS)

$(TESTDIR)/t-circ_cache: $(BUILDDIR)/circ/trott.o
$(TESTDIR)/t-circ_trott: $(BUILDDIR)/circ/trott.o
$(TESTDIR)/t-circ_trott2: $(BUILDDIR)/circ/trott2.o
$(TESTDIR)/t-circ_trott_coeff: $(BUILDDIR)/circ/trott.o
$(TESTDIR)/t-circ_trott2_coeff: $(BUILDDIR)/circ/trott2.o

build-test: check-tests-coverage $(TESTS) $(TESTDIR)/run

# Parallel cargo-style runner.  Standalone C binary, no
# phase2 / MPI / HDF5 dependencies.  Build with the test
# CFLAGS so any DEBUG / sanitiser flags apply consistently.
$(TESTDIR)/run: $(TESTDIR)/run.c
	$(CC) -std=c11 -Wall -Wextra -O2 -o $@ $<

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

# Drift guard for the subsystem source dirs (phase2,
# circ, lib, ph2run, bench).  Enumerate every .c /
# .cu in those directories and verify each appears
# in DECLARED_SRCS.  A new source dropped into one
# of them but not wired into an OBJS variable fails
# the build immediately.
.PHONY: check-srcs-coverage
check-srcs-coverage:
	@tmp=$$(mktemp -d);					\
	ls phase2/*.c phase2/*.cu				\
	   circ/*.c lib/*.c ph2run/*.c bench/*.c		\
	   2>/dev/null | sort > $$tmp/expected;			\
	for s in $(DECLARED_SRCS); do echo $$s; done		\
		| sort -u > $$tmp/declared;			\
	missing=$$(comm -23 $$tmp/expected $$tmp/declared);	\
	rc=0;							\
	if [ -n "$$missing" ]; then				\
		echo "build: .c sources not in DECLARED_SRCS:";	\
		echo "$$missing";				\
		rc=1;						\
	fi;							\
	rm -rf $$tmp;						\
	exit $$rc

# Tests are always built with -DDEBUG so log_trace /
# log_debug in include/log.h expand to real emits.
# t-log_release cancels this so the release-build strip
# of trace/debug macros can be verified end-to-end --
# the override must come after the blanket -DDEBUG so
# the trailing -UDEBUG wins on the gcc command line.
$(TESTS): CFLAGS += -DDEBUG -g -Og
$(TESTDIR)/t-log_release: CFLAGS += -UDEBUG

# All check targets delegate to test/run; the runner
# fans the suite out across cores, captures per-test
# stdout/stderr, and prints a cargo-style summary.
# See test/run.c for the harness contract.
.PHONY: check
check: build-test
	@./$(TESTDIR)/run -- $(TESTS)				\
		$(TESTDIR)/t-ref-coeff_matrix.py

# Slow tests: built but not part of the default check
# target.
.PHONY: build-test-slow
build-test-slow: $(TESTS_SLOW)

.PHONY: check-slow
check-slow: build-test-slow $(TESTDIR)/run
	@./$(TESTDIR)/run -- $(TESTS_SLOW)

.PHONY: test-slow
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

.PHONY: check-mpi
check-mpi: build-test
	@./$(TESTDIR)/run --mpiranks=$(MPIRANKS) -- $(TESTS)

# Pattern target: `make check-<filter>` runs the
# subset of the suite whose effective names match
# `<filter>*` (fnmatch glob).  The effective name is
# the basename without the `t-` prefix and the `.py`
# suffix.  Examples:
#   make check-data    -> all t-data_*
#   make check-paulis  -> just t-paulis
#   make check-circ    -> all t-circ*
# Explicit targets (check-mpi, check-slow,
# check-tests-coverage) win over this pattern.
check-%: build-test
	@./$(TESTDIR)/run --filter='$**' -- $(TESTS)		\
		$(TESTDIR)/t-ref-coeff_matrix.py


# --------------------------------------------------------------------------- #
# Clean up.                                                                   #
# --------------------------------------------------------------------------- #

clean:
	@$(RM) -r $(BUILDDIR)
	@$(RM) $(TESTDIR)/*.d

distclean: clean
	@$(RM) $(BENCHES)
	@$(RM) $(TESTS)
	@$(RM) $(TESTDIR)/run
	@$(RM) $(RUNFIX)
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


# --------------------------------------------------------------------------- #
# Auto-dep includes.  The -MMD -MP in CFLAGS makes gcc                        #
# emit a <obj>.d alongside every <obj>.o it produces.                         #
# Each .d names the .o target and lists every header                          #
# the .c #include'd, with -MP adding empty fallback                           #
# rules for each header so a renamed header doesn't                           #
# break the build.  The wildcard is silently empty on                         #
# the first build.                                                            #
# --------------------------------------------------------------------------- #

-include $(shell find $(BUILDDIR) -name '*.d' 2>/dev/null)
-include $(wildcard $(TESTDIR)/*.d)
