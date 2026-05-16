# ========================================================================== #
#                                                                            #
#   phase2 -- top-level build orchestrator                                   #
#                                                                            #
#   Three sections below:                                                    #
#     1. USER CONFIGURATION -- toolchain, library paths, backend.            #
#     2. INTERNAL           -- subsystem layout, derived vars, exports.      #
#     3. TARGETS            -- dispatch + high-level user-facing targets.    #
#                                                                            #
#   Per-subsystem build rules live in $(SUBSYS)/Makefile.  Shared compile    #
#   patterns live in build-rules.mk.  A user porting phase2 to a new        #
#   machine edits Section 1 only.                                           #
#                                                                            #
# ========================================================================== #


# ========================================================================== #
# Section 1: USER CONFIGURATION                                              #
# ========================================================================== #
# Edit values below to match your toolchain and library install paths.       #
# No edits should be needed elsewhere in this Makefile (or in the           #
# per-subsystem Makefiles) for a typical port.                              #
# -------------------------------------------------------------------------- #

# --- Toolchain ------------------------------------------------------------- #
CC              ?= gcc
CFLAGS          += -std=c11 -Wall -Wextra -O3 -march=native -mavx2 -MMD -MP
# EXTRA_CFLAGS / EXTRA_LDFLAGS allow command-line overrides (e.g. sanitiser
# builds) without clobbering the rest of the flag pipeline.
EXTRA_CFLAGS    ?=
EXTRA_LDFLAGS   ?=

# --- System library root --------------------------------------------------- #
# Most Debian/Ubuntu installs put MPI and HDF5 dynamic libs under
# /usr/lib/x86_64-linux-gnu.  Adjust if your distro differs.
LIB64           := /usr/lib/x86_64-linux-gnu

# --- MPI ------------------------------------------------------------------- #
# Default: OpenMPI from the system package.  Query paths via `mpicc -showme`
# if your install lives elsewhere.
MPI_CFLAGS      := -I$(LIB64)/openmpi/include
MPI_LDFLAGS     := -L$(LIB64)/openmpi/lib
MPI_LDLIBS      := -lmpi
MPIRUN          := mpirun
MPIFLAGS        :=
MPIRANKS        ?= 2

# --- HDF5 (serial) --------------------------------------------------------- #
# Query paths via `h5cc -shlib -show` if unsure.
HDF5_CFLAGS     := -I/usr/include/hdf5/serial
HDF5_LDFLAGS    := -L$(LIB64)/hdf5/serial -Wl,-rpath -Wl,$(LIB64)/hdf5/serial
HDF5_LDLIBS     := -lhdf5 -lhdf5_hl -lcurl -lsz -lz -ldl -lm

# --- Backend (qreg = CPU, cuda = NVIDIA GPU) ------------------------------- #
BACKEND         := qreg
# CUDA toolkit path (only relevant when BACKEND=cuda).
CUDA_PREFIX     ?= /usr/local/cuda
NVCC            ?= nvcc
NVCCFLAGS       += -O3 -dopt=on -arch=native
# Compile for NVIDIA H100 instead of the build host:
#NVCCFLAGS       += -O3 -dopt=on -arch=sm_90a

# --- Build output ---------------------------------------------------------- #
# Every .o and .d lands under $(BUILDDIR)/<srcdir>/, mirroring the source
# tree.  Final binaries stay next to their sources (ph2run/ph2run, ...).
BUILDDIR        := $(CURDIR)/build

# --- Version --------------------------------------------------------------- #
VERSION_MAJOR   := 1
VERSION_MINOR   := 1
VERSION_PATCH   := 0
VERSION         := $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)


# ========================================================================== #
# Section 2: INTERNAL                                                        #
# ========================================================================== #
# Derived variables, subsystem layout, and the export contract that makes   #
# values from Section 1 visible to per-subsystem Makefiles.  Edit only when #
# adding or removing a subsystem.                                           #
# -------------------------------------------------------------------------- #

# --- Tree layout ----------------------------------------------------------- #
TOPDIR          := $(CURDIR)
INCLUDE         := $(TOPDIR)/include
CIRCDIR         := $(TOPDIR)/circ
PHASE2DIR       := $(TOPDIR)/phase2
LIBDIR          := $(TOPDIR)/lib
PH2RUNDIR       := $(TOPDIR)/ph2run
BENCHDIR        := $(TOPDIR)/bench
TESTDIR         := $(TOPDIR)/test

# --- Backend-conditional objects + flags ----------------------------------- #
BACKEND_OBJS    :=
BACKEND_CFLAGS  :=
BACKEND_LDFLAGS :=
BACKEND_LDLIBS  :=

ifeq ($(BACKEND),qreg)
BACKEND_N       := 0
BACKEND_OBJS    += $(BUILDDIR)/phase2/qreg_qreg.o
endif

ifeq ($(BACKEND),cuda)
BACKEND_N       := 2
CUDA_INCLUDE    := $(CUDA_PREFIX)/include
CUDA_LIBDIR     := $(CUDA_PREFIX)/lib64
BACKEND_OBJS    += $(BUILDDIR)/phase2/qreg_cuda.o		\
                   $(BUILDDIR)/phase2/qreg_cuda_lo.o		\
                   $(BUILDDIR)/phase2/qreg_cuda_lo_dlink.o	\
                   $(BUILDDIR)/phase2/world_cuda.o
BACKEND_CFLAGS  += -I$(CUDA_INCLUDE)
BACKEND_LDFLAGS += -L$(CUDA_LIBDIR) -Wl,-rpath -Wl,$(CUDA_LIBDIR)
BACKEND_LDLIBS  += -lcudart -lstdc++
NVCCFLAGS       += $(MPI_CFLAGS) $(HDF5_CFLAGS)
endif

BACKEND_CFLAGS  += -DPHASE2_BACKEND=$(BACKEND_N)

# --- Object set unions ----------------------------------------------------- #
# Each subsys's Makefile derives its own OBJS from its SRCS list; these
# union aliases let downstream subsystems (ph2run, bench, test) reference
# all upstream objects in one place.  Kept in sync with the SRCS lists in
# each subsys Makefile by the per-subsys check-srcs-coverage drift guard.
PHASE2OBJS      := $(BUILDDIR)/phase2/circ.o			\
                   $(BUILDDIR)/phase2/circ_cache.o		\
                   $(BUILDDIR)/phase2/paulis.o			\
                   $(BUILDDIR)/phase2/prob.o			\
                   $(BUILDDIR)/phase2/qreg.o			\
                   $(BUILDDIR)/phase2/state_prep_coeff.o	\
                   $(BUILDDIR)/phase2/world.o			\
                   $(BACKEND_OBJS)
CIRCOBJS        := $(BUILDDIR)/circ/cmpsit.o			\
                   $(BUILDDIR)/circ/qdrift.o			\
                   $(BUILDDIR)/circ/trott.o			\
                   $(BUILDDIR)/circ/trott2.o
LIBOBJS         := $(BUILDDIR)/lib/combinations.o		\
                   $(BUILDDIR)/lib/det_small.o			\
                   $(BUILDDIR)/lib/log.o			\
                   $(BUILDDIR)/lib/xoshiro256ss.o
PH2RUN_DATA_OBJS := $(BUILDDIR)/ph2run/data.o

# --- CFLAGS / LDFLAGS / LDLIBS pipeline ------------------------------------ #
CFLAGS          += -I$(INCLUDE) $(MPI_CFLAGS) $(HDF5_CFLAGS)	\
                   $(BACKEND_CFLAGS)				\
                   -DPHASE2_VERSION=\"$(VERSION)\"		\
                   $(EXTRA_CFLAGS)
LDFLAGS         += $(MPI_LDFLAGS) $(HDF5_LDFLAGS)		\
                   $(BACKEND_LDFLAGS) $(EXTRA_LDFLAGS)
LDLIBS          += -lm $(MPI_LDLIBS) $(HDF5_LDLIBS) $(BACKEND_LDLIBS)

# --- Helpers --------------------------------------------------------------- #
RM              := rm -fv
MKDIR           := mkdir -p

# --- Export contract -- visible to every sub-make -------------------------- #
export CC CFLAGS LDFLAGS LDLIBS NVCC NVCCFLAGS
export TOPDIR INCLUDE BUILDDIR
export PHASE2DIR CIRCDIR LIBDIR PH2RUNDIR BENCHDIR TESTDIR
export PHASE2OBJS CIRCOBJS LIBOBJS PH2RUN_DATA_OBJS BACKEND_OBJS BACKEND_CFLAGS
export BACKEND BACKEND_N
export MPI_CFLAGS MPI_LDFLAGS MPI_LDLIBS MPIRUN MPIFLAGS MPIRANKS
export HDF5_CFLAGS HDF5_LDFLAGS HDF5_LDLIBS
export RM MKDIR


# ========================================================================== #
# Section 3: TARGETS                                                         #
# ========================================================================== #

.DEFAULT_GOAL := all

.PHONY: all build debug shared						\
	phase2 circ lib ph2run bench test-build				\
	build-bench build-test build-test-slow				\
	check check-mpi check-slow check-% check-srcs-coverage		\
	check-tests-coverage						\
	bench bench-build bench-mpi bench-clean				\
	clean distclean format						\
	test-asan test-valgrind test-mpi-asan

all: build build-bench build-test

debug: CFLAGS  += -DDEBUG -g -Og
debug: all

# --- Per-subsystem dispatch ------------------------------------------------ #
# Each phony target invokes the matching subsys Makefile.  Dep edges
# between them (`ph2run: phase2 circ lib`) serialise across subsystem
# boundaries; -j parallelises inside each subsys and between independent
# subsystems (phase2, circ, lib).
phase2:
	+@$(MAKE) -C $(PHASE2DIR)

circ:
	+@$(MAKE) -C $(CIRCDIR)

lib:
	+@$(MAKE) -C $(LIBDIR)

ph2run: phase2 circ lib
	+@$(MAKE) -C $(PH2RUNDIR)

bench-build: phase2 circ lib ph2run-data
	+@$(MAKE) -C $(BENCHDIR)

test-build: phase2 circ lib ph2run-data
	+@$(MAKE) -C $(TESTDIR) build

# Internal: build ph2run/data.o (used by tests and benches) without
# linking the ph2run binary.
.PHONY: ph2run-data
ph2run-data: phase2 circ lib
	+@$(MAKE) -C $(PH2RUNDIR) data

# --- High-level aliases ---------------------------------------------------- #
build:       check-srcs-coverage ph2run
build-bench: bench-build
build-test:  test-build

build-test-slow: phase2 circ lib ph2run-data
	+@$(MAKE) -C $(TESTDIR) build-slow

# --- Shared library (libphase2.so) ----------------------------------------- #
# Pure compute surface, no HDF5; for Python callers via ctypes.  -fPIC
# objects build under $(BUILDDIR)/shared/<srcdir>/ on a disjoint path
# from the regular tree, so `make build` and `make shared` cannot mix
# PIC and non-PIC objects.
SHARED_LDFLAGS := $(MPI_LDFLAGS) $(BACKEND_LDFLAGS)
SHARED_LDLIBS  := $(MPI_LDLIBS) $(BACKEND_LDLIBS)
SHARED_OBJS    := $(BUILDDIR)/shared/phase2/phase2_run.o		\
                  $(patsubst $(BUILDDIR)/%,$(BUILDDIR)/shared/%,	\
                             $(PHASE2OBJS) $(LIBOBJS))

shared: $(SHARED_OBJS)
	@$(MKDIR) $(BUILDDIR)
	$(CC) -shared -o $(BUILDDIR)/libphase2.so $^		\
	      $(SHARED_LDFLAGS) $(SHARED_LDLIBS)

# Dispatch the -fPIC compiles into the relevant subsys's `shared` target.
$(SHARED_OBJS): | phase2-shared lib-shared

.PHONY: phase2-shared lib-shared
phase2-shared:
	+@$(MAKE) -C $(PHASE2DIR) shared
lib-shared:
	+@$(MAKE) -C $(LIBDIR) shared

# --- Tests ----------------------------------------------------------------- #
check: test-build
	+@$(MAKE) -C $(TESTDIR) check

check-mpi: test-build
	+@$(MAKE) -C $(TESTDIR) check-mpi

check-slow: build-test-slow
	+@$(MAKE) -C $(TESTDIR) check-slow

# Pattern: `make check-data` -> runner --filter='data*'.
# Explicit targets above win over this pattern.
check-%: test-build
	+@$(MAKE) -C $(TESTDIR) check-$*

# --- Benchmarks ------------------------------------------------------------ #
# `bench` is the user-facing target: build then run.  `bench-build` is the
# internal dispatch that only builds.  `bench-mpi` runs under mpirun.
bench: bench-build
	@for bb in $(BUILDDIR)/bench/b-paulis $(BUILDDIR)/bench/b-qreg; do \
		"$$bb" || { echo "$$bb: FAIL"; exit 1; };		\
	done

bench-mpi: bench-build
	@for bb in $(BUILDDIR)/bench/b-paulis $(BUILDDIR)/bench/b-qreg; do \
		$(MPIRUN) -n $(MPIRANKS) $(MPIFLAGS) "$$bb" ||		\
			{ echo "$$bb: FAIL"; exit 1; };			\
	done

# Reset the host's bench baseline.  Useful after an intentional
# perf-affecting change has landed: the next `make bench` will show
# `--` in the delta column instead of pretending the new numbers
# regressed against the old ones.
bench-clean:
	@$(RM) bench/runs/$$(hostname).jsonl

# --- Drift guards ---------------------------------------------------------- #
# Each subsys checks its own *.c/*.cu against its own SRCS list.
check-srcs-coverage:
	+@$(MAKE) -C $(PHASE2DIR) check-srcs-coverage
	+@$(MAKE) -C $(CIRCDIR)   check-srcs-coverage
	+@$(MAKE) -C $(LIBDIR)    check-srcs-coverage
	+@$(MAKE) -C $(PH2RUNDIR) check-srcs-coverage
	+@$(MAKE) -C $(BENCHDIR)  check-srcs-coverage

check-tests-coverage:
	+@$(MAKE) -C $(TESTDIR) check-tests-coverage

# --- Sanitiser / valgrind -------------------------------------------------- #
ASAN_FLAGS       := -fsanitize=address,undefined -fno-omit-frame-pointer -g
ASAN_OPTIONS_VAL := detect_leaks=0:halt_on_error=1
export ASAN_OPTIONS_VAL

test-asan:
	$(MAKE) distclean
	ASAN_OPTIONS=$(ASAN_OPTIONS_VAL)				\
	$(MAKE) check							\
	      EXTRA_CFLAGS="$(ASAN_FLAGS)"				\
	      EXTRA_LDFLAGS="$(ASAN_FLAGS)"

test-valgrind: test-build
	+@$(MAKE) -C $(TESTDIR) test-valgrind

test-mpi-asan:
	$(MAKE) distclean
	$(MAKE) build-test						\
	      EXTRA_CFLAGS="$(ASAN_FLAGS)"				\
	      EXTRA_LDFLAGS="$(ASAN_FLAGS)"
	+@$(MAKE) -C $(TESTDIR) test-mpi-asan

# --- Clean ----------------------------------------------------------------- #
# Every compile artefact (objects, .d, binaries, libphase2.so) lives under
# $(BUILDDIR), so `clean` is a single recursive remove.  `distclean` is an
# alias kept for habit.
clean:
	@$(RM) -r $(BUILDDIR)

distclean: clean

# --- Misc ------------------------------------------------------------------ #
format:
	@find $(TOPDIR) -name '*.c' -o -name '*.h' -o			\
	                -name '*.cpp' -o -name '*.cu' |			\
		while read f; do clang-format --style=file -i "$$f"; done


# ========================================================================== #
# Auto-dep includes.  -MMD -MP in CFLAGS makes gcc emit <obj>.d alongside    #
# every <obj>.o; the .d names the .o target and lists every header the .c   #
# #include'd transitively.  Wildcard is silently empty on first build.      #
# ========================================================================== #

-include $(shell find $(BUILDDIR) -name '*.d' 2>/dev/null)
-include $(wildcard $(TESTDIR)/*.d)
