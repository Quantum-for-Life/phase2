# ========================================================================== #
# test -- C test suite, Python harness, cargo-style parallel runner          #
#                                                                            #
# Builds every t-*.c into a binary at $(BUILDDIR)/test/t-* and the           #
# standalone runner at $(BUILDDIR)/test/run.  Owns the check / check-mpi /   #
# check-slow / check-<filter> targets that delegate to the runner.  Two     #
# drift guards live here: check-tests-coverage (every t-*.c is declared in  #
# TESTS) and the per-subsys check-srcs-coverage (handled by sibling         #
# Makefiles).                                                                #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := test
include $(TOPDIR)/build-rules.mk

# Default-suite tests.  Names listed bare; the
# $(BUILDDIR)/test/ prefix is added below.
TESTS_BASE := \
	t-bench t-bitstring_index t-circ_cache			\
	t-circ_prepst_coeff					\
	t-circ_trott t-circ_trott2 t-circ_trott2_coeff		\
	t-circ_trott_coeff t-circ t-combinations		\
	t-data_attr t-data_coeff_matrix t-data_hamil		\
	t-data_hamil_validate t-data_multidet t-data_open	\
	t-data_dets_validate t-data_idempotence t-data_mpi	\
	t-data_trott_steps t-det_small t-log t-log_release	\
	t-paulis t-prob t-qreg t-ref-bendazzoli t-run		\
	t-state_prep_coeff_csf t-state_prep_coeff_expand	\
	t-world

TBIN       := $(BUILDDIR)/test
TESTS      := $(TESTS_BASE:%=$(TBIN)/%)
TESTS_SLOW := $(TBIN)/t-state_prep_coeff_large
PYTEST     := $(TESTDIR)/t-ref-coeff_matrix.py

# Cargo-style parallel runner + synthetic fixtures
# driven by t-run.
RUNNER     := $(TBIN)/run
RUNFIXDIR  := $(TBIN)/run-fixtures
RUNFIX     := $(addprefix $(RUNFIXDIR)/, pass fail sleep abort banner)

# ----- Targets -------------------------------------------------------------
.PHONY: build build-slow check check-mpi check-slow check-% \
	check-tests-coverage clean				\
	test-valgrind test-mpi-asan

build: check-tests-coverage $(TESTS) $(RUNNER) $(RUNFIX)

build-slow: $(TESTS_SLOW)

# Test binary build.  Each $(BUILDDIR)/test/t-X builds
# from $(TESTDIR)/t-X.c plus the upstream subsys
# objects.  Inter-subsys link prereqs (PHASE2OBJS,
# CIRCOBJS, LIBOBJS, PH2RUN_DATA_OBJS) come from the
# top-level export.  -DDEBUG is added via target-
# specific CFLAGS so the more-specific -UDEBUG on
# t-log_release lands later on the gcc command line
# and wins (the C preprocessor's last-wins rule).
$(TESTS) $(TESTS_SLOW): CFLAGS += -I$(TESTDIR) -I$(PHASE2DIR)		\
                                  -DPH2_TESTDIR=\"./test\"		\
                                  -DDEBUG -g -Og

$(TESTS) $(TESTS_SLOW): $(TBIN)/%: $(TESTDIR)/%.c			\
                        $(TESTDIR)/test.h $(TESTDIR)/t-data.h		\
                        $(PHASE2OBJS) $(CIRCOBJS) $(LIBOBJS)		\
                        $(PH2RUN_DATA_OBJS)
	@$(MKDIR) $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) $< $(filter %.o, $^) $(LDLIBS) -o $@

# Tests that need a specific algorithm object linked
# in directly (the test file references the algorithm
# symbols without going through phase2_run dispatch).
$(TBIN)/t-circ_cache:        $(BUILDDIR)/circ/trott.o
$(TBIN)/t-circ_trott:        $(BUILDDIR)/circ/trott.o
$(TBIN)/t-circ_trott2:       $(BUILDDIR)/circ/trott2.o
$(TBIN)/t-circ_trott_coeff:  $(BUILDDIR)/circ/trott.o
$(TBIN)/t-circ_trott2_coeff: $(BUILDDIR)/circ/trott2.o

# t-bench tests bench/bench.h directly; the header
# lives outside the default test include path.
$(TBIN)/t-bench: CFLAGS += -I$(BENCHDIR)

# t-log_release verifies the release-build strip of
# trace/debug macros: -UDEBUG cancels the -DDEBUG
# inherited from the pattern target-specific above.
$(TBIN)/t-log_release: CFLAGS += -UDEBUG

# Cargo-style parallel runner: standalone C program,
# libc + POSIX only.  No phase2 / MPI / HDF5 deps.
$(RUNNER): $(TESTDIR)/run.c
	@$(MKDIR) $(@D)
	$(CC) -std=c11 -Wall -Wextra -O2 -o $@ $<

# Synthetic fixtures used by t-run to exercise the
# runner against controlled outcomes.
$(RUNFIXDIR)/%: $(TESTDIR)/run-fixtures/%.c
	@$(MKDIR) $(@D)
	$(CC) -std=c11 -Wall -Wextra -O2 -o $@ $<

# t-run is a meta-test: it shells out to the runner
# against the fixtures.  Order-only prereqs ensure
# they exist at runtime but stay out of the implicit
# %: %.c link line.
$(TBIN)/t-run: | $(RUNNER) $(RUNFIX)

# ----- Drift guard: every t-*.c is in TESTS or TESTS_SLOW -------------------
check-tests-coverage:
	@tmp=$$(mktemp -d);					\
	ls $(TESTDIR)/t-*.c 2>/dev/null				\
		| sed 's|^$(TESTDIR)/||; s|\.c$$||'		\
		| sort > $$tmp/expected;			\
	for t in $(TESTS) $(TESTS_SLOW); do			\
		echo $${t##*/};					\
	done | sort > $$tmp/declared;				\
	missing=$$(comm -23 $$tmp/expected $$tmp/declared);	\
	rc=0;							\
	if [ -n "$$missing" ]; then				\
		echo "test/: t-*.c files not in TESTS:";	\
		echo "$$missing";				\
		rc=1;						\
	fi;							\
	rm -rf $$tmp;						\
	exit $$rc

# ----- check targets (delegate to the runner) -------------------------------
check: build
	@cd $(TOPDIR) && $(RUNNER) -- $(TESTS) $(PYTEST)

check-mpi: build
	@cd $(TOPDIR) && $(RUNNER) --mpiranks=$(MPIRANKS) -- $(TESTS)

check-slow: build-slow $(RUNNER)
	@cd $(TOPDIR) && $(RUNNER) -- $(TESTS_SLOW)

# Pattern: make check-data, check-paulis, check-circ,
# etc. -- threads the stem as <stem>* fnmatch into
# the runner.
check-%: build
	@cd $(TOPDIR) && $(RUNNER) --filter='$**' -- $(TESTS) $(PYTEST)

# ----- Valgrind and MPI-ASan (do not go through the runner) -----------------
test-valgrind: build
	@cd $(TOPDIR) && for tt in $(TESTS); do				\
		valgrind --quiet --error-exitcode=1			\
			--leak-check=full				\
			--errors-for-leak-kinds=all			\
			--track-origins=yes $$tt &&			\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

test-mpi-asan: build
	@cd $(TOPDIR) && for tt in $(TBIN)/t-circ_trott_coeff; do	\
		$(MPIRUN) -n 4 $(MPIFLAGS)				\
			-x ASAN_OPTIONS=$(ASAN_OPTIONS_VAL) $$tt &&	\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

# ----- Clean ----------------------------------------------------------------
# Test binaries + .d files all live under $(BUILDDIR)/test/; top-level
# `make clean` wipes $(BUILDDIR) wholesale.  Nothing to clean here.
clean:
	@:
