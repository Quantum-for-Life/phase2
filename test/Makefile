# ========================================================================== #
# test -- C test suite, Python harness, cargo-style parallel runner          #
#                                                                            #
# Builds every t-*.c into a binary at $(TESTDIR)/t-* and the standalone      #
# runner at $(TESTDIR)/run.  Owns the check / check-mpi / check-slow /       #
# check-<filter> targets that delegate to the runner.  Two drift guards     #
# live here: check-tests-coverage (every t-*.c is declared in TESTS) and    #
# the per-subsys check-srcs-coverage (handled by sibling Makefiles).        #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := test
include $(TOPDIR)/build-rules.mk

# Default-suite tests.  Names are listed bare; the
# $(TESTDIR)/ prefix is added below.
TESTS_BASE := \
	t-bitstring_index t-circ_cache t-circ_prepst_coeff	\
	t-circ_trott t-circ_trott2 t-circ_trott2_coeff		\
	t-circ_trott_coeff t-circ t-combinations		\
	t-data_attr t-data_coeff_matrix t-data_hamil		\
	t-data_hamil_validate t-data_multidet t-data_open	\
	t-data_dets_validate t-data_idempotence t-data_mpi	\
	t-data_trott_steps t-det_small t-log t-log_release	\
	t-paulis t-prob t-qreg t-ref-bendazzoli t-run		\
	t-state_prep_coeff_csf t-state_prep_coeff_expand	\
	t-world

TESTS      := $(TESTS_BASE:%=$(TESTDIR)/%)
TESTS_SLOW := $(TESTDIR)/t-state_prep_coeff_large
PYTEST     := $(TESTDIR)/t-ref-coeff_matrix.py

# Cargo-style parallel runner + synthetic fixtures
# driven by t-run.
RUNNER     := $(TESTDIR)/run
RUNFIXDIR  := $(TESTDIR)/run-fixtures
RUNFIX     := $(addprefix $(RUNFIXDIR)/, pass fail sleep abort banner)

# Tests compile with -DDEBUG (so log_trace / log_debug
# in include/log.h expand to real emits).  Applied
# below via target-specific CFLAGS so the more-specific
# `-UDEBUG` on t-log_release lands later on the gcc
# command line and wins (`last wins` rule for the C
# preprocessor).

# ----- Targets -------------------------------------------------------------
.PHONY: build build-slow check check-mpi check-slow check-% \
	check-tests-coverage clean distclean			\
	test-valgrind test-mpi-asan

build: check-tests-coverage $(TESTS) $(RUNNER) $(RUNFIX)

build-slow: $(TESTS_SLOW)

# Test binary build.  Each t-X builds from t-X.c plus
# the upstream subsys objects.  Inter-subsys link
# prereqs (PHASE2OBJS, CIRCOBJS, LIBOBJS,
# PH2RUN_DATA_OBJS) come from the top-level export.
$(TESTS) $(TESTS_SLOW): CFLAGS += -I$(TESTDIR) -I$(PHASE2DIR)		\
                                  -DPH2_TESTDIR=\"./test\"		\
                                  -DDEBUG -g -Og

$(TESTS) $(TESTS_SLOW): $(TESTDIR)/%: $(TESTDIR)/%.c			\
                        $(TESTDIR)/test.h $(TESTDIR)/t-data.h		\
                        $(PHASE2OBJS) $(CIRCOBJS) $(LIBOBJS)		\
                        $(PH2RUN_DATA_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $< $(filter %.o, $^) $(LDLIBS) -o $@

# Tests that need a specific algorithm object linked in
# directly (the test file references the algorithm
# symbols without going through phase2_run dispatch).
$(TESTDIR)/t-circ_cache:        $(BUILDDIR)/circ/trott.o
$(TESTDIR)/t-circ_trott:        $(BUILDDIR)/circ/trott.o
$(TESTDIR)/t-circ_trott2:       $(BUILDDIR)/circ/trott2.o
$(TESTDIR)/t-circ_trott_coeff:  $(BUILDDIR)/circ/trott.o
$(TESTDIR)/t-circ_trott2_coeff: $(BUILDDIR)/circ/trott2.o

# t-log_release cancels -DDEBUG to verify the release
# build's strip of trace/debug macros.  The override
# fires after TEST_CFLAGS via target-specific CFLAGS.
$(TESTDIR)/t-log_release: CFLAGS += -UDEBUG

# Cargo-style parallel runner: standalone C program,
# libc + POSIX only.  No phase2 / MPI / HDF5 deps.
$(RUNNER): $(TESTDIR)/run.c
	$(CC) -std=c11 -Wall -Wextra -O2 -o $@ $<

# Synthetic fixtures used by t-run to exercise the
# runner against controlled outcomes.
$(RUNFIXDIR)/%: $(RUNFIXDIR)/%.c
	$(CC) -std=c11 -Wall -Wextra -O2 -o $@ $<

# t-run is a meta-test: it shells out to ./test/run
# against the fixtures.  Order-only prereqs ensure
# they exist at runtime but stay out of the implicit
# %: %.c link line (which would try to feed binaries
# to ld).
$(TESTDIR)/t-run: | $(RUNNER) $(RUNFIX)

# ----- Drift guard: every t-*.c is in TESTS or TESTS_SLOW -------------------
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
	@for tt in $(TESTS); do						\
		valgrind --quiet --error-exitcode=1			\
			--leak-check=full				\
			--errors-for-leak-kinds=all			\
			--track-origins=yes ./$$tt &&			\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

test-mpi-asan: build
	@for tt in $(TESTDIR)/t-circ_trott_coeff; do			\
		$(MPIRUN) -n 4 $(MPIFLAGS)				\
			-x ASAN_OPTIONS=$(ASAN_OPTIONS_VAL) ./$$tt &&	\
			echo "$$tt: OK" ||				\
			( echo "$$tt: FAIL"; exit 1 );			\
	done

# ----- Clean ----------------------------------------------------------------
clean:
	@$(RM) -r $(BUILDDIR)/test
	@$(RM) $(TESTDIR)/*.d

distclean: clean
	@$(RM) $(TESTS) $(TESTS_SLOW)
	@$(RM) $(RUNNER) $(RUNFIX)
