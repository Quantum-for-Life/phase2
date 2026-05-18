# ========================================================================== #
# bench -- micro-benchmarks (b-paulis, b-qreg, b-log, b-circ, b-prob,        #
#                            b-algos)                                        #
#                                                                            #
# Each benchmark binary is built from a single b-NAME.c that includes the    #
# header-only bench.h framework.  No shared bench.o.                         #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := bench
include $(TOPDIR)/build-rules.mk

# Per-benchmark source files; each builds to a binary
# at $(BUILDDIR)/bench/<name>.
PROG_NAMES := b-paulis b-qreg b-log b-circ b-prob b-algos
PROGS      := $(PROG_NAMES:%=$(BUILDDIR)/bench/%)

# bench's headers live alongside its sources; tell gcc.
CFLAGS += -I.

.PHONY: all check-srcs-coverage clean
all: $(PROGS)

# Per-binary link rule.  $(BUILDDIR)/bench/b-NAME built
# from b-NAME.c + every upstream subsys object set
# (CIRCOBJS reaches the circ algorithms for b-algos; the
# others go through PHASE2OBJS / LIBOBJS).
$(PROGS): $(BUILDDIR)/bench/%: %.c bench.h				\
                               $(PHASE2OBJS) $(CIRCOBJS)		\
                               $(LIBOBJS) $(PH2RUN_DATA_OBJS)
	@$(MKDIR) $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) $<					\
	      $(PHASE2OBJS) $(CIRCOBJS) $(LIBOBJS)			\
	      $(PH2RUN_DATA_OBJS) $(LDLIBS) -o $@

check-srcs-coverage:
	@declared=" $(PROG_NAMES:%=%.c) ";				\
	rc=0;								\
	for s in *.c; do						\
		[ -f "$$s" ] || continue;				\
		case "$$declared" in *" $$s "*) ;;			\
		*) echo "bench: $$s missing from SRCS"; rc=1 ;;		\
		esac;							\
	done;								\
	exit $$rc

# Bench binaries live under $(BUILDDIR)/bench/; top-level
# `make clean` wipes $(BUILDDIR) wholesale.
clean:
	@:
