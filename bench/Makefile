# ========================================================================== #
# bench -- micro-benchmarks (b-paulis, b-qreg)                               #
#                                                                            #
# Compiles bench/bench.c into $(BUILDDIR)/bench/bench.o (a shared support    #
# unit used by every b-*) and links the benchmark binaries from b-*.c plus  #
# the prebuilt subsys objects.                                              #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := bench
include $(TOPDIR)/build-rules.mk

# Per-benchmark source files; each builds to a binary
# at $(BENCHDIR)/<name>.
PROGS := b-paulis b-qreg

# Shared support code linked into every benchmark.
SHARED_SRCS := bench.c
SHARED_OBJS := $(SHARED_SRCS:%.c=$(BUILDDIR)/bench/%.o)

# bench's headers live alongside its sources; tell gcc.
CFLAGS += -I.

.PHONY: all check-srcs-coverage clean
all: $(PROGS)

# Per-binary link rule.  b-NAME built from b-NAME.c +
# bench.o + every upstream subsys object set.
$(PROGS): %: %.c bench.h $(SHARED_OBJS) \
              $(PHASE2OBJS) $(LIBOBJS) $(PH2RUN_DATA_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $< $(SHARED_OBJS) \
	      $(PHASE2OBJS) $(LIBOBJS) $(PH2RUN_DATA_OBJS) \
	      $(LDLIBS) -o $@

check-srcs-coverage:
	@declared=" $(SHARED_SRCS) $(PROGS:%=%.c) ";			\
	rc=0;								\
	for s in *.c; do						\
		[ -f "$$s" ] || continue;				\
		case "$$declared" in *" $$s "*) ;;			\
		*) echo "bench: $$s missing from SRCS"; rc=1 ;;		\
		esac;							\
	done;								\
	exit $$rc

clean:
	@$(RM) -r $(BUILDDIR)/bench
	@$(RM) $(PROGS)
