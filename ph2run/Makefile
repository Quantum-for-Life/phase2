# ========================================================================== #
# ph2run -- driver CLI (ph2run binary) and the HDF5 data loader              #
#                                                                            #
# Compiles ph2run/data.c into $(BUILDDIR)/ph2run/data.o and links the        #
# ph2run binary from ph2run.c + the prebuilt phase2 / circ / lib /           #
# ph2run-data objects.  The top-level Makefile ensures the upstream subsys  #
# objects exist before invoking this sub-make.                              #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := ph2run
include $(TOPDIR)/build-rules.mk

# Source of the data loader (compiled to an .o reused
# by tests and benches).
SRCS  := data.c
OBJS  := $(SRCS:%.c=$(BUILDDIR)/ph2run/%.o)
SOBJS :=

# Program source (linked into the ph2run binary).
PROG_SRC := ph2run.c
PROG_BIN := $(BUILDDIR)/ph2run/ph2run

.PHONY: all data check-srcs-coverage clean
all:  $(PROG_BIN)
data: $(OBJS)

# Link ph2run: ph2run.c plus the data.o + every upstream
# subsys object.  Inter-subsys link prereqs come from
# the env (PHASE2OBJS, CIRCOBJS, LIBOBJS) exported by
# the top-level Makefile.
$(PROG_BIN): $(PROG_SRC) $(OBJS) $(PHASE2OBJS) $(CIRCOBJS) $(LIBOBJS)
	@$(MKDIR) $(@D)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

check-srcs-coverage:
	@declared=" $(SRCS) $(PROG_SRC) ";				\
	rc=0;								\
	for s in *.c; do						\
		[ -f "$$s" ] || continue;				\
		case "$$declared" in *" $$s "*) ;;			\
		*) echo "ph2run: $$s missing from SRCS"; rc=1 ;;	\
		esac;							\
	done;								\
	exit $$rc

# ph2run binary lives under $(BUILDDIR)/ph2run/; top-level
# `make clean` wipes $(BUILDDIR) wholesale.
clean:
	@:
