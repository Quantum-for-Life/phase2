# ========================================================================== #
# lib -- utilities (combinations, det_small, log, prob, xoshiro256ss)        #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := lib
include $(TOPDIR)/build-rules.mk

SRCS  := combinations.c det_small.c log.c prob.c xoshiro256ss.c
OBJS  := $(SRCS:%.c=$(BUILDDIR)/lib/%.o)
SOBJS := $(SRCS:%.c=$(BUILDDIR)/shared/lib/%.o)

.PHONY: all shared check-srcs-coverage clean
all:    $(OBJS)
shared: $(SOBJS)

check-srcs-coverage:
	@declared=" $(SRCS) ";						\
	rc=0;								\
	for s in *.c; do						\
		[ -f "$$s" ] || continue;				\
		case "$$declared" in *" $$s "*) ;;			\
		*) echo "lib: $$s missing from SRCS"; rc=1 ;;		\
		esac;							\
	done;								\
	exit $$rc

clean:
	@$(RM) -r $(BUILDDIR)/lib $(BUILDDIR)/shared/lib
