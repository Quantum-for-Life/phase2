# ========================================================================== #
# circ -- Trotter / Strang-Trotter / qDRIFT / composite circuit algorithms   #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := circ
include $(TOPDIR)/build-rules.mk

SRCS  := cmpsit.c qdrift.c trott.c trott2.c
OBJS  := $(SRCS:%.c=$(BUILDDIR)/circ/%.o)
SOBJS := $(SRCS:%.c=$(BUILDDIR)/shared/circ/%.o)

.PHONY: all shared check-srcs-coverage clean
all:    $(OBJS)
shared: $(SOBJS)

check-srcs-coverage:
	@declared=" $(SRCS) ";						\
	rc=0;								\
	for s in *.c; do						\
		[ -f "$$s" ] || continue;				\
		case "$$declared" in *" $$s "*) ;;			\
		*) echo "circ: $$s missing from SRCS"; rc=1 ;;		\
		esac;							\
	done;								\
	exit $$rc

clean:
	@$(RM) -r $(BUILDDIR)/circ $(BUILDDIR)/shared/circ
