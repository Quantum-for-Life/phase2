# ========================================================================== #
# phase2 -- state-vector kernel, Pauli arithmetic, register layout, world    #
#                                                                            #
# Compiles phase2/*.c into $(BUILDDIR)/phase2/*.o.  Backend-conditional      #
# sources (qreg_qreg.c for CPU, qreg_cuda*.c + world_cuda.c for GPU) are     #
# pulled in via $(BACKEND_OBJS) exported by the top-level Makefile.  The     #
# CUDA NVCC rule for qreg_cuda_lo.cu lives here too.                         #
# ========================================================================== #

ifndef TOPDIR
$(error TOPDIR not set -- invoke via top-level Makefile, e.g. `make -C ..`)
endif

SUBSYS := phase2
include $(TOPDIR)/build-rules.mk

# Always-compiled sources.
SRCS  := circ.c circ_cache.c paulis.c prob.c qreg.c \
         state_prep_coeff.c world.c

# Extra source compiled only for the shared library
# (libphase2.so), not for ph2run.
SRCS_SHARED_EXTRA := phase2_run.c

# Backend-conditional sources (only some are built per
# the active BACKEND, but all are *declared* so the
# drift guard recognises them as known files).
SRCS_BACKEND := qreg_qreg.c qreg_cuda.c \
                qreg_cuda_lo.cu world_cuda.c

OBJS          := $(SRCS:%.c=$(BUILDDIR)/phase2/%.o)
SOBJS         := $(SRCS:%.c=$(BUILDDIR)/shared/phase2/%.o)	\
                 $(BUILDDIR)/shared/phase2/phase2_run.o
BACKEND_SOBJS := $(patsubst $(BUILDDIR)/%,$(BUILDDIR)/shared/%, $(BACKEND_OBJS))

.PHONY: all shared check-srcs-coverage clean
all:    $(OBJS) $(BACKEND_OBJS)
shared: $(SOBJS) $(BACKEND_SOBJS)

# NVCC compile rule for the one .cu source.  Only fires
# when BACKEND=cuda places qreg_cuda_lo.o in
# $(BACKEND_OBJS).
$(BUILDDIR)/phase2/qreg_cuda_lo.o: qreg_cuda_lo.cu
	@$(MKDIR) $(@D)
	$(NVCC) $(NVCCFLAGS) $(BACKEND_CFLAGS) -I$(INCLUDE) -c $< -o $@

$(BUILDDIR)/phase2/qreg_cuda_lo_dlink.o: $(BUILDDIR)/phase2/qreg_cuda_lo.o
	$(NVCC) $(NVCCFLAGS) -dlink $< -o $@

# Drift guard: every .c/.cu in this dir must appear in
# SRCS, SRCS_SHARED_EXTRA, or SRCS_BACKEND.  Fails the
# build if a new source slips in without being declared.
check-srcs-coverage:
	@declared=" $(SRCS) $(SRCS_SHARED_EXTRA) $(SRCS_BACKEND) ";	\
	rc=0;								\
	for s in *.c *.cu; do						\
		[ -f "$$s" ] || continue;				\
		case "$$declared" in *" $$s "*) ;;			\
		*) echo "phase2: $$s missing from SRCS"; rc=1 ;;	\
		esac;							\
	done;								\
	exit $$rc

clean:
	@$(RM) -r $(BUILDDIR)/phase2 $(BUILDDIR)/shared/phase2
