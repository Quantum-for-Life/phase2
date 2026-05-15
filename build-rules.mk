# ========================================================================== #
# build-rules.mk -- shared compile rules for the phase2 build                #
#                                                                            #
# Included by every per-subsystem Makefile after that Makefile sets SUBSYS   #
# to its own directory name.  Two pattern rules emit objects under           #
# $(BUILDDIR)/<subsys>/ (regular) and $(BUILDDIR)/shared/<subsys>/ (-fPIC,   #
# for the shared library).  Header deps are auto-tracked via gcc -MMD -MP    #
# (set in CFLAGS at the top level); the .d files are pulled in at the foot  #
# of the top-level Makefile.                                                 #
# ========================================================================== #

$(BUILDDIR)/$(SUBSYS)/%.o: %.c
	@$(MKDIR) $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/shared/$(SUBSYS)/%.o: %.c
	@$(MKDIR) $(@D)
	$(CC) $(CFLAGS) -fPIC -c $< -o $@
