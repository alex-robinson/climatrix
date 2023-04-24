.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir   = src
objdir   = librembo/include
bindir   = librembo/bin
libdir   = libs

# Command-line options at make call
debug ?= 0

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##							
## List of rembo rules and source files
##
###############################################

include config/Makefile_rembo.mk

# GLOBAL RULES ####
$(objdir)/nml.o : $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
	
$(objdir)/ncio.o : $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/ncio_transpose.o : $(libdir)/ncio_transpose.f90 $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp1D.o : $(libdir)/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

rembo_libs = 	$(objdir)/nml.o \
				$(objdir)/ncio.o \
				$(objdir)/ncio_transpose.o \
				$(objdir)/interp1D.o

###############################################
##
## Compilation of complete programs
##
###############################################

# MAIN PROGRAMS ####

# Static library compilation
rembo-static: $(rembo_libs) $(rembo_base)
	ar rc $(objdir)/librembo.a $(rembo_libs) $(rembo_base)
	ranlib $(objdir)/librembo.a
	@echo " "
	@echo "    $(objdir)/librembo.a is ready."
	@echo " "

rembo: rembo-static
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_REMBO) -o librembo/bin/test_rembo.x \
		test_rembo.f90 $(LIB_REMBO) $(LFLAGS)
	@echo " "
	@echo "    librembo/bin/test_rembo.x is ready."
	@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make rembo : compiles rembo.x"
	@echo " make clean : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
