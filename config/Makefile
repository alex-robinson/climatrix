.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir   = src
objdir   = libclimatrix/include
bindir   = libclimatrix/bin
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
## List of rules and source files
##
###############################################

include config/Makefile_files.mk

# EXTRA RULES ####
$(objdir)/nml.o : $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
	
$(objdir)/ncio.o : $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/interp1D.o : $(libdir)/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

extra_libs = 	$(objdir)/nml.o \
				$(objdir)/ncio.o \
				$(objdir)/interp1D.o

###############################################
##
## Compilation of complete programs
##
###############################################

CLIMATRIXROOT = ${CURDIR}
INC_CLIMATRIX = -I${CLIMATRIXROOT}/libclimatrix/include 
LIB_CLIMATRIX = -L${CLIMATRIXROOT}/libclimatrix/include -lclimatrix

# MAIN PROGRAMS ####

# Static library compilation
climatrix-static: $(extra_libs) $(climatrix_base)
	ar rc $(objdir)/libclimatrix.a $(extra_libs) $(climatrix_base)
	ranlib $(objdir)/libclimatrix.a
	@echo " "
	@echo "    $(objdir)/libclimatrix.a is ready."
	@echo " "

test: climatrix-static
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_CLIMATRIX) -o libclimatrix/bin/test_climatrix.x \
		tests/test_climatrix.f90 $(LIB_CLIMATRIX) $(LFLAGS)
	@echo " "
	@echo "    libclimatrix/bin/test_climatrix.x is ready."
	@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make test  : compiles test program"
	@echo " make clean : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
