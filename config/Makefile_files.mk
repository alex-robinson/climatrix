
#### CLIMATRIX RULES ####

$(objdir)/climatrix_defs.o : $(srcdir)/climatrix_defs.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $< 

$(objdir)/climatrix.o : $(srcdir)/climatrix.f90 $(objdir)/climatrix_defs.o $(objdir)/nml.o \
						$(objdir)/climate_interpolation.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $< 

$(objdir)/climate_interpolation.o : $(srcdir)/climate_interpolation.f90 $(objdir)/climatrix_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $< 

#### CLIMATRIX LIST ####

climatrix_base =  	$(objdir)/climatrix_defs.o \
					$(objdir)/climatrix.o \
					$(objdir)/climate_interpolation.o



