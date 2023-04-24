
#### REMBO RULES ####

$(objdir)/exchange.o : $(srcdir)/exchange.f90 $(objdir)/parameters.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/parameters.o : $(srcdir)/parameters.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/emb_global.o : $(srcdir)/emb_global.f90 $(objdir)/ncio.o $(objdir)/exchange.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $< 
$(objdir)/emb_pdesolver.o : $(srcdir)/emb_pdesolver.f90 $(objdir)/emb_global.o 
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
$(objdir)/projector.o : $(srcdir)/projector.f90 $(objdir)/emb_global.o 
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
$(objdir)/emb_functions.o : $(srcdir)/emb_functions.f90 $(objdir)/emb_global.o $(objdir)/ncio.o 
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
$(objdir)/rembo_functions.o : $(srcdir)/rembo_functions.f90 $(objdir)/emb_global.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
$(objdir)/smb_itm.o : $(srcdir)/smb_itm.f90 $(objdir)/emb_global.o $(objdir)/emb_functions.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
$(objdir)/rembo_main.o : $(srcdir)/rembo_main.f90 $(objdir)/rembo_functions.o $(objdir)/emb_global.o \
$(objdir)/emb_functions.o $(objdir)/emb_pdesolver.o $(objdir)/projector.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<
$(objdir)/climate.o : $(srcdir)/climate.f90 $(objdir)/emb_global.o $(objdir)/emb_functions.o \
 							$(objdir)/rembo_main.o 
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/sinsol_orbit.o : $(srcdir)/sinsol_orbit.f
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<



#############################################################
##                          
## List of rembo files
##
#############################################################

rembo_base =           $(objdir)/exchange.o \
					   $(objdir)/parameters.o \
					   $(objdir)/emb_global.o \
					   $(objdir)/emb_pdesolver.o \
					   $(objdir)/projector.o \
					   $(objdir)/emb_functions.o \
					   $(objdir)/rembo_functions.o \
					   $(objdir)/smb_itm.o \
					   $(objdir)/rembo_main.o \
					   $(objdir)/climate.o \
					   $(objdir)/sinsol_orbit.o



