#  This makefile builds the "atrem" program from source files in the
#    current directory.  All objects files on which the final executable
#    depends are listed in the OBJECTS macro.  Any new files needed to
#    build the executable should have their object-file form listed there
#    as well.  To actually build the executable version of "atrem",
#     simply enter "make" at the shell prompt while in this directory.  All
#    outdated object files will be recompiled and the executable relinked
#    if this is necessary.  To recompile everything, first type "make clean"
#    at the shell prompt and this will remove all ".o" files and the "atrem"
#    executable.  Then type "make" at the shell prompt.
#
# Define these for your system
#F90 = /home/mmontes/bin/g95
F90 = /sw/bin/gfortran
COPTS = -c -O2  
#COPTS = -c -O2 -OPT:Olimit=0 -static
#COPTS = -c -O2 -OPT:Olimit=0
LOPTS = -O2 -static-libgfortran
#LOPTS = -O2 -static
CP = /bin/cp
CHMOD = /bin/chmod
MV = /bin/mv
RM = /bin/rm
LN = /bin/ln
MKDIR = /bin/mkdir
TAR = /usr/bin/tar
GZIP = /usr/bin/gzip -v --best
# result of DATE should be something like 20011120
DATE = /bin/date +%Y%m%d
# result of FDATE should be something like 2001Nov20 
FDATE = /bin/date +%Y%b%d
#  Define a macro containing the names of all needed object files.
# Is there someway to not have to define F_OBJECTS and OBJECTS?
# I only use F_OBJECTS in directory maintainence, tarfiles, etc.
F_OBJECTS= SUNCOR.f90 SOLCOR.f90 \
JULIAN.f90 HAZEL.f90 file_io_new.f90 write_envi_header.f90 \
envi_header.f90 check_range.f90 errors.f90 auxtools.f90 chartools.f90 \
default_character_length.f90 generic_info.f90 ASTMATH.FOR \
ecef2latlon.f90 solar_geometry.f90 ASTREDUC.FOR ASTTIME.f q_to_r.f90 params.f90 \
proc_hico.f90 get_astro_data.f90 intercept.f90 bore_site.f90
F_CMN= ASTMATH.CMN ASTREDUC.CMN
OBJECTS= SUNCOR.o SOLCOR.o JULIAN.o \
HAZEL.o file_io_new.o write_envi_header.o \
envi_header.o check_range.o errors.o auxtools.o ASTMATH.o \
chartools.o default_character_length.o generic_info.o \
ecef2latlon.o solar_geometry.o ASTREDUC.o ASTTIME.o q_to_r.o params.o \
proc_hico.o get_astro_data.o intercept.o bore_site.o
MODULES= bore_site.mod intercept.mod quat_to_rot.mod auxtools.mod chartools.mod \
check_range.mod default_character_length.mod errors_mod.mod hazel_mod.mod julian_mod.mod \
solcor_mod.mod suncod_mod.mod sol_geo_mod get_ast_data.mod ecef2llh.mod generic_info.mod \
file_io_mod.mod envi_header.mod envi_header_mod.mod suncor_mod.mod check_range_mod.mod \
sol_geo.mod params.mod pos_vel_quat.mod
#  Define a macro containing the name of the executable file.
EXEFIL= proc_hico_mjm
#$(EXEFIL):	$(OBJECTS) $(F_OBJECTS)
$(EXEFIL):	$(OBJECTS)
	$(F90) -o $(EXEFIL)_tmp $(LOPTS) $(OBJECTS)
	$(CHMOD) 774 $(EXEFIL)_tmp
	$(MV) $(EXEFIL)_tmp $(EXEFIL)_`$(FDATE)`
	$(RM) -f $(EXEFIL)
	$(LN) -s $(EXEFIL)_`$(FDATE)` $(EXEFIL)
#	mv $(EXEFIL)_tmp $(EXEFIL)_`/sbin/date +%Y%m%d`
#dependencies (takes care of changes to compile options this way):
$(OBJECTS): makefile
# Normal file-to-file dependencies
# generated with g95 -M *.f90
# My generic statements were not working, so explicitly specify them....
bore_site.o bore_site.mod: bore_site.f90
	$(F90) $(COPTS) bore_site.f90
intercept.o intercept.mod: intercept.f90 ecef2llh.mod
	$(F90) $(COPTS) intercept.f90
params.o params.mod: params.f90
	$(F90) $(COPTS) params.f90
q_to_r.o quat_to_rot.mod: q_to_r.f90
	$(F90) $(COPTS) q_to_r.f90
ASTREDUC.o: ASTREDUC.FOR ASTREDUC.CMN
	$(F90) $(COPTS) ASTREDUC.FOR
ASTTIME.o: ASTTIME.f
	$(F90) $(COPTS) ASTTIME.f
ASTMATH.o: ASTMATH.FOR ASTMATH.CMN
	$(F90) $(COPTS) ASTMATH.FOR
auxtools.o auxtools.mod: auxtools.f90
	$(F90) $(COPTS) auxtools.f90
chartools.o chartools.mod: chartools.f90
	$(F90) $(COPTS) chartools.f90
check_range.o check_range_mod.mod: check_range.f90
	$(F90) $(COPTS) check_range.f90
default_character_length.o default_character_length.mod: default_character_length.f90
	$(F90) $(COPTS) default_character_length.f90
errors.o errors_mod.mod: errors.f90
	$(F90) $(COPTS) errors.f90
HAZEL.o hazel_mod.mod: HAZEL.f90
	$(F90) $(COPTS) HAZEL.f90
JULIAN.o julian_mod.mod: JULIAN.f90
	$(F90) $(COPTS) JULIAN.f90
SOLCOR.o solcor_mod.mod: SOLCOR.f90
	$(F90) $(COPTS) SOLCOR.f90
SUNCOR.o suncor_mod.mod: julian_mod.mod solcor_mod.mod SUNCOR.f90
	$(F90) $(COPTS) SUNCOR.f90
solar_geometry.o sol_geo.mod: suncor_mod.mod hazel_mod.mod solar_geometry.f90
	$(F90) $(COPTS) solar_geometry.f90
get_astro_data.o get_ast_data.mod: get_astro_data.f90
	$(F90) $(COPTS) get_astro_data.f90
ecef2latlon.o ecef2llh.mod: ecef2latlon.f90 ASTMATH.CMN
	$(F90) $(COPTS) ecef2latlon.f90
generic_info.o generic_info.mod: default_character_length.mod generic_info.f90
	$(F90) $(COPTS) generic_info.f90
file_io_new.o file_io_mod.mod: errors_mod.mod file_io_new.f90
	$(F90) $(COPTS) file_io_new.f90
envi_header.o envi_header.mod: errors_mod.mod default_character_length.mod envi_header.f90 chartools.mod
	$(F90) $(COPTS) envi_header.f90
write_envi_header.o envi_header_mod.mod: errors_mod.mod file_io_mod.mod default_character_length.mod check_range_mod.mod envi_header.mod write_envi_header.f90
	$(F90) $(COPTS) write_envi_header.f90
proc_hico.o: proc_hico.f90 quat_to_rot.mod sol_geo.mod ecef2llh.mod envi_header_mod.mod auxtools.mod params.mod quat_to_rot.mod ASTTIME.o ASTMATH.o ASTREDUC.o params.mod ASTMATH.CMN generic_info.mod get_ast_data.mod intercept.mod bore_site.mod
	$(F90) $(COPTS) proc_hico.f90



# Include .f90 as a suffix, especially if  make -fp - </dev/null  (on
#an SGI) does not have  .f90 in the list of suffixes. On the SGI, the
#following *ADDS* it to the list without  clobbering the pre-existing
#list. Experience shows that this line needs to occur before  .f90.o
#rule that is below it. Thus, a "file.o" means we can find a "file.f",
#"file.f90",  "file.c" etc, in order to apply various builtin rules, as
#well as those below.
# For gnu make, for whatever reason I need to clear the suffix list first. 
# Some of the files were being seen as *.mod.0 meaning m2c -o *.o *.mod. 
# I cannot figure out why only some and not all. I don't want it to behave 
# that way, but I must be missing something. So, clearout all previous rules,
# and use the only rules that we need. 
.SUFFIXES: 
.SUFFIXES: .f90 .o .mod
# Some macros including the arguements that we want. 
.f90.o:
	$(F90)  $(COPTS) $<
#	$(F90)  -static -c -O2 $<
.f90.mod:
	$(F90)  $(COPTS) $<
#.mod:
#	$(F90)  $(COPTS) $<
#	$(F90)  -static -c -O2 $<
clean:
	$(RM) -f $(OBJECTS) $(MODULES) $(EXEFIL)
frozen_dir:
	if test ! -d frozen_$(EXEFIL)_`$(FDATE)`; then $(MKDIR) frozen_$(EXEFIL)_`$(FDATE)`; fi
	$(CP) -p -f $(F_OBJECTS) $(F_CMN) frozen_$(EXEFIL)_`$(FDATE)`
	$(CP) -p -f makefile frozen_$(EXEFIL)_`$(FDATE)`
tarfile:
	$(TAR) cvf frozen_$(EXEFIL)_`$(FDATE)`.tar $(F_OBJECTS) $(F_CMN) makefile
gztarfile:
	$(TAR) cvf frozen_$(EXEFIL)_`$(FDATE)`.tar $(F_OBJECTS) $(F_CMN) makefile
	$(GZIP) frozen_$(EXEFIL)_`$(FDATE)`.tar
#
