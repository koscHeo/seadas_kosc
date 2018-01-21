cmake_minimum_required (VERSION 2.6)

project (HEALPIX)

ADD_LIBRARY(healpix STATIC
    healpix_types.F90
#    fits_template.f90
    cgetEnvironment.c
    indmed.f90
#    indmed_part1.f90
#    indmed_part2.f90
    utilities.f90
    bit_manipulation.F90
    extension.F90
    long_intrinsic.F90
    obsolete.f90
    ran_tools_dist.f90
    rngmod.f90
    healpix_fft.F90
    misc_utils.F90
    head_fits.F90
    num_rec.F90
    statistics.f90
    paramfile_io.F90
    pix_tools.F90
    coord_v_convert.f90
    fitstools.F90
    udgrade_nr.f90
    mask_tools.F90
    alm_tools.F90
    healpix_modules.f90
  )
# objsharp = ../sharp/libsharp_all.o

target_link_libraries(healpix libcfitsio)

install (TARGETS healpix ARCHIVE DESTINATION $ENV{LIB3_DIR}/lib)

