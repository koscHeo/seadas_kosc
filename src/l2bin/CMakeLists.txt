cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l2bin
  l2bin.c
  l2bin_input.c
  query_disc_wrapper.f90
  dataday.c
)

target_link_libraries(l2bin
  l2
  bin
  gsl
  gslcblas
  chealpix
  healpix
  cfitsio
  chealpix
  healpix
  cfitsio
  nav
)

set_target_properties(l2bin PROPERTIES LINKER_LANGUAGE Fortran)
install (TARGETS l2bin DESTINATION bin/${OCSSW_ARCH})
