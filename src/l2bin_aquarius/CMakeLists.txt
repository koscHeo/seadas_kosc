cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l2bin_aquarius
  l2bin_aquarius.cpp
)

target_link_libraries(l2bin_aquarius
  bin++
  aquarius
  gsl
  gslcblas
  netcdf
  hdf5_hl
  ${EXTRA_OS_LIBS}
)

install (TARGETS l2bin_aquarius DESTINATION bin/${OCSSW_ARCH})

