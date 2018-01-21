cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l3bin
  l3bin.cpp
  l3bin_input.c
)

add_executable(l3binmerge
  l3binmerge.cpp
  l3bin_input.c
)

add_executable(l3rebin_meris
  l3rebin_meris.cpp
)

target_link_libraries(l3bin 
  piutils
  bin++
  bin
  netcdf
  ${HDF5_LIBRARIES}
  chealpix
  healpix
  cfitsio
  genutils
  ${EXTRA_OS_LIBS}
)

target_link_libraries(l3binmerge 
  bin++
  bin
  genutils
  netcdf
  ${HDF5_LIBRARIES}
  ${EXTRA_OS_LIBS}
)

target_link_libraries(l3rebin_meris 
  bin++
  bin
  genutils
  netcdf
  ${HDF5_LIBRARIES}
  ${EXTRA_OS_LIBS}
)

install (TARGETS l3bin l3binmerge l3rebin_meris DESTINATION bin/${OCSSW_ARCH})
