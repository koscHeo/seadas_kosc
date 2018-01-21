cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_library (netcdfutils STATIC
  ncdf_utils.c
  nc_gridutils.c
  nc4utils.c
  cdl_utils.cpp
)

target_link_libraries(netcdfutils
  piutils
  netcdf
  ${HDF5_LIBRARIES}
  z
)
