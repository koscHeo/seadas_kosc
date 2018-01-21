cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_library (hdf5utils STATIC
  hdf5util.cpp
  hdf5_utils.f
)

target_link_libraries(hdf5utils
  piutils
  hdf5_fortran
  ${HDF5_LIBRARIES}
)
