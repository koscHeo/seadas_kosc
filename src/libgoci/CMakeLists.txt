cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(goci STATIC
  goci.c
  goci_slot.c
)

target_link_libraries(goci
  proj
  hdf5_fortran
  ${HDF5_LIBRARIES}
  m
)
