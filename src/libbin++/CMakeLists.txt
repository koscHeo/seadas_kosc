cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(bin++ STATIC
  bin_io.cpp
  bin_util.cpp
  L3BinShape.cpp
  L3File.cpp
  L3BinShapeSMI.cpp
  L3FileSMI.cpp
)

target_link_libraries(bin++
  netcdf_c++4 
  netcdf
  ${HDF5_LIBRARIES}
  bin
)
