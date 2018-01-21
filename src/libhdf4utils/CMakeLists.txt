cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_library (hdf4utils STATIC
  hdf_utils.c
  hdf4_utils.c
)

target_link_libraries(hdf4utils
  piutils
  ${HDF4_LIBRARIES}
)
