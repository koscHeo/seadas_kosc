cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_library (dfutils STATIC
  wrapper.c
  scale_sds.c
)

target_link_libraries(dfutils
  hdf4utils
  netcdfutils
  hdf5utils
)
