cmake_minimum_required(VERSION 2.6)

project(OCSSW)

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSTDC_HEADERS -DHAVE_STRING_H")
add_executable(nccmp
  opt.c
  strlist.c
  xmalloc.c
  ncinfo.c
  nccmp.cpp
)

target_link_libraries(nccmp
  netcdf
  ${HDF5_LIBRARIES}
)

install (TARGETS nccmp DESTINATION bin3/${OCSSW_ARCH})
