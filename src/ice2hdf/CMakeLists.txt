cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(ice2hdf
  ice2hdf.c
)

target_link_libraries(ice2hdf 
  genutils
  hdfutils
)

install (TARGETS ice2hdf DESTINATION bin/${OCSSW_ARCH})
