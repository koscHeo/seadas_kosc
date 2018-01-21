cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l1asubset_modis
  l1asubset_modis.c
)

target_link_libraries(l1asubset_modis
  genutils
  ${HDF4_LIBRARIES}
)

install (TARGETS l1asubset_modis DESTINATION bin/${OCSSW_ARCH})
