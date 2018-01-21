cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l1aextract_modis
  l1aextract_modis.c
)

target_link_libraries(l1aextract_modis
  ${HDF4_LIBRARIES}
)

install (TARGETS l1aextract_modis DESTINATION bin/${OCSSW_ARCH})

