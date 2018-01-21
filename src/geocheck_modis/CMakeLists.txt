cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(geocheck_modis
  geocheck_modis.c
)

target_link_libraries(geocheck_modis
  ${HDF4_LIBRARIES}
  m
)

install (TARGETS geocheck_modis DESTINATION bin/${OCSSW_ARCH})

