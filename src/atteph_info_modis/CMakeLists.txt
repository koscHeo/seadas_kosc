cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(atteph_info_modis
  atteph_info_modis.c
)

include_directories (
  ${PGS_INC}
)

target_link_libraries(atteph_info_modis
  ${PGS_LIB}
  genutils
)

install (TARGETS atteph_info_modis DESTINATION bin/${OCSSW_ARCH})

