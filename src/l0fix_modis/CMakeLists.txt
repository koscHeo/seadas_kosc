cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l0fix_modis
  l0fix_modis.c
)

include_directories (
  ${PGS_INC}
)

target_link_libraries(l0fix_modis
  ${PGS_LIB}
  m
)

install (TARGETS l0fix_modis DESTINATION bin/${OCSSW_ARCH})

