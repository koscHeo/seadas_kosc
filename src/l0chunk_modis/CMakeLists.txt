cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l0chunk_modis
  l0chunk_modis.c
)

include_directories (
  ${PGS_INC}
)

target_link_libraries(l0chunk_modis
  ${PGS_LIB}
  m
)

install (TARGETS l0chunk_modis DESTINATION bin/${OCSSW_ARCH})

