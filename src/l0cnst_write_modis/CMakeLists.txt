cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l0cnst_write_modis
  l0cnst_write_modis.c
)

include_directories (
  ${PGS_INC}
)

target_link_libraries(l0cnst_write_modis
  ${PGS_LIB}
  hdf5utils
  hdf4utils
  dl
)

install (TARGETS l0cnst_write_modis DESTINATION bin/${OCSSW_ARCH})

