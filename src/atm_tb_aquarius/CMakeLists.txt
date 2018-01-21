cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(atm_tb_aquarius
  atm_tb_aquarius.f
)

target_link_libraries(atm_tb_aquarius
  aquarius
  hdf5utils
  hdf4utils
)

set_target_properties(atm_tb_aquarius PROPERTIES LINKER_LANGUAGE Fortran)

if (${Fortran_COMPILER_NAME} MATCHES "ifort")
   set (CMAKE_Fortran_LINK_FLAGS "")
endif (${Fortran_COMPILER_NAME} MATCHES "ifort")

install (TARGETS atm_tb_aquarius DESTINATION bin/${OCSSW_ARCH})

