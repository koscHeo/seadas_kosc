cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(mk_aquarius_ancillary_data
  mk_aquarius_ancillary_data.f
)

target_link_libraries(mk_aquarius_ancillary_data
  aqsss
  swfnav
  netcdff
  hdf5utils
  hdf4utils
  hdf5_hl
)

set_target_properties(mk_aquarius_ancillary_data PROPERTIES LINKER_LANGUAGE Fortran)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,-stack_size -Wl,0x20000000")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if (${Fortran_COMPILER_NAME} MATCHES "ifort")
  set (CMAKE_Fortran_LINK_FLAGS "")
endif (${Fortran_COMPILER_NAME} MATCHES "ifort")

install (TARGETS mk_aquarius_ancillary_data DESTINATION bin/${OCSSW_ARCH})

