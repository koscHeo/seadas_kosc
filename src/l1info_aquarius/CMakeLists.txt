cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l1info_aquarius
  l1info_aquarius.cpp
)

target_link_libraries(l1info_aquarius
  aqgeo
  gsl
  gslcblas
  ${IFORT_CXX_LIBS}
  ${EXTRA_OS_LIBS}
)

if (${Fortran_COMPILER_NAME} MATCHES "ifort")
    set_target_properties(l1info_aquarius PROPERTIES LINKER_LANGUAGE Fortran)
endif (${Fortran_COMPILER_NAME} MATCHES "ifort")

install (TARGETS l1info_aquarius DESTINATION bin/${OCSSW_ARCH})

