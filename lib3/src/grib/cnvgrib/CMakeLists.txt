cmake_minimum_required(VERSION 2.6)

project(GRIB)

add_executable(cnvgrib
  gds2gdt.f
  pds2pdt.f
  pds2pdtens.f
  cnvgrib.f
  cnv12.f
  cnv22.f
  cnv21.f
  gdt2gds.f
  makepds.f
  putgbexn.f
  makepdsens.f
  setbit.f
)

target_link_libraries(cnvgrib
  g2
)

set_target_properties(cnvgrib PROPERTIES LINKER_LANGUAGE Fortran)

if (${Fortran_COMPILER_NAME} MATCHES "ifort")
  set (CMAKE_Fortran_LINK_FLAGS "")
endif (${Fortran_COMPILER_NAME} MATCHES "ifort")

install (TARGETS cnvgrib DESTINATION bin3/$ENV{OCSSW_ARCH})

