cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l2gen_aquarius
  l2gen_aquarius.cpp
  get_prelaunch_calibration_coeffs.f
  count_to_ta.f
  l2gen_aquarius_module.f
  get_rfi_auxiliary_data.f
  rfi_detection.f
  SolarFluxFlag.c
  gsw_saar.c
  gsw_oceanographic_toolbox.c
)

add_executable(l2gen_sm
  l2gen_sm.cpp
)

add_executable(aq_sca_main
  aq_sca_main.f90
  tb2vsm.f90
)

include_directories (
  ${CMAKE_SOURCE_DIR}/include/aquarius
)

target_link_libraries(l2gen_aquarius
  aqgeo
  aqsss
  aquarius
#  aqscatter
#  scatterutils
  gsl
  gslcblas 
  hdf5utils
  ${IFORT_CXX_LIBS}
  ${EXTRA_OS_LIBS}
)

target_link_libraries(l2gen_sm
   aquarius 
   hdf5utils 
   z
   ${EXTRA_OS_LIBS}
)

target_link_libraries(aq_sca_main
  hdf5_fortran 
  ${HDF5_LIBRARIES}
)

if(BUILD_JPL)
  add_executable(l2gen_scatter
    l2gen_scatter.f
  )

  target_link_libraries(l2gen_scatter
    aqscatter
    scatterutils
    aqsss
    hdf5utils
    hdf4utils
  )

  set_target_properties(l2gen_scatter PROPERTIES LINKER_LANGUAGE Fortran)
endif(BUILD_JPL)

set_target_properties(aq_sca_main PROPERTIES LINKER_LANGUAGE Fortran)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,-stack_size -Wl,0x20000000")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if (${Fortran_COMPILER_NAME} MATCHES "ifort")
  set_target_properties(l2gen_aquarius PROPERTIES LINKER_LANGUAGE Fortran)
  set (CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -DIFORT_MODULE_NAME_MANGLING")
endif (${Fortran_COMPILER_NAME} MATCHES "ifort")

install (TARGETS l2gen_aquarius l2gen_sm aq_sca_main DESTINATION bin/${OCSSW_ARCH})
if(BUILD_JPL)
  install (TARGETS l2gen_scatter DESTINATION bin/${OCSSW_ARCH})
endif(BUILD_JPL)
