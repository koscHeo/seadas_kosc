cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l1bgen_modisa
  L1B.c
  Preprocess.c
  L1B_Tables.c
  Granule.c
  HDF_Lib.c
  Metadata.c
  Reflective_Cal.c
  Emissive_Cal.c
  L1B_Setup.c
)

add_definitions(-DALLOW_MISSING_GRANULES -DNOCHECKLUT -DSEADAS)

include_directories (
  ${PGS_INC}
)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
  set(CMAKE_EXE_LINKER_FLAGS "-Wl,-stack_size -Wl,0x8000000")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
 

target_link_libraries(l1bgen_modisa
  ${PGS_LIB}
  ${HDFEOS_LIB}
  hdf5utils
  hdf4utils
  dl
)

install (TARGETS l1bgen_modisa DESTINATION bin/${OCSSW_ARCH})

