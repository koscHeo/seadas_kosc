cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(smitoppm
  smitoppm.c
)

target_link_libraries(smitoppm
  map
  netcdf
  ${HDF5_LIBRARIES}
)

install (TARGETS smitoppm DESTINATION bin/${OCSSW_ARCH})

