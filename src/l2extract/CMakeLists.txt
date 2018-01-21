cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l2extract
  main_l2extract.c
  l2extract_netcdf.c
)

include_directories (
  ${CMAKE_SOURCE_DIR}/src/l2gen
)

target_link_libraries(l2extract
  l2
)

install (TARGETS l2extract DESTINATION bin/${OCSSW_ARCH})

