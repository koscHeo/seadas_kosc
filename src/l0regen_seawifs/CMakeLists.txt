cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l0regen_seawifs
  main_l0regen.c
)

target_link_libraries(l0regen_seawifs
  hdf4utils
)

install (TARGETS l0regen_seawifs DESTINATION bin/${OCSSW_ARCH})
