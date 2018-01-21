cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l1aextract_seawifs
  main_l1aextract.c
  extract_sub.c
)

target_link_libraries(l1aextract_seawifs
  seawifs
  hdf4utils
)

install (TARGETS l1aextract_seawifs DESTINATION bin/${OCSSW_ARCH})
