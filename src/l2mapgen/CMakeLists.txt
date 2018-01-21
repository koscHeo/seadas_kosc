cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l2mapgen
  main_l2mapgen.c
  l2mapgen_input.c
  mipolyutil.c
)

target_link_libraries(l2mapgen
  l2
  geotiff
  tiff
  png
)

install (TARGETS l2mapgen DESTINATION bin/${OCSSW_ARCH})

