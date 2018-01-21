cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l2brsgen
  main_l2brsgen.c
  l2brsgen_input.c
  put_l2brs.c
)

target_link_libraries(l2brsgen
  l2
  png
)

install (TARGETS l2brsgen DESTINATION bin/${OCSSW_ARCH})
