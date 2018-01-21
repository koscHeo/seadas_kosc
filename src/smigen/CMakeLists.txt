cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(smigen
  smigen.cpp
  smigen_input.c
  put_smi.cpp
)

add_executable(fixrgb
  fixrgb.c
)

target_link_libraries(smigen
  bin++
  bin
  dfutils
  ${EXTRA_OS_LIBS}
)

install (TARGETS smigen fixrgb DESTINATION bin/${OCSSW_ARCH})
