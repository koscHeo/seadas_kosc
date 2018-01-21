cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(hico_L0_hdr_echo
  hico_L0_hdr_echo.c
)

#if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#  include_directories (
#    /usr/include/machine
#  )
#endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")                                    

add_executable(time_offset
  time_offset.c
)

install (TARGETS hico_L0_hdr_echo time_offset DESTINATION bin/${OCSSW_ARCH})

