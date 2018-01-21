cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(map STATIC
  get_l3m.c
  put_l3m.c
)

target_link_libraries(map
  dfutils
)
