cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(piutils STATIC
  productInfo.c
  sensorInfo.c
  rdsensorinfo.c
)

target_link_libraries(piutils
  genutils
)

