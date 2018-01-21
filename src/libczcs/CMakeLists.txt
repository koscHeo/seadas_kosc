cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(czcs STATIC
  czcscaleg.f
  czsuba.f
  satang.f
  subsatxy.f
)

target_link_libraries(czcs
  nav
  m
)
