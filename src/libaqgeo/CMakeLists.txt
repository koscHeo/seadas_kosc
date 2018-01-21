cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(aqgeo STATIC
  moon_vector.f
  geolocation.f
  geolocation_routines.f
  fd_date_2000.f
  pre_nut_routines.f
  sun_vector.f
)

target_link_libraries(aqgeo
  aquarius
  m
)
