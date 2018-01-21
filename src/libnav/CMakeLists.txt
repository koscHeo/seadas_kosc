cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(nav STATIC
  cdata.f
  crossp.f
  esdist.f
  filenv.f
  gha2000.f
  jd.f
  l_sun.f
  ocorient.f
  sunangs.f
  get_zenaz.c
  compute_alpha.c
)

target_link_libraries(nav
  gsl
  gslcblas
  m
)
