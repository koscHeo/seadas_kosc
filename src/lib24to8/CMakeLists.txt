cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(24to8 STATIC
  error.c
  image.c
  quantize.c
  signature.c
)

