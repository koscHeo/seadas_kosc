cmake_minimum_required(VERSION 2.8)

PROJECT(OCSSW)

ADD_LIBRARY(timeutils STATIC
  addmsec.c
  date2ydmsec.c
  gmt_offset.c
  now.c
  unix2yds.c
  unix2ydhmsf.c
  unix2ymds.c
  unix2ymdhms.c
  yd2md.c
  ydhms2ydmsec.c
  ydhmsf.c
  yds2unix.c
  ymds2unix.c
  zulu2unix.c
  isleap.c
  leapsecond.c
  get_time.c
  jdate.c
  ccsds_to_yds.c
  yds2tai.c
  jday.c
)

target_link_libraries(timeutils
  m
)
