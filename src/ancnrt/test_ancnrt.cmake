# This test directive runs the regression tests for the ancnrt program
# The location of the source files is $OCSSWROOT/test/ancnrt
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
#
#  Note that on the hdiff tests, only data (-d option) will return a 0 exit 
#  status, to be sure the global attributes have not changed exceot for the 
#  processing time and inputs, run hdiff without the -d on.
#
if (OCSSWTestDir_FOUND)
add_test(NAME "ancnrt-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
	COMMAND $ENV{OCSSWROOT}/test/ancnrt/ancnrt/met_ing_forecast1.scr ../ancnrt/gdas1.2015072.t00z.pgrb2.1p00.f000)

add_test(NAME "ancnrt-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND hdiff -d ancnrt/N201507200_MET_NCEP_0360x0181_f000.hdf output/N201507200_MET_NCEP_0360x0181_f000.hdf)
#
#  test the forecast time feature
add_test(NAME "ancnrt-gfs_f12-test"
        WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
        COMMAND $ENV{OCSSWROOT}/test/ancnrt/ancnrt/met_ing_forecast1.scr ../ancnrt/gfs.2015072.t00z.pgrb2.1p00.f012)

add_test(NAME "ancnrt-gfs_f12-check"
        WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
        COMMAND hdiff -d ancnrt/N201507212_MET_NCEP_0360x0181_f012.hdf output/N201507212_MET_NCEP_0360x0181_f012.hdf)
#
#  test the standard way of getting a MET file, just an analysis and no ozone
add_test(NAME "ancnrt-std_anl-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
        COMMAND $ENV{OCSSWROOT}/test/ancnrt/ancnrt/met_ingest2.scr ../ancnrt/gdas1.2015072.t00z.pgrb2.1p00.f000)

add_test(NAME "ancnrt-std_anl-check"
        WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
        COMMAND hdiff -d ancnrt/N201507200_MET_NCEP_6h.hdf output/N201507200_MET_NCEP_6h.hdf)
#
#  test ancnrt with the toast ozone
add_test(NAME "ancnrt-toast_oz-test"
        WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
        COMMAND $ENV{OCSSWROOT}/test/ancnrt/ancnrt/tovs_ingest.scr ../ancnrt/TOAST16_050417.GRB)

add_test(NAME "ancnrt-toast_oz-check"
        WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
        COMMAND hdiff -d ancnrt/S20051070001223_TOVS.OZONE output/S20051070001223_TOVS.OZONE)
endif(OCSSWTestDir_FOUND)
