# This test directive runs the regression tests for the ancnrt_aquarius program
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
add_test(NAME "ancnrt_aquarius-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
	COMMAND $ENV{OCSSWROOT}/test/ancnrt/ancnrt_aquarius/aquarius_ing.scr ../ancnrt_aquarius/gdas1.t12z.pgrb2.1p00.f000 ../ancnrt_aquarius/gfs.t12z.pgrb2.0p50.f000)

add_test(NAME "ancnrt_aquarius-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND hdiff -d ancnrt_aquarius/N201509312_QMET_NCEP_6h output/N201509312_QMET_NCEP_6h)
#
endif(OCSSWTestDir_FOUND)
