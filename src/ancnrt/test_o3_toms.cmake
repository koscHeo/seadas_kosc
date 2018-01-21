# This test directive runs the regression tests for the o3_tome program
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
add_test(NAME "o3_toms-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/o3_toms
	COMMAND o3_toms L3_ozone_omi_20150101.txt ../output L3_ozone_omi_20141231.txt)
	#COMMAND $ENV{OCSSWROOT}/cbuild/src/ancnrt/o3_toms L3_ozone_omi_20150101.txt ../output L3_ozone_omi_20141231.txt)

add_test(NAME "o3_toms-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND hdiff -d o3_toms/N201500100_O3_AURAOMI_24h.hdf output/N201500100_O3_AURAOMI_24h.hdf)
#
endif(OCSSWTestDir_FOUND)
