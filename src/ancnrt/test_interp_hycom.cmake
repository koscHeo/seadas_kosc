# This test directive runs the regression tests for the interp_hycom program
# The location of the source files is $OCSSWROOT/test/ancnrt
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
#
#
if (OCSSWTestDir_FOUND)
add_test(NAME "interp_hycom-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/interp_hycom
#	COMMAND $ENV{OCSSWROOT}/cbuild/src/ancnrt/interp_hycom  archv.2015_001_00_3zs.nc ../output)
	COMMAND interp_hycom archv.2015_001_00_3zs.nc ../output)

add_test(NAME "interp_hycom-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND nccmp -m -g -d -f -C 10 $ENV{OCTEST_TOLERANCE} interp_hycom/N201500100_SALINITY_HYCOM_24h.h5 output/N201500100_SALINITY_HYCOM_24h.h5)
#
endif(OCSSWTestDir_FOUND)
