# This test directive runs the regression tests for the solar_xray program
# The location of the source files is $OCSSWROOT/test/ancnrt
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
#
#
if (OCSSWTestDir_FOUND)
add_test(NAME "solar_xray-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
#	COMMAND $ENV{OCSSWROOT}/cbuild/src/ancnrt/solar_xray ../solar_xray/20150301_Gs_xr_1m.txt)
	COMMAND solar_xray ../solar_xray/20150301_Gs_xr_1m.txt)

add_test(NAME "solar_xray-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND nccmp -m -g -d -f -C 10 $ENV{OCTEST_TOLERANCE} solar_xray/N201506000_XRAY_GOES_24h.h5 output/N201506000_XRAY_GOES_24h.h5)
#
endif(OCSSWTestDir_FOUND)
