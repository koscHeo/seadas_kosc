# This test directive runs the regression tests for the swh program
# The location of the source files is $OCSSWROOT/test/ancnrt
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
#
#
if (OCSSWTestDir_FOUND)
add_test(NAME "swh-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/swh
#	COMMAND $ENV{OCSSWROOT}/cbuild/src/ancnrt/swh swh.dat ../output/N201509900_SWH_NCEP_6h.h5)
	COMMAND swh swh.dat ../output/N201509900_SWH_NCEP_6h.h5)

add_test(NAME "swh-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND nccmp -m -g -d -f -C 10 $ENV{OCTEST_TOLERANCE} swh/N201509900_SWH_NCEP_6h.h5 output/N201509900_SWH_NCEP_6h.h5)
#
endif(OCSSWTestDir_FOUND)
