# This test directive runs the regression tests for the tec program
# The location of the source files is $OCSSWROOT/test/ancnrt
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
#
#
if (OCSSWTestDir_FOUND)
add_test(NAME "tec-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
#	COMMAND $ENV{OCSSWROOT}/cbuild/src/ancnrt/tec ../tec/igsg0800.15i)
	COMMAND tec ../tec/igsg0800.15i)

add_test(NAME "tec-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND nccmp -m -g -d -f -C 10 $ENV{OCTEST_TOLERANCE} tec/N201508000_TEC_IGS_24h.h5 output/N201508000_TEC_IGS_24h.h5)
#
endif(OCSSWTestDir_FOUND)
