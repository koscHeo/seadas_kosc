# This test directive runs the regression tests for the anc_seaice program
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
#file(REMOVE $ENV{OCSSWROOT}/test/ancnrt/output/N201417500_SEAICE_NCEP_24h.hdf)
add_test(NAME "anc_seaice-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
	COMMAND $ENV{OCSSWROOT}/test/ancnrt/anc_seaice/ice_ncep_ingest2.scr ../anc_seaice/eng5min.grib2.20140624)

add_test(NAME "anc_seaice-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND hdiff -d anc_seaice/N201417500_SEAICE_NCEP_24h.hdf output/N201417500_SEAICE_NCEP_24h.hdf)
endif(OCSSWTestDir_FOUND)
