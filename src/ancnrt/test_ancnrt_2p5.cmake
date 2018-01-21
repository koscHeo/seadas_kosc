# This test directive runs the regression tests for the ancnrt_2p5 program
#  it currently deals with the reanalysis 2 2.5 deg grid and synthesizes 
#  the 10 m wind
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
add_test(NAME "ancnrt_2p5-test"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt/output
	COMMAND $ENV{OCSSWROOT}/test/ancnrt/ancnrt/met_ingest_ncep_reanal2.scr ../ancnrt/pgb.201406 NCEPpreR2)

add_test(NAME "ancnrt_2p5-check"
	WORKING_DIRECTORY $ENV{OCSSWROOT}/test/ancnrt
	COMMAND hdiff -d ancnrt/N201415200_MET_NCEPpreR2_6h.hdf output/N201415200_MET_NCEPpreR2_6h.hdf)
#
endif(OCSSWTestDir_FOUND)
