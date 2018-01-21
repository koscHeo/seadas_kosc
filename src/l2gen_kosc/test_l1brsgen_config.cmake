# This test directive runs the regression tests for the l1brsgen program
# The location of the source files is $OCSSWROOT/test/l1brsgen
#
################################################################################
# Only include these test directives if the OCSSW test tree exists

if (OCSSWTestDir_FOUND)

################################################################################
# PNG creation - Aqua
################################################################################
set(FUZZ "0%")
if(DEFINED ENV{OCTEST_TOLERANCE})
  set(FUZZ "6%")
endif(DEFINED ENV{OCTEST_TOLERANCE})

add_test(NAME "l1brsgen_A2008080195500.L1B_BRS_PNG-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND l1brsgen par=A2008080195500.L1B_BRS_PNG.par)

add_test(NAME "l1brsgen_A2008080195500.L1B_BRS_PNG-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND compare -metric AE -fuzz ${FUZZ} A2008080195500.L1B_BRS.png output/A2008080195500.L1B_BRS.png /dev/null )


################################################################################
#  HDF creation - Aqua
################################################################################
# modify the $ENV{OCTEST_TOLERANCE} for hdiff output

add_test(NAME "l1brsgen_A2008080195500.L1B_BRS-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND l1brsgen par=A2008080195500.L1B_BRS.par)

set(TOLERANCE '')
if(DEFINED ENV{OCTEST_TOLERANCE})
  STRING(REGEX REPLACE "-T" "-p" TOLERANCE $ENV{OCTEST_TOLERANCE})
  add_test(NAME "l1brsgen_A2008080195500.L1B_BRS-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND hdiff -e 10 ${TOLERANCE} A2008080195500.L1B_BRS.hdf output/A2008080195500.L1B_BRS.hdf)
else(DEFINED ENV{OCTEST_TOLERANCE})
  add_test(NAME "l1brsgen_A2008080195500.L1B_BRS-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND hdiff -e 10 A2008080195500.L1B_BRS.hdf output/A2008080195500.L1B_BRS.hdf)
endif(DEFINED ENV{OCTEST_TOLERANCE})


endif(OCSSWTestDir_FOUND)
