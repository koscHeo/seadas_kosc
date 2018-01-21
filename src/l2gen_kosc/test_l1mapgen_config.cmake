# This test directive runs the regression tests for the l1mapgen program
# The location of the source files is $OCSSWROOT/test/l1mapgen
#
################################################################################
# Only inclue these test directives if the OCSSW test tree exists

if (OCSSWTestDir_FOUND)

set(FUZZ "0%")
if(DEFINED ENV{OCTEST_TOLERANCE})
  set(FUZZ "6%")
endif(DEFINED ENV{OCTEST_TOLERANCE})

################################################################################
# PPM creation - Aqua
################################################################################

add_test(NAME "l1mapgen_A2008080195500.L1B_MAP_PNM-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND l1mapgen par=A2008080195500.L1B_MAP_PNM.par)

add_test(NAME "l1mapgen_A2008080195500.L1B_MAP_PNM-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND compare -metric AE -fuzz ${FUZZ} A2008080195500.L1B_MAP.ppm output/A2008080195500.L1B_MAP.ppm output/A2002365234500.L1B_MAP.diff.ppm )

################################################################################
# PNG creation - Aqua
################################################################################
add_test(NAME "l1mapgen_A2008080195500.L1B_MAP_PNG-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND l1mapgen par=A2008080195500.L1B_MAP_PNG.par)

add_test(NAME "l1mapgen_A2008080195500.L1B_MAP_PNG-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND compare -metric AE -fuzz ${FUZZ} A2008080195500.L1B_MAP.png output/A2008080195500.L1B_MAP.png output/A2002365234500.L1B_MAP.diff.png )

################################################################################
# GeoTIFF creation - Aqua
################################################################################
add_test(NAME "l1mapgen_A2008080195500.L1B_MAP_TIF-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND l1mapgen par=A2008080195500.L1B_MAP_TIF.par)

add_test(NAME "l1mapgen_A2008080195500.L1B_MAP_TIF-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND compare -metric AE -fuzz ${FUZZ} A2008080195500.L1B_MAP.tiff output/A2008080195500.L1B_MAP.tiff output/A2002365234500.L1B_MAP.diff.tiff)

################################################################################
# PNG creation - SeaWiFS
################################################################################
add_test(NAME "l1mapgen_S2002079071209.L1A_GAC_MAP_PNG-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND l1mapgen par=S2002079071209.L1A_GAC_MAP_PNG.par)

add_test(NAME "l1mapgen_S2002079071209.L1A_GAC_MAP_PNG-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1brsgen
    COMMAND compare -metric AE -fuzz ${FUZZ} S2002079071209.L1A_GAC_MAP.png output/S2002079071209.L1A_GAC_MAP.png output/S1998199173926.L1A_GAC_MAP.diff.png )

endif(OCSSWTestDir_FOUND)
