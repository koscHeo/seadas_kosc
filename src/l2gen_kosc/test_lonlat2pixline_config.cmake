# This test directive runs the regression tests for the lonlat2pixline program
# The location of the source files is $OCSSWROOT/test/lonlat2pixline
#
################################################################################
# Only inclue these test directives if the OCSSW test tree exists

if (OCSSWTestDir_FOUND)

################################################################################
# SeaWiFS find point
################################################################################

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_point-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND lonlat2pixline -o output/S2002079134718.L1A_GAC.middle_point.txt ../l2gen/seawifs/S2002079134718.L1A_GAC.middle -34 -8 )

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_point-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u S2002079134718.L1A_GAC.middle_point.txt output/S2002079134718.L1A_GAC.middle_point.txt )

################################################################################
# SeaWiFS find box
################################################################################

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_box-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND lonlat2pixline -o output/S2002079134718.L1A_GAC.middle_box.txt ../l2gen/seawifs/S2002079134718.L1A_GAC.middle -34 -13 -30 -8 )

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_box-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u S2002079134718.L1A_GAC.middle_box.txt output/S2002079134718.L1A_GAC.middle_box.txt )

################################################################################
# SeaWiFS find point failure, exit=1
################################################################################

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_point-fail"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND test_exit_code.sh 1 lonlat2pixline ../l2gen/seawifs/S2002079134718.L1A_GAC.middle 0 0 )

################################################################################
# SeaWiFS find box failure, exit=1
################################################################################

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_box-fail"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND test_exit_code.sh 1 lonlat2pixline ../l2gen/seawifs/S2002079134718.L1A_GAC.middle 0 0 2 2 )

################################################################################
# SeaWiFS partial box extracted, exit=0
################################################################################

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_partial-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND lonlat2pixline -o output/S2002079134718.L1A_GAC.middle_partial.txt ../l2gen/seawifs/S2002079134718.L1A_GAC.middle -34 -8 50 50 )

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_partial-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u S2002079134718.L1A_GAC.middle_partial.txt output/S2002079134718.L1A_GAC.middle_partial.txt )

################################################################################
# SeaWiFS full box not extracted, exit=110
################################################################################

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_full-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND test_exit_code.sh 110 lonlat2pixline -F -o output/S2002079134718.L1A_GAC.middle_full.txt ../l2gen/seawifs/S2002079134718.L1A_GAC.middle -34 -8 50 50 )

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_full-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u S2002079134718.L1A_GAC.middle_full.txt output/S2002079134718.L1A_GAC.middle_full.txt )

################################################################################
# SeaWiFS whole file is inside box, exit=120
################################################################################

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_whole-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND test_exit_code.sh 120 lonlat2pixline -F -o output/S2002079134718.L1A_GAC.middle_whole.txt ../l2gen/seawifs/S2002079134718.L1A_GAC.middle -50 -50 50 50 )

add_test(NAME "lonlat2pixline_S2002079134718.L1A_GAC.middle_whole-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u S2002079134718.L1A_GAC.middle_whole.txt output/S2002079134718.L1A_GAC.middle_whole.txt )

################################################################################
# Aqua find point
################################################################################

add_test(NAME "lonlat2pixline_A2008080195500.L1B_LAC.subline_point-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND lonlat2pixline -o output/A2008080195500.L1B_LAC.subline_point.txt ../l2gen/aqua/A2008080195500.L1B_LAC.subline ../l2gen/aqua/A2008080195500.GEO.subline -92 -47 )

add_test(NAME "lonlat2pixline_A2008080195500.L1B_LAC.subline_point-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u A2008080195500.L1B_LAC.subline_point.txt output/A2008080195500.L1B_LAC.subline_point.txt )

################################################################################
# Aqua find box
################################################################################

add_test(NAME "lonlat2pixline_A2008080195500.L1B_LAC.subline_box-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND lonlat2pixline -o output/A2008080195500.L1B_LAC.subline_box.txt ../l2gen/aqua/A2008080195500.L1B_LAC.subline ../l2gen/aqua/A2008080195500.GEO.subline -92 -47 -91 -46 )

add_test(NAME "lonlat2pixline_A2008080195500.L1B_LAC.subline_box-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u A2008080195500.L1B_LAC.subline_box.txt output/A2008080195500.L1B_LAC.subline_box.txt )

################################################################################
# Meris find point
################################################################################

add_test(NAME "lonlat2pixline_MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_point-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND lonlat2pixline -o output/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_point.txt ../l2gen/meris/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix -40 42 )

add_test(NAME "lonlat2pixline_MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_point-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_point.txt output/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_point.txt )

################################################################################
# Meris find box
################################################################################

add_test(NAME "lonlat2pixline_MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_box-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND lonlat2pixline -o output/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_box.txt ../l2gen/meris/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix -40 41.2 -36.5 43 )

add_test(NAME "lonlat2pixline_MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_box-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/lonlat2pixline
    COMMAND diff -u MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_box.txt output/MER_RR__1PRACR20080917_123930_000026292072_00124_34246_0000.N1.subpix_box.txt )


endif(OCSSWTestDir_FOUND)


