# This test directive runs the regression tests for the l1agen_seawifs program
# The location of the source files is $OCSSWROOT/test/l1agen_seawifs
#
################################################################################
# Only inclue these test directives if the OCSSW test tree exists

if (OCSSWTestDir_FOUND)

################################################################################
# GAC
################################################################################

add_test(NAME "l1agen_S2000016180835.L0_GAC-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1agen_seawifs/output
    COMMAND l1agen_seawifs -n sdps ../S2000016180835.L0_GAC)

add_test(NAME "l1agen_S2000016180835.L0_GAC-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1agen_seawifs
    COMMAND hdiff -e 10 -s S2000016180835.L0_GAC output/S2000016180835.L0_GAC )

################################################################################
# HRPT
################################################################################
add_test(NAME "l1agen_S2000011192121.L0_HNSG-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1agen_seawifs/output
    COMMAND l1agen_seawifs -n sdps ../S2000011192121.L0_HNSG)

add_test(NAME "l1agen_S2000011192121.L0_HNSG-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l1agen_seawifs
    COMMAND hdiff -e 10 -s S2000011192121.L1A_HNSG output/S2000011192121.L1A_HNSG )

endif(OCSSWTestDir_FOUND)
