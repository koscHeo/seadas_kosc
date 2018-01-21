# This test directive runs the regression tests for the l3bindump program
# The location of the source files is $OCSSWROOT/test/l3bindump
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
# 
################################################################################
# Only include these test directives if the OCSSW test tree exists

if (OCSSWTestDir_FOUND)

################################################################################
# l3bindump tests
################################################################################
    set(testdir $ENV{OCSSWROOT}/test/l3bindump)
    file(GLOB files "${testdir}/*.par")

    foreach(filename ${files})

        GET_FILENAME_COMPONENT(parfile ${filename} NAME)
        STRING(REGEX REPLACE ".par" "" outfile ${parfile})

        add_test(NAME "l3bindump_${outfile}-test"
                 WORKING_DIRECTORY ${testdir}
                 COMMAND l3bindump par=${parfile})

        add_test(NAME "l3bindump_${outfile}-check"
                 WORKING_DIRECTORY ${testdir}
                 COMMAND diff ${outfile} output/${outfile})

    endforeach(filename)

endif(OCSSWTestDir_FOUND)
