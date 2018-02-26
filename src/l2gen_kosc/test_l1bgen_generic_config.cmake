# This test directive runs the regression tests for the l1bgen_generic program
# The location of the source files is $OCSSWROOT/test/l2gen
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
#
# There is a directory for each mission, and this script will process
# the parameter files that match the regex "<mission letter>*.L1B*.par" 
# contained in those directories 
# 
# The comparision is made with the nccmp program using the following options:
# -m : compare variable attribute metadata
# -g : compare global attribute metadata
# -d : compare the data
# -f : continue checks (default is to stop on any difference)
# -C 10 : stop comparisons after 10 differences are found
# --globalex='date_created,software_version' : ignore global attribute 
#             differences we expect
#
# Optionally, include a tolerance on data comparisions if the $OCTEST_TOLERANCE 
#        environment variable was set prior to running cmake.  The form of this
#         variable is '-T %diff' or '-t abs diff'
# 
#
#  To add a new test, include the source file(s) and parameter file in the 
#  appropriate test directory.  Currently, the following mission directories 
#  are searched:
#     czcs
#     octs
#     seawifs
#     aqua
#     terra
#     meris
#     viirs
#     hico
#     goci
#     oli
#     mos
#     ocm
#     ocm2
#     osmi
#
################################################################################
# Only inclue these test directives if the OCSSW test tree exists

if (OCSSWTestDir_FOUND)

################################################################################
# loop through l1bgen_generic tests
################################################################################
set(MISSIONS 
     czcs
     octs
     seawifs
     aqua
     terra
     meris
     viirs
     hico
     goci
     oli
     mos
     ocm
     ocm2
     osmi
)
foreach(mission ${MISSIONS})
    file(GLOB files "$ENV{OCSSWROOT}/test/l2gen/${mission}/*L1B*.par")
    foreach(filename ${files})

      GET_FILENAME_COMPONENT(L1BPAR ${filename} NAME)
      STRING(REGEX REPLACE ".par" ".nc" L1BFILE ${L1BPAR})

      add_test(NAME "l1bgen_generic_${L1BFILE}-test"
             WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l2gen/${mission}
             COMMAND l1bgen_generic par=${L1BPAR})

      add_test(NAME "l1bgen_generic_${L1BFILE}-check"
             WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l2gen
             COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE} ${mission}/${L1BFILE} output/${L1BFILE})

    endforeach(filename)
endforeach(mission)

endif(OCSSWTestDir_FOUND)
