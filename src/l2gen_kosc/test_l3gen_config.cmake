# This test directive runs the regression tests for the l3gen program
# The location of the source files is $OCSSWROOT/test/l3gen
#
# Each test consists of two checks: 
#    1) process the parameter file to generate the output file 
#    2) compare the new output file to the baseline version via nccmp
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
################################################################################
# Only include these test directives if the OCSSW test tree exists

if (OCSSWTestDir_FOUND)

################################################################################
# chl creation - Aqua
################################################################################
add_test(NAME "l3gen_A2003060.L3b_DAY_CHL-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND l3gen par=A2003060.L3b_DAY_CHL.par)

add_test(NAME "l3gen_A2003060.L3b_DAY_CHL-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE}  A2003060.L3b_DAY_CHL.nc output/A2003060.L3b_DAY_CHL.nc)

################################################################################
# Hu chl creation - SeaWiFS
################################################################################
add_test(NAME "l3gen_S2002002.L3b_DAY_RRS_CHLHU-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND l3gen par=S2002002.L3b_DAY_RRS_CHLHU.par)

add_test(NAME "l3gen_S2002002.L3b_DAY_RRS_CHLHU-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE}  S2002002.L3b_DAY_RRS_CHLHU.nc output/S2002002.L3b_DAY_RRS_CHLHU.nc)

################################################################################
# (G)IOP creation - SeaWiFS
################################################################################
add_test(NAME "l3gen_S2002002.L3b_DAY_RRS_GIOP-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND l3gen par=S2002002.L3b_DAY_RRS_GIOP.par)

add_test(NAME "l3gen_S2002002.L3b_DAY_GIOP-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE}  S2002002.L3b_DAY_RRS_GIOP.nc output/S2002002.L3b_DAY_RRS_GIOP.nc)

################################################################################
# test out the Virtual Constellation stuff - HDF4
################################################################################

# modify the $ENV{OCTEST_TOLERANCE} for hdiff output
set(TOLERANCE '')
if(DEFINED ENV{OCTEST_TOLERANCE})
  STRING(REGEX REPLACE "-T" "-p" TOLERANCE $ENV{OCTEST_TOLERANCE})
endif(DEFINED ENV{OCTEST_TOLERANCE})

# Generate the Rrs_vc product
add_test(NAME "l3gen_A2006001.L3b_DAY_RRS_small_vc-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND l3gen par=A2006001.L3b_DAY_RRS_small_vc.par)

if(DEFINED ENV{OCTEST_TOLERANCE})
  add_test(NAME "l3gen_A2006001.L3b_DAY_RRS_small_vc-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND hdiff -e 10 ${TOLERANCE} A2006001.L3b_DAY_RRS_small_vc.L3b output/A2006001.L3b_DAY_RRS_small_vc.L3b)
else(DEFINED ENV{OCTEST_TOLERANCE})
  add_test(NAME "l3gen_A2006001.L3b_DAY_RRS_small_vc-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND hdiff -e 10 A2006001.L3b_DAY_RRS_small_vc.L3b output/A2006001.L3b_DAY_RRS_small_vc.L3b)
endif(DEFINED ENV{OCTEST_TOLERANCE})

# Generate the VC product
add_test(NAME "l3gen_A2006001.L3b_DAY_RRS_small_vc2-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND l3gen par=A2006001.L3b_DAY_RRS_small_vc2.par)

if(DEFINED ENV{OCTEST_TOLERANCE})
  add_test(NAME "l3gen_A2006001.L3b_DAY_RRS_small_vc2-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND hdiff -e 10 ${TOLERANCE}  A2006001.L3b_DAY_RRS_small_vc2.L3b output/A2006001.L3b_DAY_RRS_small_vc2.L3b)
else(DEFINED ENV{OCTEST_TOLERANCE})
  add_test(NAME "l3gen_A2006001.L3b_DAY_RRS_small_vc2-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND hdiff -e 10 A2006001.L3b_DAY_RRS_small_vc2.L3b output/A2006001.L3b_DAY_RRS_small_vc2.L3b)
endif(DEFINED ENV{OCTEST_TOLERANCE})

################################################################################
# test out the Virtual Constellation stuff - netCDF4
################################################################################

# Generate the Rrs_vc product
add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc.nc-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND l3gen par=A2003060.L3b_DAY_RRS_small_vc.nc.par)

add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc.nc-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE}  A2003060.L3b_DAY_RRS_small_vc.L3b.nc output/A2003060.L3b_DAY_RRS_small_vc.L3b.nc )

# Generate the VC product
add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc2.nc-test"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND l3gen par=A2003060.L3b_DAY_RRS_small_vc2.nc.par)

add_test(NAME "l3gen_A2003060.L3b_DAY_RRS_small_vc2.nc-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3gen
    COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version $ENV{OCTEST_TOLERANCE}  A2003060.L3b_DAY_RRS_small_vc.L3b.nc output/A2003060.L3b_DAY_RRS_small_vc.L3b.nc )


endif(OCSSWTestDir_FOUND)
