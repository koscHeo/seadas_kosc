# This test directive runs the regression tests for the l3mapgen program
# The location of the source files is $OCSSWROOT/test/l3mapgen
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
# l3mapgen tests
################################################################################
set(FUZZ "0%")
if(DEFINED ENV{OCTEST_TOLERANCE})
  set(FUZZ "6%")
endif(DEFINED ENV{OCTEST_TOLERANCE})

# test for NetCDF output
add_test(NAME "l3mapgen_netcdf-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=out.nc.par)

add_test(NAME "l3mapgen_netcdf-check"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} out.nc output/out.nc)


# test for HDF4 SMI output
add_test(NAME "l3mapgen_hdf4-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=out.hdf.par)

add_test(NAME "l3mapgen_hdf4-check"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND hdiff -e 10 out.hdf output/out.hdf)


# test PNG output
add_test(NAME "l3mapgen_png-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=out.png.par)

add_test(NAME "l3mapgen_png-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
    COMMAND compare -metric AE -fuzz ${FUZZ} out.png output/out.png /dev/null )


# test PPM output
add_test(NAME "l3mapgen_ppm-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=out.ppm.par)

add_test(NAME "l3mapgen_ppm-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
    COMMAND compare -metric AE -fuzz ${FUZZ} out.ppm output/out.ppm /dev/null )


# need to add a test for the geotiff metadata tags
#
add_test(NAME "l3mapgen_tiff-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=out.tiff.par)

add_test(NAME "l3mapgen_tiff-check"
    WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
    COMMAND tiffcmp out.tiff output/out.tiff )


# test a different resolution
add_test(NAME "l3mapgen_resolution-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=resolution.nc.par)

add_test(NAME "l3mapgen_resolution-check"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} resolution.nc output/resolution.nc)


# test a mollweide projection projection
# rotate the central meridian to 180
add_test(NAME "l3mapgen_projection-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=projection.nc.par)

add_test(NAME "l3mapgen_projection-check"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} projection.nc output/projection.nc)


# test multiple product output
# area interpolation
add_test(NAME "l3mapgen_multiple-test"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND l3mapgen par=multiple.nc.par)

add_test(NAME "l3mapgen_multiple-check"
         WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l3mapgen
         COMMAND nccmp -m -g -d -f -C 10 -G date_created,software_version,_lastModified $ENV{OCTEST_TOLERANCE} multiple.nc output/multiple.nc)




endif(OCSSWTestDir_FOUND)
