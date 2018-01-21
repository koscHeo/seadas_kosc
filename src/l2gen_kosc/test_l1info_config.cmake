# This test directive runs the regression tests for the l1info program
# The location of the source files is $OCSSWROOT/test/l2gen
#
# Each test consists of two checks: 
#    1) run l1info on the L1 file piped to a text file 
#    2) run a simple diff command on the text output versus the repository copy
#
# There is a directory for each mission, and this script will process
# the L1 files that match the regex "<mission letter>*.L1*" 
# contained in those directories 
# 
#  Currently, the following mission directories are searched:
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
# loop through l1info tests
################################################################################
set(MISSIONS 
     czcs
     octs
     seawifs
     aqua
     terra
     meris
     hico
     goci
#     mos
#     ocm
#     ocm2
     osmi
     viirs
#     oli
)
foreach(mission ${MISSIONS})
    set(L1TYPE "L1A")
    if(${mission} MATCHES aqua|terra)
      set(L1TYPE "L1B_LAC")
    endif(${mission} MATCHES aqua|terra)

    if(${mission} MATCHES meris)
      set(L1TYPE "N1")
    endif(${mission} MATCHES meris)

    if(${mission} MATCHES hico)
      set(L1TYPE "L1B_ISS")
    endif(${mission} MATCHES hico)

    if(${mission} MATCHES goci)
      set(L1TYPE "L1B")
    endif(${mission} MATCHES goci)

    if(${mission} MATCHES viirs)
      set(L1TYPE "SVM01")
    endif(${mission} MATCHES viirs)

    file(GLOB files "$ENV{OCSSWROOT}/test/l2gen/${mission}/*${L1TYPE}*")

    file(GLOB parfiles "$ENV{OCSSWROOT}/test/l2gen/${mission}/*${L1TYPE}*par")
    foreach(filename ${parfiles})
        list(REMOVE_ITEM files ${filename})
    endforeach(filename)

    file(GLOB infofiles "$ENV{OCSSWROOT}/test/l2gen/${mission}/*${L1TYPE}*info")
    foreach(filename ${infofiles})
        list(REMOVE_ITEM files ${filename})
    endforeach(filename)

    
    foreach(filename ${files})
      GET_FILENAME_COMPONENT(L1FILE ${filename} NAME)
      set(GEOFILE "")
      if(${mission} MATCHES aqua|terra)
        STRING(REGEX REPLACE "L1B_LAC" "GEO" GEOFILE ${L1FILE})
      endif(${mission} MATCHES aqua|terra) 

      if(${mission} MATCHES viirs)
        STRING(REGEX REPLACE "SVM01" "GMTCO" GEOFILE ${L1FILE})
      endif(${mission} MATCHES viirs) 

      add_test(NAME "l1info_${L1FILE}-test"
             WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l2gen/${mission}
             COMMAND l1info -i 50 -o ../output/${L1FILE}.l1info ${L1FILE} ${GEOFILE} )

      add_test(NAME "l1info_${L1FILE}-check"
             WORKING_DIRECTORY $ENV{OCSSWROOT}/test/l2gen
             COMMAND diff ${mission}/${L1FILE}.l1info output/${L1FILE}.l1info)

    endforeach(filename)
endforeach(mission)

endif(OCSSWTestDir_FOUND)
