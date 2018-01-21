#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <palette.h>
#include <genutils.h>
#include <clo.h>
#include <dfutils.h>
#include <sensorInfo.h>

#include "l3mapgen.h"

// should the par file processing descend into other par files
static int enableFileDescending = 1;

/** CLO callback function for the "par" option.  Loads the parameter file
 into the list that is stored in option->cb_data */
void par_option_cb(struct clo_option_t *option) {
    if (enableFileDescending)
        clo_readFile((clo_optionList_t*) option->cb_data, option->valStr);
}

/** add all of the accepted command line options to list */
int l3mapgen_init_options(clo_optionList_t* list, const char* softwareVersion) {
    char tmpStr[2048];
    clo_option_t* option;
    int i;

    clo_setVersion2("l3mapgen", softwareVersion);

    sprintf(tmpStr, "Usage: l3mapgen argument-list\n\n");

    strcat(tmpStr, "  This program takes a product (or products if netCDF output) from an L3 bin\n");
    strcat(tmpStr, "  or SMI file, reprojects the data using Proj.4 and writes a mapped file in\n");
    strcat(tmpStr, "  the requested output format.\n\n");

    strcat(tmpStr, "  Return values\n");
    strcat(tmpStr, "    0 = All Good\n");
    strcat(tmpStr, "    1 = Error\n");
    strcat(tmpStr, "    110 = north/south or east/west extents identical\n\n");

    strcat(tmpStr, "  The argument-list is a set of keyword=value pairs. The arguments can\n");
    strcat(tmpStr, "  be specified on the commandline, or put into a parameter file, or the\n");
    strcat(tmpStr, "  two methods can be used together, with commandline over-riding.\n\n");
    strcat(tmpStr, "The list of valid keywords follows:\n");
    clo_setHelpStr(tmpStr);

    // add the par option and add the callback function
    option = clo_addOption(list, "par", CLO_TYPE_STRING, NULL, "input parameter file");
    option->cb_data = (void*) list;
    option->cb = par_option_cb;

    clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, "input L3 bin file name");
    clo_addOption(list, "ofile", CLO_TYPE_OFILE, "output", "output filename");
    strcpy(tmpStr, "output file format\n");
    strcat(tmpStr, "        netcdf4: netCDF4 file, can contain more than one product\n");
    strcat(tmpStr, "        hdf4:    HDF4 file (old SMI format)\n");
    strcat(tmpStr, "        png:     PNG image file\n");
    strcat(tmpStr, "        ppm:     PPM image file\n");
    strcat(tmpStr, "        tiff:    TIFF file with georeference tags");
    clo_addOption(list, "oformat", CLO_TYPE_STRING, "netcdf4", tmpStr);

    clo_addOption(list, "ofile2", CLO_TYPE_OFILE, NULL, "second output filename");
    strcpy(tmpStr, "second output file format\n");
    strcat(tmpStr, "        netcdf4: netCDF4 file, can contain more than one product\n");
    strcat(tmpStr, "        hdf4:    HDF4 file (old SMI format)\n");
    strcat(tmpStr, "        png:     PNG image file\n");
    strcat(tmpStr, "        ppm:     PPM image file\n");
    strcat(tmpStr, "        tiff:    TIFF file with georeference tags");
    clo_addOption(list, "oformat2", CLO_TYPE_STRING, "png", tmpStr);

    clo_addOption(list, "deflate", CLO_TYPE_INT, "4", "deflation level");

    strcpy(tmpStr, "comma separated list of products.\n");
    strcat(tmpStr, "        Each product can have an optional colon and modifier appended.\n");
    strcat(tmpStr, "        For example  product=chlor_a,chlor_a:stdev,Kd_490:nobs\n");
    strcat(tmpStr, "        The avaliable modifiers are:\n");
    strcat(tmpStr, "            avg       average value (default)\n");
    strcat(tmpStr, "            stdev     standard deviation\n");
    strcat(tmpStr, "            var       variance\n");
    strcat(tmpStr, "            nobs      number of observations in the bin\n");
    strcat(tmpStr, "            nscenes   number of contributing scenes\n");
    strcat(tmpStr, "            obs_time  average observation time (TAI93)\n");
    strcat(tmpStr, "            bin_num   bin ID number");
    clo_addOption(list, "product", CLO_TYPE_STRING, NULL, tmpStr);

    strcpy(tmpStr, "size of the output pixel in meters or\n");
    strcat(tmpStr, "        SMI dimensions\n");
    strcat(tmpStr, "        90km: 432 x 216 image for full globe\n");
    strcat(tmpStr, "        36km: 1080 x 540\n");
    strcat(tmpStr, "        18km: 2160 x 1080\n");
    strcat(tmpStr, "         9km: 4320 x 2160\n");
    strcat(tmpStr, "         4km: 8640 x 4320\n");
    strcat(tmpStr, "         2km: 17280 x 8640\n");
    strcat(tmpStr, "         1km: 34560 x 17280\n");
    strcat(tmpStr, "         hkm: 69120 x 34560\n");
    strcat(tmpStr, "         qkm: 138240 x 69120\n");
    strcat(tmpStr, "         smi: 4096 x 2048\n");
    strcat(tmpStr, "        smi4: 8192 x 4096\n");
    strcat(tmpStr, "        land: 8640 x 4320\n");
    strcat(tmpStr, "         #.#:  width of a pixel in meters\n");
    strcat(tmpStr, "       #.#km:  width of a pixel in kilometers\n");
    strcat(tmpStr, "      #.#deg:  width of a pixel in degrees");
    clo_addOption(list, "resolution", CLO_TYPE_STRING, "9km", tmpStr);

    strcpy(tmpStr, "proj.4 projection string or one\n");
    strcat(tmpStr, "        of the following predefined projections:\n");
    strcat(tmpStr, "        smi:       Standard Mapped image, cylindrical projection, uses\n");
    strcat(tmpStr, "                   central_meridian.  n,s,e,w default to whole globe.\n");
    strcat(tmpStr, "                   projection=\"+proj=eqc +lat_0=<central_meridian>\"\n");
    strcat(tmpStr, "        platecarree: Plate Carree image, cylindrical projection, uses\n");
    strcat(tmpStr, "                   central_meridian\n");
    strcat(tmpStr, "                   projection=\"+proj=eqc +lat_0=<central_meridian>\"\n");
    strcat(tmpStr, "        mollweide: Mollweide projection \n");
    strcat(tmpStr, "                   projection=\"+proj=moll +lat_0=<central_meridian>\"\n");
    strcat(tmpStr, "        lambert:   Lambert conformal conic projection \n");
    strcat(tmpStr, "                   projection=\"+proj=lcc +lat_0=<central_meridian>\"\n");
    strcat(tmpStr, "        albersconic: Albers equal-area conic projection \n");
    strcat(tmpStr, "                   projection=\"+proj=aea +lat_0=<central_meridian>\"\n");
    strcat(tmpStr, "        mercator:  Mercator cylindrical map projection \n");
    strcat(tmpStr, "                   projection=\"+proj=merc +lat_0=<central_meridian>\"\n");
    strcat(tmpStr, "        ease2:     Ease Grid 2 projection \n");
    strcat(tmpStr, "                   projection=\"+proj=cea +lon_0=0 +lat_ts=30 +ellps=WGS84\n");
    strcat(tmpStr, "                         +datum=WGS84 +units=m +lat_0=<central_meridian>\"\n");
    strcat(tmpStr, "        raw:       raw dump of the bin file");
    clo_addOption(list, "projection", CLO_TYPE_STRING, "platecarree", tmpStr);
    clo_addOption(list, "central_meridian", CLO_TYPE_FLOAT, "-999", "central meridian for projection in\n        deg east.  Only used for smi, mollweide and raw projection");

    strcpy(tmpStr, "interpolation method:\n");
    strcat(tmpStr, "        nearest: Nearest Neighbor\n");
//    strcat(tmpStr, "        linear:  bi-linear interpolation using the corner coords of the\n                  output pixel\n");
    strcat(tmpStr, "        bin:     bin all of the pixels that intersect the area of the\n                  output pixel\n");
    strcat(tmpStr, "        area:    bin weighted by area of all the pixels that intersect\n                  the area of the output pixel");
    clo_addOption(list, "interp", CLO_TYPE_STRING, "nearest", tmpStr);
        
    clo_addOption(list, "north", CLO_TYPE_FLOAT, "-999", "Northern most Latitude (-999=file north)");
    clo_addOption(list, "south", CLO_TYPE_FLOAT, "-999", "Southern most Latitude (-999=file south)");
    clo_addOption(list, "east", CLO_TYPE_FLOAT, "-999", "Eastern most Longitude (-999=file east)");
    clo_addOption(list, "west", CLO_TYPE_FLOAT, "-999", "Western most Longitude (-999=file west)");

    strcpy(tmpStr, "apply color A palette:\n");
    strcat(tmpStr, "        yes: color image\n");
    strcat(tmpStr, "         no: grayscale image");
    clo_addOption(list, "apply_pal", CLO_TYPE_BOOL, "yes", tmpStr);
    clo_addOption(list, "palfile", CLO_TYPE_IFILE, NULL, "palette filename.  (default = means the palette\n        is chosen using the product.xml file");
    clo_addOption(list, "palette_dir", CLO_TYPE_IFILE, "$OCDATAROOT/common/palette", "directory\n        containing the palette files");

    clo_addOption(list, "datamin", CLO_TYPE_FLOAT, NULL, "minimum value for scaling (default from product.xml)");
    clo_addOption(list, "datamax", CLO_TYPE_FLOAT, NULL, "maximum value for scaling (default from product.xml)");

    strcpy(tmpStr, "data scaling type (default from product.xml)\n");
    strcat(tmpStr, "        linear:  linear scaling\n");
    strcat(tmpStr, "        log:     logarithmic scaling\n");
    strcat(tmpStr, "        arctan:  arc tangent scaling");
    clo_addOption(list, "scale_type", CLO_TYPE_STRING, NULL, tmpStr);
    clo_addOption(list, "quiet", CLO_TYPE_BOOL, "false", "stop the status printing");
    clo_addOption(list, "pversion", CLO_TYPE_STRING, "Unspecified", "processing version string");
    clo_addOption(list, "use_quality", CLO_TYPE_BOOL, "yes", "should we do quality factor processing");
    clo_addOption(list, "use_rgb", CLO_TYPE_BOOL, "no", "should we use product_rgb to make a\n        psudo-true color image");
    clo_addOption(list, "product_rgb", CLO_TYPE_STRING, "rhos_670,rhos_555,rhos_412", "3 products to use\n        for RGB.  Default is sensor specific");
    clo_addOption(list, "fudge", CLO_TYPE_FLOAT, "1.0", "fudge factor used to modify size of L3 pixels");
    clo_addOption(list, "threshold", CLO_TYPE_FLOAT, "0", "minimum percentage of filled pixels before\n        an image is generated");

    return 0;
}

const char* getSensorDirectory(const char* fileName) {
    idDS dsId;
    int sensorId = -1;
    
    dsId = openDS(fileName);
    if (dsId.fid == FAIL) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.\n",
        __FILE__, fileName);
        exit(EXIT_FAILURE);
    }

    if (dsId.fftype == DS_NCDF) {
        char* instrumentStr = readAttrStr(dsId, "instrument");
        if (instrumentStr) {
            char* platformStr = readAttrStr(dsId, "platform");
            if (platformStr) {
                sensorId = instrumentPlatform2SensorID(instrumentStr, platformStr);
            }
        }
    } else {
        char* sensorNameStr = readAttrStr(dsId, "Sensor Name");
        if (sensorNameStr) {
            sensorId = sensorName2SensorId(sensorNameStr);
        }
    }

    if(sensorId == -1) {
        printf("-E- %s: Input file '%s' does not have a valid sensor ID.\n",
        __FILE__, fileName);
        exit(EXIT_FAILURE);
    }
           
    return sensorDir[sensorId];
}

/* 
   Read the command line option and all of the default parameter files.

   This is the order for loading the options:
    - load the main program defaults file
    - load the command line (including specified par files)
    - re-load the command line disabling file descending so command
       line arguments will over ride

 */
int l3mapgen_read_options(clo_optionList_t* list, int argc, char* argv[]) {
    char *dataRoot;
    char tmpStr[FILENAME_MAX];
    int i;

    assert(list);

    if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
        fprintf(stderr, "-E- OCDATAROOT environment variable is not defined.\n");
        return (-1);
    }

    // disable the dump option until we have read all of the files
    clo_setEnableDumpOptions(0);

    // make sure file descending is turned on
    enableFileDescending = 1;

    // load program defaults
    sprintf(tmpStr, "%s/common/l3mapgen_defaults.par", dataRoot);
    clo_readFile(list, tmpStr);

    // read all arguments including descending par files
    clo_readArgs(list, argc, argv);

    // if a file was descended, make sure the command args over-rides
    enableFileDescending = 0;
    clo_readArgs(list, argc, argv);
   
    // get sensor directory
    const char* sensorDir = getSensorDirectory(clo_getString(optionList, "ifile"));

    // load the sensor specific defaults file
    sprintf(tmpStr, "%s/%s/l3mapgen_defaults.par", dataRoot, sensorDir);
    enableFileDescending = 1;
    clo_readFile(list, tmpStr);
    
    // enable the dump option
    clo_setEnableDumpOptions(1);

    enableFileDescending = 0;
    clo_readArgs(list, argc, argv);
    enableFileDescending = 1;

    return 0;
}
