#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <netcdf.h>
#include <string.h>
#include <roxml.h>
#include <xmlUtils.h>

#include "l12_parms.h"
#include "filehandle.h"
#include "l1_mos_hdf.h"
#include "l1_octs_hdf.h"
#include "l1_osmi_hdf.h"
#include "l1_ocm_hdf.h"
#include "l1_ocm2_hdf.h"
#include "l1_meris_N1.h"
#include "l1_orca.h"
#include "l1_aviris.h"
#include "l1_prism.h"
#include "l1_viirs_h5.h"
#include "l1_hico_h5.h"
#include "l1_goci.h"
#include "l1_olci.h"

#define EOSMETALEN 32768

int32_t getFormat(filehandle *file) {
    int32_t fmt = -1;
    int32_t ocean_subset = -1;
    int32_t sd_id;
    char eosmeta[EOSMETALEN] = "";
    char tempstr[32] = "";
    int32_t chk_viirs(char *, filehandle *);
    int status;
    int32_t chk_hico(char *, filehandle *);
    int32_t chk_goci(char *, filehandle *);
    int32_t chk_oli(char *fname, filehandle *file);
    int32_t chk_oli_geo(char *fname);
    int32_t chk_aviris(char *fname, filehandle *file);
    int32_t chk_prism (char *fname, filehandle *file);
    int32_t chk_olci_xml(char *fname, filehandle *file);
    idDS ds_id;

    file->mode = READ;

    /* Does the file exist? */
    if (access(file->name, F_OK) || access(file->name, R_OK)) {
        printf("-E- %s: Input file '%s' does not exist or cannot open.\n",
        __FILE__, file->name);
        return (fmt);
    }


    // need to do HDF5 files before trying NetCDF
    if ((fmt = chk_viirs(file->name, file)) != -1) {
        if (want_verbose)
            printf("Input file %s is a VIIRS L1B HDF 5 file.\n", file->name);
        return (fmt);
    }

    if ((fmt = chk_hico(file->name, file)) != -1) {
        if (want_verbose)
            printf("Input file %s is a HICO L1B HDF 5 file.\n", file->name);
        return (fmt);
    }

    if ((fmt = chk_goci(file->name, file)) != -1) {
        if (want_verbose)
            printf("Input file %s is a GOCI L1B HDF 5 file.\n", file->name);
        return (fmt);
    }

    /* Is it netCDF */
    ds_id = openDS(file->name);
    if (ds_id.fid != FAIL) {
        if (ds_id.fftype == DS_NCDF) {
            char* titleStr = readAttrStr(ds_id, "title");

            if (titleStr) {
                lowcase(titleStr);
                if (strstr(titleStr, "viirs level-1a")){
                    fmt = FMT_VIIRSL1A;
                    file->format = fmt;
                    file->sensorID = VIIRS;
                    if (want_verbose)
                        printf("Input file %s is VIIRS L1A NetCDF4.\n", file->name);
                    free(titleStr);
                    endDS(ds_id);
                    return (fmt);
                }
                if (strstr(titleStr, "viirs m-band reflected solar band")){
                    fmt = FMT_VIIRSL1BNC;
                    file->format = fmt;
                    file->sensorID = VIIRS;
                    if (want_verbose)
                    printf("Input file %s is VIIRS L1B NetCDF4.\n", file->name);
                    free(titleStr);
                    endDS(ds_id);
                    return (fmt);
                }
                if (strstr(titleStr, "meris l1b")) {
                    fmt = FMT_MERISCC;
                    file->format = fmt;
                    file->sensorID = MERIS;
                    if (want_verbose)
                        printf("Input file %s is MERIS CC file.\n", file->name);
                    free(titleStr);
                    endDS(ds_id);
                    return (fmt);
                }
                if (strstr(titleStr, "oci level-1b")){
                    fmt = FMT_ORCA;
                    file->format = fmt;
                    file->sensorID = ORCA;
                    if (want_verbose)
                        printf("Input file %s is OCI L1B file.\n", file->name);
                    free(titleStr);
                    endDS(ds_id);
                    return (fmt);
                }
                if (strstr(titleStr, "olci level 1b product")){
//                    file->private_data = createPrivateData_olci(MAXOLCI_RADFILES);
                    fmt = FMT_OLCI;
                    file->format = fmt;
                    file->sensorID = OLCI;
                    if (want_verbose)
                        printf("Input file %s is OLCI L1B file.\n", file->name);
                    free(titleStr);
                    endDS(ds_id);
                    return (fmt);
                }

                char* platformStr = readAttrStr(ds_id, "platform");
                if (platformStr) {
                    char* instrumentStr = readAttrStr(ds_id, "instrument");
                    if (instrumentStr) {
                        char* processingLevelStr = readAttrStr(ds_id,
                                "processing_level");
                        if (processingLevelStr) {

                            if (!strcmp(processingLevelStr, "L1B")) {
                                fmt = FMT_L1BNCDF;
                            } else if (!strcmp(processingLevelStr, "L2")) {
                                fmt = FMT_L2NCDF;
                            } else if (!strcmp(processingLevelStr,
                                    "L3 Binned")) {
                                if (file->ocrvc_opt > 0)
                                    strcpy(instrumentStr, "OCRVC");
                                fmt = FMT_L3BIN;
                            } else if (!strcmp(processingLevelStr,
                                    "L3 Mapped")) {
                                fmt = FMT_L3MAP;
                            }
                            if (fmt != -1) {
                                file->format = fmt;
                                file->sensorID = instrumentPlatform2SensorID(
                                        instrumentStr, platformStr);
                                if (want_verbose)
                                    printf(
                                            "Input file %s is a NetCDF %s %s file.\n",
                                            file->name, instrumentStr,
                                            processingLevelStr);
                                free(processingLevelStr);
                                free(instrumentStr);
                                free(platformStr);
                                free(titleStr);
                                endDS(ds_id);
                                return (fmt);
                            }
                            free(processingLevelStr);
                        } // processingLevel found
                        free(instrumentStr);
                    } // instrument found
                    free(platformStr);
                } // platform found
                free(titleStr);
            } // title found
        } // is a NetCDF file
        endDS(ds_id);
    } // data set opened successfully

    /* Is it HDF? */
    sd_id = SDstart(file->name, DFACC_RDONLY);
    if (sd_id != FAIL) {

        /* File is HDF. Is it one of ours? */

        char title[255];
        char sensor[80];
        if (SDreadattr(sd_id, SDfindattr(sd_id, "Title"), (VOIDP) title) == 0) {
            if (strstr(title, "Level-3 Binned Data") != NULL) {
                if (SDreadattr(sd_id, SDfindattr(sd_id, "Sensor Name"), (VOIDP) sensor) == 0) {

                    if (file->ocrvc_opt > 0) strcpy(sensor, "OCRVC");

                    // kludge for VIIRS EDR L3
                    if (strcmp(sensor, "VIIRS") == 0) strcpy(sensor, "VIIRSN");

                    if ((file->sensorID = sensorName2SensorId(sensor)) != -1) {
                        fmt = FMT_L3BIN;
                        file->format = fmt;
                        if (want_verbose)
                            printf("Input file %s is %s.\n", file->name, title);
                    } else {
                        fprintf(stderr,
                                "-E- %s Line %d: Unknown sensor name in Level-3 file %s\n",
                                __FILE__, __LINE__, file->name);
                        return (fmt);
                    }
                } else {
                    fprintf(stderr,
                            "-E- %s Line %d: No sensor name attribute in Level-3 file %s\n",
                            __FILE__, __LINE__, file->name);
                    return (fmt);
                }
            } else if (strstr(title, "Level-2 Data") != NULL) {
                if (SDreadattr(sd_id, SDfindattr(sd_id, "Sensor Name"), (VOIDP) sensor) == 0) {
                    if ((file->sensorID = sensorName2SensorId(sensor)) != -1) {
                        fmt = FMT_L2HDF;
                        file->format = fmt;
                        if (want_verbose)
                            printf("Input file %s is %s.\n", file->name, title);
                    } else {
                        fprintf(stderr,
                                "-E- %s Line %d: Unknown sensor name in Level-2 file %s\n",
                                __FILE__, __LINE__, file->name);
                        return (fmt);
                    }
                } else {
                    fprintf(stderr,
                            "-E- %s Line %d: No sensor name attribute in Level-2 file %s\n",
                            __FILE__, __LINE__, file->name);
                    return (fmt);
                }
            } else if (strcmp(title, "SeaWiFS Level-1A Data") == 0) {
                fmt = FMT_SEAWIFSL1A;
                file->format = fmt;
                file->sensorID = SEAWIFS;
                if (SDreadattr(sd_id, SDfindattr(sd_id, "Data Type"), (VOIDP) tempstr) == 0) {
                    if (strcmp(tempstr, "GAC") == 0) {
                        file->subsensorID = SEAWIFS_GAC;
                        if (want_verbose)
                            printf("Input file %s is SeaWiFS Level-1A GAC.\n",
                                    file->name);
                    } else if (strcmp(tempstr, "LAC") == 0) {
                        file->subsensorID = SEAWIFS_LAC;
                        if (want_verbose)
                            printf("Input file %s is SeaWiFS Level-1A LAC.\n",
                                    file->name);
                    } else {
                        file->subsensorID = SEAWIFS_LAC;
                        if (want_verbose)
                            printf(
                                    "Input file %s is assumed to be SeaWiFS Level-1A LAC.\n",
                                    file->name);
                    }
                } else {
                    file->subsensorID = SEAWIFS_LAC;
                    if (want_verbose)
                        printf(
                                "Input file %s is assumed to be SeaWiFS Level-1A LAC.\n",
                                file->name);
                }
            } else if (strcmp(title, "OCTS Level-1A GAC Data") == 0) {
                fmt = FMT_OCTSL1A;
                file->format = fmt;
                file->sensorID = OCTS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "OSMI Level-1A Data") == 0) {
                fmt = FMT_OSMIL1A;
                file->format = fmt;
                file->sensorID = OSMI;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "CZCS Level-1A Data") == 0) {
                fmt = FMT_CZCSL1A;
                file->format = fmt;
                file->sensorID = CZCS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "OCM1 Level-1B (OBPG)") == 0) {
                fmt = FMT_OCML1B;
                file->format = fmt;
                file->sensorID = OCM1;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strncmp(title, "OCM Level-1B", 12) == 0) {
                fmt = FMT_OCML1BDB;
                file->format = fmt;
                file->sensorID = OCM1;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "Oceansat OCM2 Level-1B Data") == 0) {
                fmt = FMT_OCM2L1B;
                file->format = fmt;
                file->sensorID = OCM2;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);

                /* generic L1B format support */
            } else if (strcmp(title, "SeaWiFS Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = SEAWIFS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);

            } else if (strcmp(title, "MERIS Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = MERIS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "VIIRS Level-1B") == 0
                    || strcmp(title, "VIIRSN Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = VIIRS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "OCM2 Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = OCM2;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "OCTS Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = OCTS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "MOS Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = MOS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "OCTS Level-1B LAC Data") == 0) {
                fmt = FMT_OCTSL1B;
                file->format = fmt;
                file->sensorID = OCTS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "OSMI Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = OSMI;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strcmp(title, "HMODIST Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = HMODIST;
                file->subsensorID = MODIS_TERRA;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", file->name, title);
                    printf(
                            "\n\n *** WARNING: polarization can not be computed from generic format ***\n\n");
                }
            } else if (strcmp(title, "HMODISA Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = HMODISA;
                file->subsensorID = MODIS_AQUA;
                if (want_verbose) {
                    printf("Input file %s is %s.\n", file->name, title);
                    printf(
                            "\n\n *** WARNING: polarization can not be computed from generic format ***\n\n");
                }
            } else if (strcmp(title, "CZCS Level-1B") == 0) {
                fmt = FMT_L1HDF;
                file->format = fmt;
                file->sensorID = CZCS;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);

            } else if (strstr(title, "Level-1 cross-calibration pixels") != NULL) {
                if (SDreadattr(sd_id, SDfindattr(sd_id, "sensorID"),
                        (VOIDP) &(file->sensorID)) != 0) {
                    fprintf(stderr,
                            "-E- %s Line %d: Unrecognized sensor name, title %s in input HDF file %s\n",
                            __FILE__, __LINE__, title, file->name);
                    return (fmt);
                }
                fmt = FMT_L1XCAL;
                file->format = fmt;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else if (strstr(title, "AVHRR") != NULL) {
                fmt = FMT_CLASSAVHRR;
                file->format = fmt;
                file->sensorID = AVHRR;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            } else {
                fprintf(stderr,
                        "-E- %s Line %d: Unrecognized title %s in input HDF file %s\n",
                        __FILE__, __LINE__, title, file->name);
                return (fmt);
            }

        } else if (SDreadattr(sd_id, SDfindattr(sd_id, "title"), (VOIDP) title)
                == 0) {

            if (strstr(title, "AVHRR") != NULL) {
                fmt = FMT_CLASSAVHRR;
                file->format = fmt;
                file->sensorID = AVHRR;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name, title);
            }

        } else if (SDreadattr(sd_id, SDfindattr(sd_id, "satellite"),
                (VOIDP) title) == 0) {

            if (strstr(title, "oceansat-1") != NULL) {
                fmt = FMT_OCML1BDB;
                file->format = fmt;
                file->sensorID = OCM1;
                if (want_verbose)
                    printf("Input file %s is %s.\n", file->name,
                            "OCM1 DB file");
            }

            /* Is it HDF-EOS L1B format? */

        } else if (SDreadattr(sd_id, SDfindattr(sd_id, "ArchiveMetadata.0"),
                (VOIDP) eosmeta) == 0) {

            if (strstr(eosmeta, "MODIS/Terra Calibrated Radiances 5-Min L1B Swath 1km") != NULL) {

                fmt = FMT_HMODISL1B;
                file->format = fmt;
                file->sensorID = HMODIST;
//                file->subsensorID = MODIS_TERRA;

                if (want_verbose)
                    printf("Input file %s is MODIS Terra Level-1B HDF-EOS product.\n", file->name);
                return (fmt);
            } else if (strstr(eosmeta,
                    "MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 1km") != NULL) {

                fmt = FMT_HMODISL1B;
                file->format = fmt;
                file->sensorID = HMODISA;
//                file->subsensorID = MODIS_AQUA;

                if (want_verbose)
                    printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product.\n", file->name);
                return (fmt);
            } else if (strstr(eosmeta, "MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 250m") != NULL) {
                printf("Input file %s is MODIS Aqua Level-1B HDF-EOS product at 250m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=250\n",
                        file->name);
                return (fmt);
            } else if (strstr(eosmeta, "MODIS/Aqua Calibrated Radiances 5-Min L1B Swath 500m") != NULL) {
                printf( "Input file %s is MODIS Aqua Level-1B HDF-EOS product at 500m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=500\n",
                        file->name);
                return (fmt);
            } else if (strstr(eosmeta, "MODIS/Terra Calibrated Radiances 5-Min L1B Swath 250m") != NULL) {
                printf( "Input file %s is MODIS Terra Level-1B HDF-EOS product at 250m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=250\n",
                        file->name);
                return (fmt);
            } else if (strstr(eosmeta, "MODIS/Terra Calibrated Radiances 5-Min L1B Swath 500m") != NULL) {
                printf( "Input file %s is MODIS Terra Level-1B HDF-EOS product at 500m Resolution.\nThis file is no longer accepted as input.\nPlease use the 1KM (LAC) file and set resolution=500\n",
                        file->name);
                return (fmt);
            } else if (strstr(eosmeta, "MODIS/Aqua Geolocation Fields") != NULL) {
                fmt = FMT_MODISGEO;
                file->format = fmt;
                file->sensorID = HMODISA;
                file->subsensorID = MODIS_AQUA;
                if (want_verbose)
                    printf("Input file %s is MODIS Aqua Geolocation Fields.\n", file->name);
            } else if (strstr(eosmeta, "MODIS/Terra Geolocation Fields") != NULL) {
                fmt = FMT_MODISGEO;
                file->format = fmt;
                file->sensorID = HMODIST;
                file->subsensorID = MODIS_TERRA;
                if (want_verbose)
                    printf("Input file %s is MODIS Terra Geolocation Fields.\n", file->name);
            } else {
                fprintf(stderr,  "-E- %s Line %d: Unrecognized HDF-EOS file %s\n", __FILE__, __LINE__, file->name);
                return (fmt);
            }

            /* Is it MODIS IMAPP Direct Broadcast L1B */

        } else if (SDreadattr(sd_id,
                SDfindattr(sd_id, "ASSOCIATEDPLATFORMSHORTNAME"),
                (VOIDP) eosmeta) == 0) {

            if (strstr(eosmeta, "Aqua") != NULL) {
                fmt = FMT_MODISL1B;
                file->format = fmt;
                file->sensorID = HMODISA;
                file->subsensorID = MODIS_AQUA;
                if (want_verbose)
                    printf("Input file %s is MODIS Aqua Level-1B IMAPP product.\n", file->name);
            } else if (strstr(eosmeta, "Terra") != NULL) {
                fmt = FMT_MODISL1B;
                file->format = fmt;
                file->sensorID = HMODIST;
                file->subsensorID = MODIS_TERRA;
                if (want_verbose)
                    printf("Input file %s is MODIS Terra Level-1B IMAPP product.\n", file->name);
            } else {
                fprintf(stderr, "-E- %s Line %d: Unrecognized IMAPP format %s\n",  __FILE__, __LINE__, file->name);
                return (fmt);
            }

            /* Is it MOS L1B HDF standard product? */

        } else if (GetFileDesc(file->name) != NULL) {
            fmt = FMT_MOSL1B;
            file->format = fmt;
            file->sensorID = MOS;
            if (want_verbose)
                printf("Input file %s is MOS Level-1B standard product.\n", file->name);
        } else {
            fprintf(stderr, "-E- %s Line %d: Unrecognized input HDF file %s\n", __FILE__, __LINE__, file->name);
            return (fmt);
        }

        SDend(sd_id);

    } else {

        FILE *fp;
        char vfname[500];

        /* Open file */
        if ((fp = fopen(file->name, "r")) == NULL) {
            fprintf(stderr, "-E- : input file %s does not exist or is read protected.\n", file->name);
            return (fmt);
        }

        /* see if we have list of VIIRS files */
        if (fscanf(fp, "%499s", vfname) == 1) {
            if ((fmt = chk_viirs(vfname, file)) != -1) {
                if (file->format == FMT_VIIRSGEO) {
                    fprintf(stderr, "-E- : File list cannot be used for a VIIRS GEO file\n");
                    return (-1);
                }
                fclose(fp);
                return (fmt);
            }
        }

        fclose(fp);

        /* Is it MERIS? */
        {
            EPR_SProductId *product_id;

            epr_init_api(e_log_debug, NULL, NULL);
            product_id = epr_open_product(file->name);
            if (product_id != NULL) {
                if (product_id->id_string[8] == '1') { /* it is a level 1 file */
                    fmt = FMT_MERISL1B;
                    file->format = fmt;
                    file->sensorID = MERIS;
                    if (want_verbose)
                        printf("Input file %s is MERIS L1 file.\n", file->name);
                }
                epr_close_product(product_id);
                /*remember to close api (epr_close_api();)*/
                return (fmt);
            }
        }

        // check for OLI
        if ((fmt = chk_oli(file->name, file)) != -1) {
            if (want_verbose)
                printf("Input file %s is a Landsat 8 OLI L1B GEOTIFF file.\n", file->name);
            file->format = fmt;
            file->sensorID = OLI;
            return (fmt);
        }

        // check for AVIRIS
        if ((fmt = chk_aviris(file->name, file)) != -1) {
            if (want_verbose)
                printf("Input file %s is a AVIRIS file.\n", file->name);
            file->format = fmt;
            file->sensorID = AVIRIS;
            return (fmt);
        }
        // check for PRISM
        if ((fmt = chk_prism(file->name, file)) != -1) {
            if (want_verbose)
                printf("Input file %s is a PRISM file.\n", file->name);
            file->format = fmt;
            file->sensorID = PRISM;
            return (fmt);
        }
        // check for OLCI - in case they specified the xml file
        if (ds_id.fid == FAIL && (fmt = chk_olci_xml(file->name, file)) != -1) {
            if (want_verbose)
                printf("Input file %s is an OLCI file.\n", file->name);
            file->format = fmt;
            file->sensorID = OLCI;
            return (fmt);
        }
  }

    return (fmt);
}

int32_t chk_goci(char *fname, filehandle *file) {
    /* ------------------------------------------------------------------------
     chk_goci

     purpose: check a file to see if it is a GOCI file

     Returns -1 if not hdf5 and GOCI L1B or the format code

     Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char *            fname            I      file to check
     filehandle *      file             I      input file information

     Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     Wayne Robinson    16 Feb, 2013    Original development

     -----------------------------------------------------------------------*/
    int32_t fmt = -1;
    h5io_str h5fid, g_id;
    char scene_title[100];

    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

    /* Turn off error handling */
    H5Eset_auto( H5E_DEFAULT, NULL, NULL);

    if (h5io_openr(fname, 0, &h5fid) == 0) {
        if (h5io_set_grp(&h5fid, "HDFEOS/POINTS/Scene Header", &g_id) == 0) {
            if (h5io_attr_exist(&g_id, "Scene Title") == 0) {
                if (h5io_rd_attr(&g_id, "Scene Title", (void *) scene_title) == 0) {
                    if (strncmp(scene_title, "GOCI Level-1B Data", 18) == 0) {
                        fmt = FMT_GOCIL1B;
                        file->format = fmt;
                        file->sensorID = GOCI;
                    }
                }
            }
            h5io_close(&g_id);
        }
        h5io_close(&h5fid);
    }

    /* Restore previous error handler */
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data);

    return (fmt);
}

int32_t chk_hico(char *fname, filehandle *file) {
    /* ------------------------------------------------------------------------
     chk_hico

     purpose: check a file to see if it is a HICO L1B file

     Returns -1 if not hdf5 and HICO L1B or the format code

     Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char *            fname            I      file to check
     filehandle *      file             I      input file information

     Modification history:
     Programmer        Date            Description of change
     ----------        ----            ---------------------
     Sean Bailey       16 Feb, 2013    Original development

     ------------------------------------------------------------------------*/
    int32_t fmt = -1;
    char inst_name[100], proclvl[50];
    h5io_str h5fid, nameGid, levelGid;
    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

    /* Turn off error handling */
    H5Eset_auto( H5E_DEFAULT, NULL, NULL);

    if (h5io_openr(fname, 0, &h5fid) == 0) {
        if (h5io_set_grp(&h5fid,
                "metadata/FGDC/Identification_Information/Platform_and_Instrument_Identification",
                &nameGid) == 0) {
            if (h5io_rd_attr(&nameGid, "Instrument_Short_Name",
                    (void *) inst_name) == 0) {
                if (strcmp(inst_name, "hico") == 0) {
                    if (h5io_set_grp(&h5fid,
                            "metadata/FGDC/Identification_Information/Processing_Level",
                            &levelGid) == 0) {
                        if (h5io_rd_attr(&levelGid,
                                "Processing_Level_Identifier", (void *) proclvl)
                                == 0) {
                            if (strcmp(proclvl, "Level-1B") == 0) {
                                fmt = FMT_HICOL1B;
                                file->format = fmt;
                                file->sensorID = HICO;
                            } // if level = Level-1B
                        } // if Processing_Level_Identifier read
                        h5io_close(&levelGid);
                    } // set group ...Processing_Level...
                } // if name is hico
            } // if Instrument_Short_Name read
            h5io_close(&nameGid);
        } // set group for short name
        h5io_close(&h5fid);
    } // open file

    /* Restore previous error handler */
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data);

    return (fmt);
}

int32_t chk_viirs(char *fname, filehandle *file)
/* ------------------------------------------------------------------------
 chk_viirs

 purpose: check a file to see if it is a VIIRS file

 Returns -1 if not hdf5 and VIIRS or the format code

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 char *            fname            I      file to check
 filehandle *      file             I      input file information

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson, SAIC,  2 Aug, 2010   Original development

 ------------------------------------------------------------------------*/{
    int32_t fmt = -1;
    char inst_name[200], miss_name[200], geoloc_ref[200], **grp_obj_nm;
    h5io_str h5fid, g_id;
    int is_viirs = 0, *grp_obj_typ, n_obj;

    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);

    /* Turn off error handling */
    H5Eset_auto( H5E_DEFAULT, NULL, NULL);

    if (h5io_openr(fname, 0, &h5fid) == 0) {
        /*
         *  If the Instrument_Short_Name is in the top, make sure it is VIIRS
         */
        if (h5io_set_grp(&h5fid, "Data_Products/VIIRS-M1-SDR", &g_id) == 0) {
            if (h5io_attr_exist(&g_id, "Instrument_Short_Name") == 0) {
                if (h5io_rd_attr(&g_id, "Instrument_Short_Name",
                        (void *) inst_name) == 0) {
                    h5io_close(&g_id);
                    if (strcmp(inst_name, "VIIRS") == 0) is_viirs = 1;
                }
            }
        }
        //    				if (h5io_attr_exist(&h5fid, "Data_Products/VIIRS-M1-SDR/Instrument_Short_Name") == 0) {
        //    						}
        //            if (h5io_rd_attr(&h5fid, "Data_Products/VIIRS-M1-SDR/Instrument_Short_Name", (void *) inst_name) == 0){
        //            	if (strcmp(inst_name, "VIIRS") == 0) is_viirs = 1;
        //            }
        //        }
        if (is_viirs && (h5io_rd_attr(&h5fid, "Mission_Name", (void *) miss_name) != 0))
            is_viirs = 0;

        /*
         *  To accomodate the observed (in datasets) and possible (in CDFCB
         *  metadata V5 doc) mission names, look for either NPP or NPOESS
         *  at the start of the string
         */
        if (is_viirs && ((strncmp(miss_name, "NPP", 3) != 0) && (strncmp(miss_name, "NPOESS", 6) != 0)))
            is_viirs = 0;
        /* it is a VIIRS file - only need to know if its band or geo */
        /*
         *  set to the All_Data group and look for a geo or rad group
         */
        if (is_viirs && (h5io_set_grp(&h5fid, "All_Data", &g_id) != 0))
            is_viirs = 0;

        if (is_viirs && (h5io_grp_contents(&g_id, &n_obj, &grp_obj_nm, &grp_obj_typ) != 0))
            is_viirs = 0;

        if (is_viirs) {
            /*
             *  This is set up only for un-bundled bands / geo - if otherwise, note it
             */
            if (n_obj != 1) {
                printf("%s, %d: Warning: the file:\n%s\nIs bundled\n", __FILE__, __LINE__, fname);
                printf("\twith %d items\n", n_obj);
            }
            if ((strcmp(grp_obj_nm[0], "VIIRS-MOD-GEO_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-MOD-GEO-TC_All") == 0)) {
                fmt = FMT_VIIRSGEO;
                file->format = fmt;
                file->sensorID = VIIRS;
                if (want_verbose)
                    printf("Input file %s is VIIRS GEO(LOCATION) file.\n", fname);
            } else if ((strcmp(grp_obj_nm[0], "VIIRS-M1-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M2-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M3-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M4-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M5-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M6-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M7-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M8-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M9-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M10-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M11-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M12-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M13-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M14-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M15-SDR_All") == 0)
                    || (strcmp(grp_obj_nm[0], "VIIRS-M16-SDR_All") == 0)) {
                fmt = FMT_VIIRSL1B;
                file->format = fmt;
                file->sensorID = VIIRS;
                if (want_verbose)
                    printf("Input file %s is VIIRS Level-1B standard product.\n", fname);
            }
        }
        h5io_close(&h5fid);
    }

    /* Restore previous error handler */
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data);

    return (fmt);
}

int32_t chk_oli(char *fname, filehandle *file) {
    /* ------------------------------------------------------------------------
     chk_oli

     purpose: check a file to see if it is an OLI Landsat8 file

     Returns -1 if not OLI L1B or the format code

     Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char *            fname            I      file to check
     filehandle *      file             I      input file information

     -----------------------------------------------------------------------*/
    const int lineSize = 500;
    int i;
    int32_t fmt = -1;
    FILE *fp;
    char line[lineSize + 1];
    char* result;

    /* Open file */
    if ((fp = fopen(file->name, "r")) == NULL) {
        return (fmt);
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return fmt;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "GROUP = L1_METADATA_FILE"
    if (strstr(line, "L1_METADATA_FILE") == NULL) {
        fclose(fp);
        return fmt;
    }

    // within 20 lines look for:
    //   SPACECRAFT_ID = "LANDSAT_8"
    //   SENSOR_ID = "OLI_TIRS"
    int foundSpacecraft = 0;
    int foundSensor = 0;
    for (i = 0; i < 20; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return fmt;
        }
        trimBlanks(line);
        if (strstr(line, "SPACECRAFT_ID")) {
            if (strstr(line, "LANDSAT_8")) foundSpacecraft = 1;
        } else if (strstr(line, "SENSOR_ID"))
            if (strstr(line, "OLI_TIRS")) foundSensor = 1;

        if (foundSpacecraft && foundSensor) {
            fmt = FMT_OLIL1B;
            break;
        }
    }

    fclose(fp);
    return fmt;
}

int32_t chk_oli_geo(char *fname) {
    /* ------------------------------------------------------------------------
     chk_oli_geo

     purpose: check a file to see if it is an OLI Landsat8 GEO file

     Returns -1 if not OLI L1B or the format code

     Parameters: (in calling order)
     Type              Name            I/O     Description
     ----              ----            ---     -----------
     char *            fname            I      file to check
     filehandle *      file             I      input file information

     -----------------------------------------------------------------------*/
    const int lineSize = 500;
    int i;
    int32_t fmt = -1;
    FILE *fp;
    char line[lineSize + 1];
    char* result;

    /* Open file */
    if ((fp = fopen(fname, "r")) == NULL) {
        return (fmt);
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return fmt;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "GROUP = FILE_HEADER"
    if (strstr(line, "FILE_HEADER") == NULL) {
        fclose(fp);
        return fmt;
    }

    // within 20 lines look for:
    //   SATELLITE = "LANDSAT_8"
    //   BAND_LIST = ...
    int foundSpacecraft = 0;
    int foundBand = 0;
    for (i = 0; i < 20; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return fmt;
        }
        trimBlanks(line);
        if (strstr(line, "SATELLITE")) {
            if (strstr(line, "LANDSAT_8")) foundSpacecraft = 1;
        } else if (strstr(line, "BAND_LIST")) foundBand = 1;

        if (foundSpacecraft && foundBand) {
            fmt = FMT_OLIL1B;
            break;
        }
    }

    fclose(fp);
    return fmt;
}
/**
 *   chk_prism
 *   purpose: check a file to see if it is PRISM hdr file
 *   Returns -1 if not PRISM or the format code
 *
 * @param char *fname      - name of file to check
 * @param filehandle *file - input file information
 */

int32_t chk_prism(char *fname, filehandle *file) {
    const int lineSize = 500;
    int i;
    int32_t fmt = -1;
    FILE *fp;
    char line[lineSize + 1];
    char* result;

    if (strstr(file->name,"prm") == NULL) {
         return fmt;
    }

    /* Open file */
    if ((fp = fopen(file->name, "r")) == NULL) {
        return (fmt);
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return fmt;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "ENVI"
    if (strstr(line, "ENVI") == NULL) {
        fclose(fp);
        return fmt;
    }

    fmt = FMT_PRISM;

    fclose(fp);
    return fmt;
}

/**
 *   chk_aviris
 *   purpose: check a file to see if it is an AVIRIS hdr file
 *   Returns -1 if not AVIRIS or the format code
 *
 * @param char *fname      - name of file to check
 * @param filehandle *file - input file information
 */

int32_t chk_aviris(char *fname, filehandle *file) {
    const int lineSize = 500;
    int i;
    int32_t fmt = -1;
    FILE *fp;
    char line[lineSize + 1];
    char* result;

    /* Open file */
    if ((fp = fopen(file->name, "r")) == NULL) {
        return (fmt);
    }

    // skip blank lines
    do {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return fmt;
        }
        trimBlanks(line);
    } while (strlen(line) == 0);

    // first line needs to be "ENVI"
    if (strstr(line, "ENVI") == NULL) {
        fclose(fp);
        return fmt;
    }

    // within 20 lines look for:
    //   AVIRIS orthocorrected
    int foundAviris = 0;
    int foundFormat = 0;
    int foundOrtho = 0;
    for (i = 0; i < 20; i++) {
        result = fgets(line, lineSize, fp);
        if (result == NULL) {
            fclose(fp);
            return fmt;
        }
        trimBlanks(line);
        if (strstr(line, "AVIRIS"))         foundAviris = 1;
        if (strstr(line, "orthocorrected")) foundOrtho = 1;
        if (strstr(line, "interleave"))     foundFormat = 1;

        if (foundAviris && foundFormat && foundOrtho) {
            fmt = FMT_AVIRIS;
            break;
        }
    }

    fclose(fp);
    return fmt;
}

/**
 *   chk_olci
 *   purpose: check a file to see if it is an OLCI xml file
 *   If it is, populate the data pointer with filenames
 *   associated with Lt's, instrument and geometry data
 *   Returns -1 if not OLCI or the format code
 *
 * @param char *fname      - name of file to check
 * @param filehandle *file - input file information
 *
 * Author: R. Healy (SAIC) 7/29/2015
 * W. Robinson, SAIC, 18 Feb 2016 pad out the olci file name storage by 1 
 *   to accomodate trailing null, allow xml files with alternate product 
 *   name to be accepted as olci files
 * s
 */

int32_t chk_olci_xml(char *fname, filehandle *file) {
    const int lineSize = 500;
    node_t* rootNode;
//    node_t* xfduNodes;
    node_t* xfduNode;
    node_t* packageNode;
    node_t* dataNode;
    node_t* dataObjNode;
    node_t* byteNode;
    node_t* fileNode;
    node_t* xfduc;
    node_t* radNode;
    node_t* radFNode;
    olci_t *data = file->private_data = createPrivateData_olci(MAXOLCI_RADFILES);

    char productName[TITLELEN],dataName[TITLELEN];
    int i, fln, nbands=0;
    int32_t fmt = -1;
    FILE *fp;
    char line[lineSize + 1];
    char* result;

    rootNode = roxml_load_doc(fname);
    if(rootNode == NULL) {
        roxml_close(rootNode);
       return fmt;
    }

    int all_nodes_1 = roxml_get_nodes_nb(rootNode, ROXML_ELM_NODE | ROXML_CMT_NODE | ROXML_PI_NODE | ROXML_TXT_NODE | ROXML_ATTR_NODE);
    int all_nodes_2 = roxml_get_nodes_nb(rootNode, ROXML_ALL_NODES);
    if(all_nodes_1 == all_nodes_2) {
        printf("%d Nodes are contained in root\n", all_nodes_1);
     }
//    printf("Number of children=%d\n",roxml_get_chld_nb  (rootNode));
//    xfduNodes = roxml_get_nodes(rootNode, ROXML_ALL_NODES, NULL,0);
    xfduNode = roxml_get_chld(rootNode, "XFDU",0);
    packageNode = roxml_get_chld(xfduNode, "informationPackageMap",0);
    dataObjNode = roxml_get_chld(xfduNode, "dataObjectSection",0);
    xfduc = roxml_get_chld(packageNode, NULL,0);
//    printf("Number of xfdu children=%d\n",roxml_get_chld_nb  (xfduNode));
//    printf("Number of package children=%d\n",roxml_get_chld_nb  (packageNode));
//    printf("Number of xfduc children=%d\n",roxml_get_chld_nb  (xfduc));
     xmlGetAttributeStr(xfduc, "textInfo", productName, TITLELEN);
    if ( strstr(productName,"OLCI" ) &&
         strstr(productName,"SENTINEL-3" )  &&
         strstr(productName,"Level 1" ) ) {
                printf("%s\n",productName);
                fmt = FMT_OLCI;
    } else {
        printf("OLCI productName does not match OLCI && SENTINEL-3 && Level 1:" );
        if (strstr(productName,"OLCI Level 2")) {
            printf("\n%s\nThis appears to be an OLCI Level 2 product.  Sorry.\n",productName);
            exit(EXIT_FAILURE);
        }
        roxml_close(rootNode);
        return fmt;
    }
    dataNode = roxml_get_chld(dataObjNode, NULL, 0);
    int foundGeoCoord = 0;
    int foundtieGeoCoord = 0;
    int foundtieGeometries = 0;
    int foundinstrumentDat = 0;
    int foundtimeCoord = 0;
    int foundTieMeteo = 0;

    while(dataNode) {
        byteNode = roxml_get_chld(dataNode, NULL, 0);
        fileNode = roxml_get_chld(byteNode, NULL, 0);
        xmlGetAttributeStr(fileNode, "href", productName, TITLELEN);
        xmlGetAttributeStr(dataNode, "ID", dataName, TITLELEN);
        //printf("Dataname=%s\n",dataName);
        if (strstr(dataName,"radianceData")) {
            if (nbands > MAXOLCI_RADFILES) {
                printf("%s, %d - E - Maximum number of radiance files (%d) reached\n",
                __FILE__, __LINE__,MAXOLCI_RADFILES);
                exit(EXIT_FAILURE);
            }
            fln = strlen(productName);
            if ((data->olci_radfiles[nbands] = (char *) malloc((fln+1)*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate radiance filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcpy(data->olci_radfiles[nbands],productName);
            //printf("dataObjectID=%s len=%d\n",data->olci_radfiles[nbands],fln);
            nbands++;
        } else if (strcmp(dataName,"geoCoordinatesData") == 0) {
            fln = strlen(productName);
            if ((data->geoCoordinatesFile = (char *) malloc((fln+1)*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate geo coordinates filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcpy(data->geoCoordinatesFile,productName);
            //printf("Geo Coordinate file=%s\n",data->geoCoordinatesFile);
            foundGeoCoord = 1;
        } else if (strcmp(dataName,"tieGeoCoordinatesData") == 0) {
            fln = strlen(productName);
            if ((data->tieGeoCoordinatesFile = (char *) malloc((fln+1)*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate tie geo coordinates filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcpy(data->tieGeoCoordinatesFile,productName);
            //printf("Geo Coordinate file=%s\n",data->tieGeoCoordinatesFile);
            foundtieGeoCoord = 1;
        }else if (strstr(dataName,"tieGeometriesData")) {
            fln = strlen(productName);
            if ((data->tieGeometriesFile = (char *) malloc((fln+1)*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate tie geometries filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcpy(data->tieGeometriesFile,productName);
            //printf("Tie Geometries file=%s\n",data->tieGeometriesFile);
            foundtieGeometries = 1;
        } else if (strstr(dataName,"instrumentDataData")) {
            fln = strlen(productName);
            if ((data->instrumentFile = (char *) malloc((fln+1)*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate instrument filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcpy(data->instrumentFile,productName);

            //printf("Instrument File file=%s\n",data->instrumentFile);
            foundinstrumentDat = 1;
        } else if (strstr(dataName,"timeCoordinatesData")) {
            fln = strlen(productName);
            if ((data->time_coordinatesFile = (char *) malloc((fln+1)*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate time_coordinates filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcpy(data->time_coordinatesFile,productName);
            //printf("Time Coordinates file=%s\n",data->time_coordinatesFile);
            foundtimeCoord = 1;
        } else if (strstr(dataName,"tieMeteoData")) {
            fln = strlen(productName);
            if ((data->tieMeteoFile = (char *) malloc((fln+1)*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate tieMeteoFile filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            strcpy(data->tieMeteoFile,productName);
            //printf("Tie Meteo file=%s\n",data->tieMeteoFile );
            foundTieMeteo = 1;
        }
        //        printf("dataObjectID=%s len=%d\n",productName,fln);
        dataNode = roxml_get_next_sibling(dataNode);
    }
    data->numRadFiles = nbands;
    if (foundtimeCoord && foundinstrumentDat && foundtieGeometries && foundtieGeoCoord && foundtieGeoCoord && foundGeoCoord  && foundTieMeteo && nbands >= data->numRadFiles) {
        printf("OLCI:  Found all meta data information\n");
    }else{
        printf("%s, %d - E - missing meta data information in file %s \n",
                __FILE__, __LINE__,fname);
        exit(EXIT_FAILURE);
    }
    roxml_close(rootNode);
    return fmt;
}

