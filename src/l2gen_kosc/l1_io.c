/* =========================================================== */
/* Module l1_io.c                                              */
/*                                                             */
/* Functions to open, close, read, and write  a level-1b file, */
/* with the format determined by the file handle content.      */
/*                                                             */
/* Written By:                                                 */
/*                                                             */
/*     B. A. Franz                                             */
/*     SAIC General Sciences Corp.                             */
/*     NASA/SIMBIOS Project                                    */
/*     April 1998                                              */
/*                                                             */
/* Modifications By:                                           */
/*    J. Gales                                                 */
/*    Futuretech                                               */
/*    NASA/SIMBIOS Project                                     */
/*    10/00                                                    */
/*                                                             */
/*    Add support for OCTSL1A                                  */
/*    W. Robinson, SAIC, 10 Dec 2004  add CZCS, VIIRS support  */
/* =========================================================== */

#include <stdio.h>
#include "l12_proto.h"
#include "l1_mos_hdf.h"
#include "l1_hdf_generic_read.h"
#include "l1a_seawifs.h"
#include "l1_octs_hdf.h"
#include "l1_osmi_hdf.h"
#include "l1_hmodis_hdf.h"
#include "l1_czcs_hdf.h"
#include "l1_xcal_hdf.h"
#include "l1_aci_hdf.h"
#include "l1_ocm_hdf.h"
#include "l1_ocm2_hdf.h"
#include "l1_ocmdb_hdf.h"
#include "l1_meris_N1.h"
#include "l1_meris_CC.h"
#include "l1_viirs_h5.h"
#include "l1b_viirs_nc.h"
#include "l1_hico_h5.h"
#include "l1_goci.h"
#include "l1_oli.h"
#include "l1_viirs_nc.h"
#include "l1_orca.h"
#include "l1_aviris.h"
#include "l1_prism.h"
#include "l1_olci.h"
#include "l1_nc_generic_read.h"

/* ---------------------------------------------------------------- */
/* Close the level 1 file associated with the input file handle.    */
/* ---------------------------------------------------------------- */
void closel1(filehandle *l1file)
{
    switch (l1file->format) {
        case FMT_L1HDF: 
            if (l1file->mode == READ)
                closel1_hdf_g(l1file);
            else
                closel1_generic(l1file);
            break;
        case FMT_L1BNCDF: 
            closel1_generic(l1file);
            break;
        case FMT_MOSL1B:
            closel1_mos_hdf(l1file);
            break;
        case FMT_SEAWIFSL1A:
            closel1a_seawifs(l1file);
            break;
        case FMT_OCTSL1B:
            closel1_octs_hdf(l1file);
            break;
        case FMT_OCTSL1A:
            closel1_octs_hdf(l1file);
            break;
        case FMT_OSMIL1A:
            closel1a_osmi(l1file);
            break;
        /* case FMT_MODISL1B: */
        /*     closel1_modis_hdf(l1file); */
        /*     break; */
        case FMT_L1XCAL:
            closel1_xcal_hdf(l1file);
            break;
        case FMT_CZCSL1A:
            closel1_czcs(l1file);
            break;
        case FMT_HMODISL1B:
            closel1_hmodis_hdf();
            break;
        case FMT_CLASSAVHRR:
             closel1_aci_hdf(l1file);
             break;
        case FMT_OCML1B:
             closel1_ocm_hdf(l1file);
             break;
        case FMT_OCM2L1B:
             closel1_ocm2_hdf(l1file);
             break;
        case FMT_OCML1BDB:
             closel1_ocmdb_hdf(l1file);
             break;
        case FMT_MERISL1B:
             closel1_meris_N1(l1file);
             break;
        case FMT_MERISCC:
             closel1_meris_CC(l1file);
             break;
        case FMT_VIIRSL1B:
             closel1_viirs_h5(l1file);
             break;
        case FMT_VIIRSL1BNC:
             closel1b_viirs_nc();
             break;
        case FMT_VIIRSL1A:
             closel1_viirs_nc(l1file);
             break;
        case FMT_HICOL1B:
             closel1_hico_h5(l1file);
             break;
        case FMT_GOCIL1B:
             closel1_goci(l1file);
             break;
        case FMT_OLIL1B:
             closel1_oli(l1file);
             break;
        case FMT_ORCA:
             closel1_orca(l1file);
             break;
        case FMT_AVIRIS:
             closel1_aviris(l1file);
             break;
        case FMT_PRISM:
             closel1_prism(l1file);
             break;
        case FMT_OLCI:
             closel1_olci(l1file);
             break;
        default:
            fprintf(stderr,
            "-E- %s Line %d: l1close - Unknown L1 file format specifier: %d\n",
            __FILE__,__LINE__,l1file->format);
            break;
      };

    return;
}

/* ---------------------------------------------------------------- */
/* Open a level 1 file for reading, load file handle with metadata  */
/* from file header. At a minimum, the name and format fields of    */
/* file handle must be loaded before this function is called.       */
/* ---------------------------------------------------------------- */
int openl1(filehandle *l1file, initstr *initrec)
{
    int status = 1;
    int32_t *wave;

    /* Init filehandle stats */
    l1file->percent_cloud  = -1;
    l1file->percent_land   = -1;
    l1file->percent_water  = -1;


    /* Get number of bands and band indexing from sensor table */
    l1file->nbands = rdsensorinfo(l1file->sensorID,l1file->input->evalmask,NULL,NULL);
    if (l1file->nbands < 0) {
        printf("-E- %s line %d: Error reading sensor table\n",
               __FILE__,__LINE__);
        return(status);
    }
    rdsensorinfo(l1file->sensorID,l1file->input->evalmask,"Bindx",(void **) &l1file->bindx);
    l1file->nbandsir = rdsensorinfo(l1file->sensorID,l1file->input->evalmask,"NbandsIR",NULL);

    /* set wavelength index */
    rdsensorinfo(l1file->sensorID,l1file->input->evalmask,"Lambda",(void **) &wave);
    bindex_set(wave,l1file->nbands+l1file->nbandsir,BANDW);

    /* Open L1 file for reading or writing, as requested */
    if (l1file->mode == READ) {

        switch (l1file->format) {
            case FMT_MOSL1B:
                status = openl1_read_mos_hdf(l1file);
                break;
            case FMT_SEAWIFSL1A:
                status = openl1a_seawifs(l1file);
                break;
            case FMT_L1HDF:
                status = openl1_read_hdf_g(l1file);
                break;
            case FMT_L1BNCDF:
                status = openl1_nc_generic(l1file);
                break;
            case FMT_OCTSL1B:
                status = openl1_read_octs_hdf(l1file);
                break;
            case FMT_OCTSL1A:
                status = openl1_read_octs_hdf(l1file);
                break;
            case FMT_OSMIL1A:
                status = openl1a_osmi(l1file);
                break;
            case FMT_MODISL1B:
                /* status = openl1_modis_hdf(l1file); */
                printf("This L1B file contains only the ocean band subset; "
                       "processing is no longer supported.\n");
                break;
            case FMT_L1XCAL:
                status = openl1_xcal_hdf(l1file);
                break;
            case FMT_CZCSL1A:
                status = openl1_czcs(l1file);
                break;
            case FMT_HMODISL1B:
                status = openl1_hmodis_hdf(l1file);
                break;
            case FMT_CLASSAVHRR:
                 status = openl1_aci_hdf(l1file);
                 break;
            case FMT_OCML1B:
                 status = openl1_ocm_hdf(l1file);
                 break;
            case FMT_OCM2L1B:
                 status = openl1_ocm2_hdf(l1file);
                 break;
            case FMT_OCML1BDB:
                 status = openl1_ocmdb_hdf(l1file);
                 break;
            case FMT_MERISL1B:
                 status = openl1_meris_N1(l1file);
                 break;
            case FMT_VIIRSL1B:
                if(l1file->geofile[0] == '\0' && strcmp(l1file->input->program_name, "l1info")) {
                    printf("-E- VIIRS processing requires a GEO file.\n");
                    status = 1;
                } else
                    status = openl1_viirs_h5(l1file);
                break;
            case FMT_VIIRSL1BNC:
                if(l1file->geofile[0] == '\0' && strcmp(l1file->input->program_name, "l1info")) {
                    printf("-E- VIIRS processing requires a GEO file.\n");
                    status = 1;
                } else
                    status = openl1b_viirs_nc(l1file);
                break;
            case FMT_VIIRSL1A:
                if(l1file->geofile[0] == '\0' && strcmp(l1file->input->program_name, "l1info")) {
                    printf("-E- VIIRS processing requires a GEO file.\n");
                    status = 1;
                } else
                     status = openl1_viirs_nc(l1file);
                break;
            case FMT_HICOL1B:
                status = openl1_hico_h5(l1file);
                break;
            case FMT_GOCIL1B:
                status = openl1_goci(l1file);
                break;
            case FMT_MERISCC:
                 status = openl1_meris_CC(l1file);
                 break;
            case FMT_OLIL1B:
                 status = openl1_oli(l1file);
                 break;
            case FMT_ORCA:
                 status = openl1_orca(l1file);
                 break;
            case FMT_AVIRIS:
                 status = openl1_aviris(l1file);
                 break;
            case FMT_PRISM:
                 status = openl1_prism(l1file);
                 break;
            case FMT_OLCI:
                status = openl1_olci(l1file);
                break;
            default:
                printf("openl1 - Unknown L1 input file format specifier: %d\n",
                       l1file->format);
                break;
          };

    } else {

        switch (l1file->format) {
            case FMT_L1HDF: 
            case FMT_L1BNCDF: 
                status = openl1_write(l1file, initrec);
                break;
            default:
                printf("Unknown L1 output file format specifier: %d\n",
                       l1file->format);
                break;
          };
    }

    return (status);
}


/* ---------------------------------------------------------------- */
/* Read a specific level 1 record from the file pointed to by the   */
/* input file handle and load the data into the l1 record structure.*/
/* ---------------------------------------------------------------- */
int readl1( filehandle *l1file, int32_t recnum, l1str *l1rec)
{
    int status;
    int32_t ip;

    /* Clear the L1 record */
    init_l1(l1rec, l1file->npix, l1file->nbands);
    l1rec->tilt   = 0.0;
    l1rec->mside  = 0;
    l1rec->detnum = 0;
    l1rec->ndets  = 1;
    /* Altitude of sensor should be set in l1_sensor.c
     * If not set it will be assumed sensor is above atmosphere
     * and no rayleigh correction will be done.
     */
    l1rec->alt = BAD_FLT;

    for (ip=0; ip<l1file->npix; ip++) {
        l1rec->pixnum[ip] = ip;
        l1rec->slot[ip]   = 0;
        l1rec->alpha[ip]  = 0.0;
    }

    l1rec->input          = l1file->input;
    l1rec->landMaskOn     = l1file->input->landmask;
    l1rec->bathMaskOn     = l1file->input->bathmask;
    l1rec->cloudMaskOn    = l1file->input->cloudmask;
    l1rec->glintMaskOn    = l1file->input->glintmask;
    l1rec->stlightMaskOn  = l1file->input->stlightmask;
    l1rec->hiltMaskOn     = l1file->input->hiltmask;
    l1rec->senzMaskOn     = l1file->input->satzenmask;
    l1rec->solzMaskOn     = l1file->input->sunzenmask;

    l1rec->sensorID = l1file->sensorID;
    l1rec->nbands   = l1file->nbands;
    l1rec->nbandsir = l1file->nbandsir;
    l1rec->bindx    = l1file->bindx;
    l1rec->ndets    = l1file->ndets;
    l1rec->nscans   = l1file->nscan;
    l1rec->n_inprods = l1file->n_inprods;
    l1rec->iscan    = recnum;

    l1rec->spix  = MAX(l1file->input->spixl - 1, 0);
    l1rec->epix  = MIN(l1file->input->epixl - 1, l1file->npix-1);
    l1rec->dpix  = MAX(l1file->input->dpixl,1);
    l1rec->sscan = MAX(l1file->input->sline - 1, 0);
    l1rec->escan = MIN(l1file->input->eline - 1, l1file->nscan-1);
    l1rec->dscan = MAX(l1file->input->dline,1);

    switch (l1file->format) {
        case FMT_L1HDF:
            status = readl1_hdf_g(l1file,recnum,l1rec);
            break;
        case FMT_L1BNCDF:
            status = readl1_nc_generic(l1file,recnum,l1rec);
            break;
        case FMT_MOSL1B:
            status = readl1_mos_hdf(l1file,recnum,l1rec);
            break;
        case FMT_SEAWIFSL1A:
            status = readl1a_seawifs(l1file,recnum,l1rec);
            break;
        case FMT_OCTSL1B:
            status = readl1_octs_hdf(l1file,recnum,l1rec);
            break;
        case FMT_OCTSL1A:
            status = readl1_octs_hdf(l1file,recnum,l1rec);
            break;
        case FMT_OSMIL1A:
            status = readl1a_osmi(l1file,recnum,l1rec);
            break;
        /* case FMT_MODISL1B: */
        /*     status = readl1_modis_hdf(l1file,recnum,l1rec); */
        /*     break; */
        case FMT_L1XCAL:
            status = readl1_xcal_hdf(l1file,recnum,l1rec);
            break;
        case FMT_CZCSL1A:
            status = readl1_czcs(l1file,recnum,l1rec);
            break;
        case FMT_HMODISL1B:
            status = readl1_hmodis_hdf(l1file,recnum,l1rec);
            break;
        case FMT_CLASSAVHRR:
             status = readl1_aci_hdf(l1file,recnum,l1rec);
             break;
        case FMT_OCML1B:
             status = readl1_ocm_hdf(l1file,recnum,l1rec);
             break;
        case FMT_OCM2L1B:
             status = readl1_ocm2_hdf(l1file,recnum,l1rec);
             break;
        case FMT_OCML1BDB:
             status = readl1_ocmdb_hdf(l1file,recnum,l1rec);
             break;
        case FMT_MERISL1B:
             status = readl1_meris_N1(l1file,recnum,l1rec);
             break;
        case FMT_MERISCC:
             status = readl1_meris_CC(l1file,recnum,l1rec);
             break;
        case FMT_VIIRSL1B:
             status = readl1_viirs_h5(l1file,recnum,l1rec, 0);
             break;
        case FMT_VIIRSL1BNC:
             status = readl1b_viirs_nc(l1file,recnum,l1rec);
             break;
        case FMT_VIIRSL1A:
             status = readl1_viirs_nc(l1file,recnum,l1rec);
             break;
        case FMT_HICOL1B:
            status = readl1_hico_h5(l1file,recnum,l1rec, 0);
            break;
        case FMT_GOCIL1B:
            status = readl1_goci(l1file,recnum,l1rec,0);
            break;
        case FMT_OLIL1B:
            status = readl1_oli(l1file,recnum,l1rec, 0);
            break;
        case FMT_ORCA:
            status = readl1_orca(l1file,recnum,l1rec);
        //    status = 0;
            break;
        case FMT_AVIRIS:
            status = readl1_aviris(l1file,recnum,l1rec);
        //    status = 0;
            break;
        case FMT_PRISM:
            status = readl1_prism(l1file,recnum,l1rec);
        //    status = 0;
            break;
        case FMT_OLCI:
            status = readl1_olci(l1file,recnum,l1rec);
        //    status = 0;
            break;
        default:
            printf("readl1 - Unknown L1 input file format specifier: %d\n",
                   l1file->format);
            break;
    };


    if (status != 0) {
        fprintf(stderr,
        "-E- %s Line %d: Error reading L1B.\n",
        __FILE__,__LINE__);
        return(1);
    }

    // TODO: we should probably just set l1rec->flags and let  setflagbits set l1rec->navfail
    for (ip=0; ip<l1file->npix; ip++) {
        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 || isnan(l1rec->lon[ip]) ||
            l1rec->lat[ip] <  -91.0 || l1rec->lat[ip] >  91.0 || isnan(l1rec->lat[ip]) ) {
            l1rec->navfail[ip] = 1;
        }
    }

    /* Reduce to sub-sample scan if requested */
    if (l1subpix(l1file,l1rec) != 0)
        return(1);

    setflagbits(0,l1rec,NULL,-1);

    /* update scene meta */
    scene_meta_put(l1rec);

    return (0);
}




/* ---------------------------------------------------------------- */
/* Read lon and lat for a specific level 1 record from the file     */
/* pointed to by the input file handle and load the data into the   */
/* l1 record structure.                                             */
/* ---------------------------------------------------------------- */
int readl1_lonlat( filehandle *l1file, int32_t recnum, l1str *l1rec)
{
    int status;

    l1rec->iscan    = recnum;
    l1rec->sensorID = l1file->sensorID;
    l1rec->nbands   = l1file->nbands;
    l1rec->nbandsir = l1file->nbandsir;

    l1rec->spix  = MAX(l1file->input->spixl - 1, 0);
    l1rec->epix  = MIN(l1file->input->epixl - 1, l1file->npix-1);
    l1rec->dpix  = MAX(l1file->input->dpixl,1);
    l1rec->sscan = MAX(l1file->input->sline - 1, 0);
    l1rec->escan = MIN(l1file->input->eline - 1, l1file->nscan-1);
    l1rec->dscan = MAX(l1file->input->dline,1);

    switch (l1file->format) {
      case FMT_SEAWIFSL1A:
        status = readl1a_lonlat_seawifs(l1file,recnum,l1rec);
        break;
      /* case FMT_MODISL1B: */
      /*   status = readl1_lonlat_modis_hdf(l1file,recnum,l1rec); */
      /*   break; */
      case FMT_HMODISL1B:
        status = readl1_lonlat_hmodis_hdf(l1file,recnum,l1rec);
        break;
      case FMT_MERISL1B:
        status = readl1_lonlat_meris_N1(l1file,recnum,l1rec);
        break;
      case FMT_VIIRSL1B:
           status = readl1_viirs_h5(l1file,recnum,l1rec, 1);
           break;
      case FMT_VIIRSL1BNC:
           status = readl1b_lonlat_viirs_nc(l1file,recnum,l1rec);
           break;
      case FMT_VIIRSL1A:
           status = readl1_lonlat_viirs_nc(l1file,recnum,l1rec);
           break;
      case FMT_HICOL1B:
        status = readl1_hico_h5(l1file,recnum,l1rec, 1);
        break;
      case FMT_GOCIL1B:
        status = readl1_goci(l1file,recnum,l1rec,1);
        break;
      case FMT_OLIL1B:
          status = readl1_oli(l1file,recnum,l1rec, 1);
          break;
      default:
        status = readl1(l1file, recnum, l1rec);
        break;
    };
    
    if (status != 0) {
        fprintf(stderr,
        "-E- %s Line %d: Error reading L1B.\n",
        __FILE__,__LINE__);
        return(1);
    }

    /* Reduce to sub-sample scan if requested */
    if (l1subpix(l1file,l1rec) != 0)
        return(1);

    return (0);
}

