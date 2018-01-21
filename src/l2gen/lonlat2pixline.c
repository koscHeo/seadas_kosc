#include "lonlat2pixline.h"

#include <readL2scan.h>

#include <stdio.h>
#include <string.h>
#include <math.h>

/* this function is used when comparing lons to see if the lon is
 * inside the requested lon range.  The lons are rotated such that
 * the first lon in the range is zero then normalized from -180 to 180 */

static float normalizeLon(float lon) {
    while (lon < -180.0)
        lon += 360.0;
    while (lon >= 180.0)
        lon -= 360.0;
    return lon;
}

static float rotateLon(float lon, float ref) {
    return normalizeLon(lon - ref);
}

static void fixModisResolution(lonlat2pixline_t *params) {
    int factor = 1;
    if(params->resolution == 500) {
        factor = 2;
        params->spixl = (params->spixl-1) * factor + 1;
        params->epixl = (params->epixl-1) * factor + 1;
        params->sline = (params->sline-1) * factor + 1;
        params->eline = (params->eline-1) * factor + 1;
    } else if(params->resolution == 250) {
        factor = 4;
        params->spixl = (params->spixl-1) * factor + 1;
        params->epixl = (params->epixl-1) * factor + 1;
        params->sline = (params->sline-1) * factor + 1;
        params->eline = (params->eline-1) * factor + 1;
    }
}

/*
 * Note that if a MODIS GEO file is given with a resolution of 500 or
 * 250 the pixl and line could be off by 1 or 3 respectively.
 */
int lonlat2pixline(lonlat2pixline_t *params) {

    int32_t iscan = 0; /* input scan number                  */
    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t escan = -1; /* end scan for subscene process      */

    int32_t npix = 0; /* Number of output pixels per scan   */
    int32_t nscan = 0; /* Number of output scans             */

    float cornerlon[2]; /* longitude corners of box */
    float cornerlat[2]; /* latitude corners of box */

    float filelon[2] = {1000, -1000}; /* longitude corners of file */
    float filelat[2] = {1000, -1000}; /* latitude corners of file */

    int32_t status = 0;

    l1str *l1rec;               /* generic level-1b scan structure */
    filehandle *l1file;         /* input file handle */
    instr *input;

    int32_t ipix;
    float lon;
    float lonMiddle;
    float lat;
    double uv[3];
    double uvpix[3];
    double maxcos = -1;
    double dot;
    float pixlon;
    float pixlat;
    float normBoxLon[2];
    float refBoxLon = 0.0; // set to the middle of the requested box

    int16 lonTest;
    int16 latTest;

    int32 minScan = 2147483647;
    int32 maxScan = -1;
    int32 minPix = 2147483647;
    int32 maxPix = -1;
    int32 pixsn = -1;
    int32 pixpx = -1;

    int32 sd_id;
    int32 attrindx;
    int32 extract_pix_off = 0;

    int32 percentDone;
    int32 percentFloor = 0;

    int32 start[3] = {0, 0, 0};
    int32 edges[3];

    static l2_prod l2_str;

    char buffer[1024];

    int32_t type;

    int32 sds_id_ll[2];
    int32 dims[2];
    int32 rank, dtype, nattrs;

    float latlon[2][2000];
    float *lonArray;
    float *latArray;
    
    // allocate structures
    l1rec = (l1str*) malloc(sizeof(l1str));
    l1file = (filehandle*) malloc(sizeof(filehandle));
    input = (instr*) malloc(sizeof(instr));
    
    
    filehandle_init(l1file);
    msl12_input_init(input);
    input->resolution = params->resolution;
    if (input->resolution == -1000) {
        input->resolution = 1000;
        input->modis_subset_compat = 1;
        l1file->modis_subset_compat = 1;
    }

    if (params->pix_srch) {

        /* Generate corner lon/lat if single pixel search */
        /* ---------------------------------------------- */

        cornerlon[0] = params->SWlon - 0.025;
        cornerlon[1] = params->SWlon + 0.025;

        cornerlat[0] = params->SWlat - 0.025;
        cornerlat[1] = params->SWlat + 0.025;
        if (cornerlat[0] < -90) cornerlat[0] = -89.99;
        if (cornerlat[1] > +90) cornerlat[1] = +89.99;

        uvpix[0] = cos(params->SWlat / RADEG) * cos(params->SWlon / RADEG);
        uvpix[1] = cos(params->SWlat / RADEG) * sin(params->SWlon / RADEG);
        uvpix[2] = sin(params->SWlat / RADEG);
    } else {
        cornerlon[0] = params->SWlon;
        cornerlon[1] = params->NElon;
        cornerlat[0] = params->SWlat;
        cornerlat[1] = params->NElat;

        if (cornerlat[0] > cornerlat[1]) {
            printf("-E- lonlat2pixline: SW lat must be south of NE lat\n");
            exit(1);
        }
    }

    /* normalize to +-180 */
    cornerlon[0] = normalizeLon(cornerlon[0]);
    cornerlon[1] = normalizeLon(cornerlon[1]);
    if(cornerlon[1] < cornerlon[0])
        cornerlon[1] += 360;

    // setup initial lon for comparing lons
    refBoxLon = normalizeLon(cornerlon[0] + (cornerlon[1] - cornerlon[0]) / 2);
    normBoxLon[0] = rotateLon(cornerlon[0], refBoxLon);
    normBoxLon[1] = rotateLon(cornerlon[1], refBoxLon);

    strcpy(l1file->name, params->input_filename);
    l1file->geofile = params->geo_filename;
    type = getFormat(l1file);

    if (type == FMT_MODISL1B && (params->resolution == 500 || params->resolution == 250)) {
        printf("-E- lonlat2pixline: Can not have resolution of %d with an Ocean Subsetted MODIS file.\n",
                params->resolution);
        exit(1);
    }

    if (type == FMT_MODISGEO) {

        sd_id = SDstart(params->input_filename, DFACC_RDONLY);

        if (sd_id == -1) {
            printf("-E- lonlat2pixline: Error opening %s for reading.\n",
                    params->input_filename);
            exit(1);
        }

        sds_id_ll[0] = SDselect(sd_id, SDnametoindex(sd_id, "Latitude"));
        if (sds_id_ll[0] == -1) {
            printf("-E- lonlat2pixline: Error opening Latitude field.\n");
            exit(1);
        }
        status = SDgetinfo(sds_id_ll[0], buffer, &rank, dims, &dtype, &nattrs);

        nscan = dims[0];
        npix = dims[1];

        status = SDreadattr(sd_id,
                SDfindattr(sd_id, "Max Earth Frames"),
                &epix);

        epix--;
        escan = nscan - 1;

        sds_id_ll[1] = SDselect(sd_id, SDnametoindex(sd_id, "Longitude"));
        if (sds_id_ll[1] == -1) {
            printf("-E- lonlat2pixline: Error opening Longitude field.\n");
            exit(1);
        }

        /* Get extract pixel offset if present */
        attrindx = SDfindattr(sd_id, "Extract Pixel Offset");
        if (attrindx != -1)
            status = SDreadattr(sd_id, attrindx, &extract_pix_off);


    } else if (type == FMT_L2HDF || type == FMT_L2NCDF || type == FMT_L1BNCDF) {

        status = openL2(params->input_filename, 0x0, &l2_str);
        if (status) {
            printf("-E- lonlat2pixline: Error opening %s for reading.\n",
                    params->input_filename);
            exit(EXIT_FAILURE);
        }

        meta_l2Type meta_l2;
        status = readL2meta(&meta_l2, l2_str.fileindex);
        if (status) {
            printf("-E- lonlat2pixline: Error opening %s for reading.\n",
                    params->input_filename);
            exit(EXIT_FAILURE);
        }

        escan = l2_str.nrec;
        epix = l2_str.nsamp;

        filelon[0] = normalizeLon(meta_l2.westlon);
        filelon[1] = normalizeLon(meta_l2.eastlon);
        if(filelon[1] < filelon[0])
            filelon[1] += 360;
        filelat[0] = meta_l2.southlat;
        filelat[1] = meta_l2.northlat;

        escan--;
        epix--;

        npix = epix - spix + 1;
        nscan = escan - sscan + 1;

        // check if the box totally contains the file
        if (normBoxLon[0] <= rotateLon(filelon[0], refBoxLon) &&
                normBoxLon[1] >= rotateLon(filelon[1], refBoxLon) &&
                cornerlat[0] <= filelat[0] &&
                cornerlat[1] >= filelat[1]) {
            params->sline = 1;
            params->eline = escan + 1;
            params->spixl = 1;
            params->epixl = epix + 1;

            // close resources
            closeL2(&l2_str, 0);
            freeL2(&l2_str);
            freeL2(NULL);

            status = 120;
            goto bail;
        }

    } else if (type > 0) {

        /*								*/
        /* Set input defaults (e.g., calfile)				*/
        /*								*/
        msl12_input_defaults(l1file, input);
        if (openl1(l1file) != 0) {
            printf("-E- lonlat2pixline: Error opening %s for reading.\n",
                    params->input_filename);
            exit(1);
        }

        /*								*/
        /* Allocate memory for L1 scan data			  	*/
        /*								*/
        if (alloc_l1(l1file->npix, l1file->nbands, NBANDSIR, l1file->n_inprods, l1rec) == 0) {
            printf("-E- lonlat2pixline: Unable to allocate L1 record.\n");
            exit(FATAL_ERROR);
        }

        epix = l1file->npix - 1;
        escan = l1file->nscan - 1;

        l1file->spix = spix;
        l1file->epix = epix;
        l1file->dpix = 1;

        npix = epix - spix + 1;
        nscan = escan - sscan + 1;

    } else {
        printf("-E- lonlat2pixline: File type not recognized for %s.\n",
                params->input_filename);
        exit(-1);
    }

    /*					 			*/
    /* Read file scan by scan                                   */
    /*								*/

    // start out assuming whole box is in the file
    params->box_failed = 0;

    /* The stderr output is intended for seadas so users can see     */
    /* some progress if they're processing a large MERIS file. (MAR) */
    if (want_verbose)
        fprintf(stderr, "\n");
    for (iscan = sscan; iscan <= escan; iscan++) {

        if (type == FMT_MODISGEO) {
            start[0] = iscan;
            start[1] = 0;
            edges[0] = 1;
            edges[1] = dims[1];
            status = SDreaddata(sds_id_ll[0], start, NULL, edges,
                    (VOIDP) latlon[0]);
            status = SDreaddata(sds_id_ll[1], start, NULL, edges,
                    (VOIDP) latlon[1]);
            lonArray = latlon[1];
            latArray = latlon[0];

        } else if (type == FMT_L2HDF || type == FMT_L2NCDF || type == FMT_L1BNCDF) {
            start[0] = iscan;
            start[1] = 0;
            edges[0] = 1;
            edges[1] = npix;
            readlonlat(&l2_str, 0, start, edges, NULL);
            lonArray = l2_str.longitude;
            latArray = l2_str.latitude;
            
        } else {
            readl1_lonlat(l1file, iscan, l1rec);
            lonArray = l1rec->lon;
            latArray = l1rec->lat;
        }

        //read midpoint of line
        lonMiddle = normalizeLon(lonArray[npix / 2]);
        
        // put min max lon into reasonable range if not set by metadata
        if(filelon[0] == 1000)
            filelon[0] = lonMiddle;
        if(filelon[1] == -1000)
            filelon[1] = lonMiddle;

        // see if endpoints of scan line are inside the box
        // this would cause less than a full box to be extracted
        if (!params->box_failed) {
            // check whole first and last line
            if (iscan == sscan || iscan == escan - 1) {
                for (ipix = 0; ipix < npix; ipix++) {
                    if (normBoxLon[0] < rotateLon(lonArray[ipix], refBoxLon) &&
                            normBoxLon[1] > rotateLon(lonArray[ipix], refBoxLon) &&
                            cornerlat[0] < latArray[ipix] &&
                            cornerlat[1] > latArray[ipix]) {
                        params->box_failed = 1;
                    }
                }
            } else {
                // check first pixel
                if (normBoxLon[0] < rotateLon(lonArray[0], refBoxLon) &&
                        normBoxLon[1] > rotateLon(lonArray[0], refBoxLon) &&
                        cornerlat[0] < latArray[0] &&
                        cornerlat[1] > latArray[0]) {
                    params->box_failed = 1;
                }
                // check last pixel
                if (normBoxLon[0] < rotateLon(lonArray[npix-1], refBoxLon) &&
                        normBoxLon[1] > rotateLon(lonArray[npix-1], refBoxLon) &&
                        cornerlat[0] < latArray[npix-1] &&
                        cornerlat[1] > latArray[npix-1]) {
                    params->box_failed = 1;
                }
            }
        }

        
        for (ipix = 0; ipix < npix; ipix++) {
            lat = latArray[ipix];
            lon = normalizeLon(lonArray[ipix]);

            // find file min,max lon and lat
            if (rotateLon(lon, lonMiddle) > 0) {
                // check against max
                if (rotateLon(lon, lonMiddle) > rotateLon(filelon[1], lonMiddle)) {
                    filelon[1] = lon;
                }
            } else {
                // check against min
                if (rotateLon(lon, lonMiddle) < rotateLon(filelon[0], lonMiddle)) {
                    filelon[0] = lon;
                }
            }

            if (lat < filelat[0]) filelat[0] = lat;
            if (lat > filelat[1]) filelat[1] = lat;

            latTest = (lat >= cornerlat[0] && lat <= cornerlat[1]);
            lonTest = (rotateLon(lon, refBoxLon) >= normBoxLon[0] && rotateLon(lon, refBoxLon) <= normBoxLon[1]);

            if (lonTest && latTest) {
                if (ipix < minPix) minPix = ipix;
                if (iscan < minScan) minScan = iscan;
                if (ipix > maxPix) maxPix = ipix;
                if (iscan > maxScan) maxScan = iscan;
            }

        } // for ipix

        if (params->pix_srch && minPix != 2147483647 && maxPix != -1) {

            for (ipix = minPix; ipix <= maxPix; ipix++) {
                lat = latArray[ipix] / RADEG;
                lon = lonArray[ipix] / RADEG;

                uv[0] = cos(lat) * cos(lon);
                uv[1] = cos(lat) * sin(lon);
                uv[2] = sin(lat);

                dot = uv[0] * uvpix[0] + uv[1] * uvpix[1] + uv[2] * uvpix[2];
                if (dot > maxcos) {
                    maxcos = dot;
                    pixsn = iscan;
                    pixpx = ipix;
                    pixlon = lon * RADEG;
                    pixlat = lat * RADEG;
                }
            }
        } // if pix search

        /* Moved this down here so that the output from readl1() wouldn't */
        /* screw up the stderr output.  Looks okay now when run both from */
        /* the unix command-line and seadas.                              */
        if (want_verbose) {
            percentDone = ((float) iscan / (float) escan)*100.0;
            if (percentDone >= percentFloor) {
                percentFloor += 5;
                fprintf(stderr, "\rlonlat2pixline %2.0lf%% done.",
                        ((float) iscan / (float) escan)*100.0);
                fflush(stderr);
            }
            if (iscan % 500 == 0) {
                fprintf(stderr, "  Reading scan: %d out of %d", iscan, escan);
                fflush(stderr);
            }
        }
    } // for iscan

    if (want_verbose) {
        fprintf(stderr, "\n");
        fflush(stderr);
    }

    if (minScan == 2147483647) minScan = -1;
    if (minPix == 2147483647) minPix = -1;

    /* If extract_pixel_offset > 0 then subtract from pixel limits */
    if (extract_pix_off > 0) {
        minPix -= extract_pix_off;
        maxPix -= extract_pix_off;
    }

    if (params->pix_srch) {
        if (pixsn != -1) {
            if (params->want_pixbox == 1) {
                params->pixLon = pixlon;
                params->pixLat = pixlat;

                ((pixsn + 1 - params->ybox) >= 1) ?
                        (params->sline = pixsn + 1 - params->ybox) :
                        (params->sline = 1);

                ((pixsn + 1 + params->ybox) <= escan) ?
                        (params->eline = pixsn + 1 + params->ybox) :
                        (params->eline = escan);

                ((pixpx + 1 - extract_pix_off - params->xbox) >= 1) ?
                        (params->spixl = pixpx + 1 - extract_pix_off - params->xbox) :
                        (params->spixl = 1);

                ((pixpx + 1 - extract_pix_off + params->xbox) <= epix) ?
                        (params->epixl = pixpx + 1 - extract_pix_off + params->xbox) :
                        (params->epixl = epix);
            } else {
                params->pixLon = pixlon;
                params->pixLat = pixlat;
                params->sline = params->eline = pixsn + 1;
                params->spixl = params->epixl = pixpx + 1 - extract_pix_off;
            }
        } else {
            params->pixLon = 0;
            params->pixLat = 0;
            params->sline = 0;
            params->eline = 0;
            params->spixl = 0;
            params->epixl = 0;
        }
    } else {
        params->sline = minScan + 1;
        params->eline = maxScan + 1;
        params->spixl = minPix + 1;
        params->epixl = maxPix + 1;
    }

    if (type == FMT_MODISGEO) {
        SDendaccess(sds_id_ll[0]);
        SDendaccess(sds_id_ll[1]);
        SDend(sd_id);
    } else if (type == FMT_L2HDF || type == FMT_L2NCDF || type == FMT_L1BNCDF) {
        closeL2(&l2_str, 0);
        freeL2(&l2_str);
        freeL2(NULL);
    } else {
        closel1(l1file);
    }

    if (minScan == -1 || maxScan == -1 || minPix == -1 || maxPix == -1) {
        // if the box is too small to have any points in it try a
        // pixel search
        if (!params->pix_srch &&
                normBoxLon[0] >= rotateLon(filelon[0], refBoxLon) &&
                normBoxLon[1] <= rotateLon(filelon[1], refBoxLon) &&
                cornerlat[0] >= filelat[0] &&
                cornerlat[1] <= filelat[1]) {
            params->pix_srch = 1;
            status = lonlat2pixline(params);
            params->pix_srch = 0;
            goto bail;
        } else {
            status = 100;
            goto bail;
        }
    }
   
    if (params->want_pixbox) {
        if (params->pix_srch) {
            if ((pixsn + 1 - params->ybox) < 1)
                params->box_failed = 1;
            if ((pixpx + 1 - extract_pix_off - params->xbox) < 1)
                params->box_failed = 1;
            if ((pixsn + 1 + params->ybox) > escan)
                params->box_failed = 1;
            if ((pixpx + 1 - extract_pix_off + params->xbox) > epix)
                params->box_failed = 1;
        } else {
            int x, y;
            int dx = (params->epixl - params->spixl) / 2;
            int dy = (params->eline - params->sline) / 2;
            if (dx < params->xbox) {
                x = params->spixl + dx;
                params->spixl = x - params->xbox;
                params->epixl = x + params->xbox;
                if (params->spixl < 1) {
                    params->spixl = 1;
                    params->box_failed = 1;
                }
                if (params->epixl > npix) {
                    params->epixl = npix;
                    params->box_failed = 1;
                }
            } // dx to small
            if (dy < params->ybox) {
                y = params->sline + dy;
                params->sline = y - params->ybox;
                params->eline = y + params->ybox;
                if (params->sline < 1) {
                    params->sline = 1;
                    params->box_failed = 1;
                }
                if (params->eline > nscan) {
                    params->eline = nscan;
                    params->box_failed = 1;
                }
            } // dy too snall
        } // not a pixel search
    } // want a pixbox

    // check if the file is totally contained in the box
    if(params->sline == 1 && params->eline == escan + 1 &&
        params->spixl == 1 && params->epixl == epix + 1) {
        status = 120;
        goto bail;
    }

    // check if the box limits goes outside of the file limits
    if (normBoxLon[0] < rotateLon(filelon[0], refBoxLon) ||
            normBoxLon[1] > rotateLon(filelon[1], refBoxLon) ||
            cornerlat[0] < filelat[0] ||
            cornerlat[1] > filelat[1]) {
        params->box_failed = 1;
    }

    // check if the box goes outside of the file
    if (params->box_failed)
        status = 110;

    // yea, goto is ugly, but it does the trick!!!
bail:

    if (type == FMT_MODISGEO)
       fixModisResolution(params);

// deallocate structures
    free(l1rec);
    free(l1file);
    free(input);

    return (status);
}

int lonlat2pixline1(char *input_filename, char *geo_filename,
        int32_t resolution, float SWlon, float SWlat, float NElon, float NElat,
        int32_t *spixl, int32_t *epixl, int32_t *sline, int32_t * eline) {
    int result;
    lonlat2pixline_t params;

    strcpy(params.input_filename, input_filename);
    strcpy(params.geo_filename, geo_filename);
    params.resolution = resolution;
    params.SWlon = SWlon;
    params.NElon = NElon;
    params.SWlat = SWlat;
    params.NElat = NElat;
    params.pix_srch = 0;
    params.want_pixbox = 0;

    result = lonlat2pixline(&params);
    if (result == 0 || result == 110 || result == 120) {
        *spixl = params.spixl;
        *epixl = params.epixl;
        *sline = params.sline;
        *eline = params.eline;
    }
    return result;
}

int lonlat2pixline2(char *input_filename, char *geo_filename,
        int32_t resolution, float lon, float lat, int32_t dx, int32_t dy,
        int32_t *spixl, int32_t *epixl, int32_t *sline, int32_t * eline) {
    int result;
    lonlat2pixline_t params;

    strcpy(params.input_filename, input_filename);
    strcpy(params.geo_filename, geo_filename);
    params.resolution = resolution;
    params.SWlon = lon;
    params.SWlat = lat;
    params.xbox = dx;
    params.ybox = dy;
    params.pix_srch = 1;
    params.want_pixbox = 1;

    result = lonlat2pixline(&params);
    if (result == 0 || result == 110 || result == 120) {
        *spixl = params.spixl;
        *epixl = params.epixl;
        *sline = params.sline;
        *eline = params.eline;
    }
    return result;
}

int lonlat2pixline3(char *input_filename, char *geo_filename,
            int32_t resolution, float lon, float lat, 
            int32_t *pixl, int32_t *line) {
    int result;
    lonlat2pixline_t params;

    strcpy(params.input_filename, input_filename);
    strcpy(params.geo_filename, geo_filename);
    params.resolution = resolution;
    params.SWlon = lon;
    params.SWlat = lat;
    params.pix_srch = 1;
    params.want_pixbox = 0;

    result = lonlat2pixline(&params);
    if (result == 0 || result == 110 || result == 120) {
        *pixl = params.spixl;
        *line = params.sline;
    }
    return result;
}

