/* =========================================================== */
/* Module read_target.c                                        */
/*                                                             */
/* Functions to open and read a recalibration target file.     */
/*                                                             */
/* Written By:                                                 */
/*                                                             */
/*     B. A. Franz                                             */
/*     SAIC General Sciences Corp.                             */
/*     NASA/SIMBIOS Project                                    */
/*     April 1998                                              */
/*                                                             */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "l12_proto.h"
#include "filehdr_struc.h"
/*
#include "target_struc.h"
#include "filehandle.h"
*/
#include "read_l3bin.h"

static FILE *fp = NULL;
static char current_file[FILENAME_MAX];


/* ----------------------------------------------------------- */
/* close_target() - close recal target file                    */
/* ----------------------------------------------------------- */
void close_target(void)
{
    fclose(fp);
}


/* ----------------------------------------------------------- */
/* read_targethdr() - loads header info                        */
/* ----------------------------------------------------------- */
int read_targethdr(filehandle *file)
{
    filehdr hdr;

    fseek(fp,0,SEEK_SET);

    if ( fread(&hdr,1,sizeof(hdr),fp) != sizeof(hdr) ) {
        printf("-E- %s: File read error.\n",__FILE__);
        return(1);
    }

    if ( hdr.length  < 0 || hdr.length > 1000000 ||
         hdr.npix    < 0 || hdr.npix   > 10000 ) {
        printf("-E- %s: Header values out of range.\n",__FILE__);
        printf("  Record length    = %d\n",hdr.length);
        printf("  Pixels per scan  = %d\n",hdr.npix);
        return(1);
    }        

    file->sensorID = hdr.sensorID;
    file->length   = hdr.length;
    file->npix     = hdr.npix;
    file->format   = hdr.format;
    file->nscan    = hdr.nscan;
    file->mode     = READ;
        
    return(0);
}


/* ----------------------------------------------------------- */
/* open_target() - opens file if not already opened            */
/* ----------------------------------------------------------- */
int open_target(filehandle *file) 
{
    if (fp == NULL || strcmp(file->name,current_file) != 0) {
        if (fp != NULL) close_target();
        if ((fp = fopen(file->name,"r")) == NULL) {
            printf("-E- %s: Error opening %s for reading.\n",
                   __FILE__,file->name);
            return(1);
        }
        strcpy(current_file,file->name);
        if ( read_targethdr(file) != 0 ) {
            printf("-E- %s: Error reading header for %s.\n",
                   __FILE__,file->name);
            return(1);
        }
    }

    return(0);
}


/* ----------------------------------------------------------- */
/* read_target() - reads one recal target record               */
/*                                                             */
/* B. A. Franz, GSC, SIMBIOS Project, March 1998               */
/* ----------------------------------------------------------- */
int read_target( filehandle *file, int32_t recnum, tgstr *tgrec)
{
    /*                                                         */
    /* Open the input file if it is not already open           */
    /*                                                         */
    if ( open_target(file) != 0 ) { 
        printf("-E- %s: File open error.\n",__FILE__);
        return(1);
    } 

    if ( feof(fp) ) {
        printf("-I- %s: End of target file %s reached.",
               __FILE__,file->name);
        return(1);
    }

    if ( fseek( fp,(recnum+1) * file->length,SEEK_SET) != 0 ) {
        printf("-E- %s: Error seeking record %d in %s.",
               __FILE__,recnum,file->name);
        return(1);
    }

    if (fread( tgrec->data,1,file->length,fp ) != file->length) {
        return(1);
    }

    tgrec->sensorID = file->sensorID;
    tgrec->length   = file->length;
    tgrec->npix     = file->npix;

    return(0);
}


/**
 * bin_match - a binary search routine to find the nearest bin in a list of bin numbers.
 * TODO: replace this function with code Don adding to libbin++
 *
 * @param nbins - number of bins
 * @param bins - list of bin numbers
 * @param bin_num - bin number to seach
 * @return
 */
int32_t bin_match(int32_t nbins, int32_t *bins, int32_t bin_num)
{
    int32_t  jl = -1, ju = nbins, jm = 0;
    int32_t  ascnd;
    
    ascnd = (bins[nbins-1] >= bins[0]);
    while ( ju - jl > 1 ) {
   	jm = (ju + jl)/2;
	if (ascnd == (bin_num >= bins[jm]))
	    jl = jm;
	else
	    ju = jm;
    }

    if (bin_num == bins[jl]) return(jl);
    if (jl+1 < nbins && bin_num == bins[jl+1]) return(jl+1);
    if (jl > 0 && bin_num == bins[jl-1]) return(jl-1);
    if (bin_num == bins[0]) return(0);
    if (bin_num == bins[nbins-1]) return(nbins-1);

    return(-1);
}


/**
 * lonlat2bin - returns a L3 bin number given a longitude and latitude
 * TODO: replace this function with code Don adding to libbin++
 * @param l3bin - l3bin structure
 * @param lon - longitude
 * @param lat - latitude
 * @return
 */
int32_t lonlat2bin(l3binstr *l3bin, float lon, float lat) 
{
    int32_t row;
    int32_t col;
    int32_t bin;

    lat = MAX(MIN(lat,90),-90);
    lon = MAX(MIN(lon,180),-180);

    row = MIN(((90 + lat)*l3bin->nrows/180.0),l3bin->nrows - 1);
    col = MIN(((lon + 180.0)*l3bin->numbin[row]/360.0),l3bin->numbin[row]);

    bin = l3bin->basebin[row] + col;
 
    return(bin);
}


/* ----------------------------------------------------------- */
/* read_target_l3() - loads one recal target record            */
/*                                                             */
/* B. A. Franz, GSC, SIMBIOS Project, March 1998               */
/* ----------------------------------------------------------- */

/**
 * read_target_l3 - loads a vicarous calibration target record,
 * applying band shifting if necessary
 *
 * @param file - filehandle structure for L3 bin target file
 * @param l1rec
 * @param nbands
 * @param tgrec
 * @return
 */
int read_target_l3(filehandle *file, l1str *l1rec, int32_t nbands, tgstr *tgrec, initstr *initrec)
{
    static int firstCall = 1;
    static l3binstr l3bin;
    static int32 nwaves;
    static float *l3wvls;
    static float *l3Fonom;
    static int *wvlmatch;
    static int nvisbands;

    int32_t ip, ib, ipb, band;
    int32_t bin, idx, nobs, nscenes;
    float chl, aot;
    static int needBandShift = 0;
    
	float interpband;	// for linear and spline interp
	float *l3valsderiv2; // for spline interp.
	
    if (firstCall) {

        read_l3bin(file->name, &l3bin, l1rec->nbands);
        nwaves = l3bin.nwave;
        if ((l3bin.sensorID != l1rec->sensorID) && (l1rec->input->band_shift_opt != 2)) {
            needBandShift = 1;
        }
        rdsensorinfo(l3bin.sensorID,0,"fwave", (void **) &l3wvls);
        wvlmatch = (int *) malloc(l1rec->nbands * sizeof(int));
        l3Fonom = (float *) malloc(l1rec->nbands * sizeof(int));
        for (ib = 0; ib < l1rec->nbands; ib++) {
            wvlmatch[ib] = windex(l1rec->fwave[ib], l3wvls, nwaves);
            get_f0_thuillier_ext(l3wvls[ib], BANDW, &l3Fonom[ib], initrec->f0rec );
            if (l1rec->iwave[ib] < 700){
                nvisbands++;
            }
        }
        firstCall = 0;
    }
    for (ip = 0; ip < l1rec->npix; ip++) {
        // initialize target rec
        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            tgrec->nLw[ipb] = BAD_FLT;
            tgrec->Lw[ipb] = BAD_FLT;
        }
        tgrec->solz[ip] = BAD_FLT;
        // check L1 masking
        if (l1rec->input->vcal_depth < 0) {
            if (l1rec->elev[ip] > l1rec->input->vcal_depth)
                continue;
        }
        // locate bin 
        bin = lonlat2bin(&l3bin, l1rec->lon[ip], l1rec->lat[ip]);
        idx = bin_match(l3bin.nbins, l3bin.bins, bin);

        // check masking and copy values 
        if (idx >= 0) {
            chl  = l3bin.chl[idx];
            aot  = l3bin.tau[idx];
            nobs = l3bin.nobs[idx];
            nscenes = l3bin.nscenes[idx];
            if (chl <= l1rec->input->chlthreshold && aot <= l1rec->input->aotthreshold
                    && nobs >= l1rec->input->vcal_min_nbin && nscenes >= l1rec->input->vcal_min_nscene) {
                for (ib = 0; ib < nvisbands; ib++) {
                    ipb = ip * nbands + ib;
                    // Only bandshift if L2 wvl is different from L3 wvl even if sensors differ
                    if (needBandShift && (l3wvls[wvlmatch[ib]] != l1rec->fwave[ib]) ) {
                        
                        float l3vals[nwaves];
                        for (band = 0; band < nwaves; band++) {
                            l3vals[band] = l3bin.data[idx][band];
                            if (!l3bin.hasRrs) // Need Rrs for bandshift
                                l3vals[band] /= l3Fonom[ib];
                        }
			            // option 0: linear interpolation. EK: added second condition
			            if (l1rec->input->band_shift_opt==0 || l1rec->fwave[ib] > 700)
                            interpband = linterp(l3wvls,l3vals,nwaves,l1rec->fwave[ib]);
                            
        	            // EK: option 1: bio-optical band shift. added second condition
			            else if (l1rec->input->band_shift_opt==1 && l1rec->fwave[ib] <= 700)
				            interpband = bioBandShift(l3wvls,l3vals,nwaves,l1rec->fwave[ib]);
                           
			            else{
		                    printf("-E- %s:%d: band_shift_opt error.\n",__FILE__,__LINE__);
				            exit(EXIT_FAILURE);
			            }
                       // if (!l3bin.hasRrs)
                         //   interpband *= l3Fonom[ib];
                       // if (l3bin.hasRrs)
                            interpband *= l3Fonom[ib];
						tgrec->nLw[ipb] = interpband;

                        //the nLws in the target files are now in W/m2/nm/sr, so a factor of 10 too big
                        if (l3bin.hasRrs == 0)
                            tgrec->nLw[ipb] /= 10.;
                    } 
					else {
                        if (l3bin.hasRrs == 1)
                            tgrec->nLw[ipb] = l3bin.data[idx][wvlmatch[ib]] * l3Fonom[ib];
                        else
                            tgrec->nLw[ipb] = l3bin.data[idx][wvlmatch[ib]] / 10.;
                    }
                }
                for (ib = nvisbands; ib < nbands; ib++) {
                    ipb = ip * nbands + ib;
                    tgrec->nLw[ipb] = 0.0;
                }
            }
        }
    }

    tgrec->sensorID = l1rec->sensorID;
    tgrec->npix = l1rec->npix;

    return (0);
}
