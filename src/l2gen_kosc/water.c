/* ============================================================================================== */
/* module water.c - functions to read and return absorption and backscatter of pure sea water     */
/*                                                                                                */
/* B. Franz, NASA/GSFC Ocean Color Discipline Processing Group, Sep. 2004                         */
/*                                                                                                */
/* ============================================================================================== */

#include "l12_proto.h"

#define MINAWTAB 200
#define MAXAWTAB 2249
#define INTAWTAB 1
#define NAWTAB  ((MAXAWTAB-MINAWTAB)/INTAWTAB + 1)

static double awtab [NAWTAB];
static double bbwtab[NAWTAB];
static int32_t   ntab   = NAWTAB;
static int    min_wl = MINAWTAB;
static int    del_wl = INTAWTAB;


/* ---------------------------------------------------------------------------------------------- */
/* read_water_spectra() - called once to load look-up table static arrays                         */
/* ---------------------------------------------------------------------------------------------- */
void read_water_spectra(void) 
{
    static int firstCall = 1;

    FILE *fp;
    char *filedir;
    char  filename[FILENAME_MAX];
    char  *line;
    line = (char *) calloc(80,sizeof(char));
    int32_t  i, imin, imax;
    int32_t  j;
    float wl=0, aw=0, bw=0;

    if (!firstCall) return;

    if ((filedir = getenv("OCDATAROOT")) == NULL) {
        printf("-E- %s: OCDATAROOT env variable undefined.\n", __FILE__);
        exit(1);
    }
    strcpy(filename, filedir);
    strcat(filename, "/common/water_spectra.dat");

    if ( (fp = fopen(filename,"r")) == NULL ) {
        fprintf(stderr,
                "-E- %s line %d: unable to open %s for reading\n",__FILE__,__LINE__,filename);
        exit(1);
    }

    i=0;
    while (i < ntab) {
        if ( fgets( line, 80, fp ) == NULL ) {
            fprintf(stderr,"-E- %s line %d: error reading %s at line %d\n",__FILE__,__LINE__,filename,i);
            exit(1);
        }
        if (line[0] != '/' && line[0] !='!') {
            sscanf(line,"%f %f %f",&wl,&aw,&bw);
            awtab [i] = aw;
            bbwtab[i] = bw/2.0;
	    i++;
	}
    }
    free(line);
    firstCall=0;
}


/* ---------------------------------------------------------------------------------------------- */
/* aw_spectra() - returns water absorption at wavelength, wl, averaged over bandwidth, width      */
/* ---------------------------------------------------------------------------------------------- */
float aw_spectra(int32_t wl, int32_t width)
{
    static int firstCall = 1;

    int32_t  itab = (wl - min_wl)/del_wl;
    int32_t  imin = MAX(itab - width/2/del_wl,    0);
    int32_t  imax = MIN(itab + width/2/del_wl, ntab);    
    float aw   = 0;
    int32_t  i;

    if (firstCall) {
        read_water_spectra();
	firstCall = 0;
    }

    if (itab < 0) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside aw table range.\n",__FILE__,__LINE__,wl);
        itab=0;
    } else if (itab > ntab-1) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside aw table range.\n",__FILE__,__LINE__,wl);
        itab=ntab-1;
    }

    for (i=imin; i<=imax; i++)
        aw += (float) awtab[i];

    aw /= (imax-imin+1);

    return(aw);
}
    

/* ---------------------------------------------------------------------------------------------- */
/* bbw_spectra() - returns water backscatter at wavelength, wl, averaged over bandwidth, width    */
/* ---------------------------------------------------------------------------------------------- */
float bbw_spectra(int32_t wl, int32_t width)
{
    static int firstCall = 1;

    int32_t  itab = (wl - min_wl)/del_wl;
    int32_t  imin = MAX(itab - width/2/del_wl,    0);
    int32_t  imax = MIN(itab + width/2/del_wl, ntab);    
    float bbw  = 0;
    int32_t  i;

    if (firstCall) {
        read_water_spectra();
	firstCall = 0;
    }

    if (itab < 0) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside bbw table range.\n",__FILE__,__LINE__,wl);
        itab=0;
    } else if (itab > ntab-1) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside bbw table range.\n",__FILE__,__LINE__,wl);
        itab=ntab-1;
    }

    for (i=imin; i<=imax; i++)
        bbw += (float) bbwtab[i];

    bbw /= (imax-imin+1);

    return(bbw);
}
    

/* ---------------------------------------------------------------------------------------------- */
/* returns aw and bbw at specified "sensor" wavelengths, appropriate to the derived nLw           */
/* ---------------------------------------------------------------------------------------------- */
void get_aw_bbw(l2str *l2rec,float wave[],int nwave,float *aw,float *bbw)
{
    int ib, iw;

    if (l2rec->input->outband_opt >= 2) {
        for (ib=0; ib<nwave; ib++) {
            aw [ib] = aw_spectra (wave[ib],BANDW);
            bbw[ib] = bbw_spectra(wave[ib],BANDW);
        }
    } else {
        float *senaw;
        float *senbbw;
        rdsensorinfo(l2rec->sensorID,l2rec->input->evalmask,"aw", (void **) &senaw );
        rdsensorinfo(l2rec->sensorID,l2rec->input->evalmask,"bbw",(void **) &senbbw);
        for (ib=0; ib<nwave; ib++) {
	    iw = bindex_get(wave[ib]);
            aw [ib] = senaw [iw];
            bbw[ib] = senbbw[iw];
	}
    }
}
