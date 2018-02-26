#include "hdf.h"
#include "mfhdf.h"
#include "l1a.h"
#include "navigation.h"
#include "l1a_proto.h"
#include "eng_qual.h"
#include "l12_proto.h"
#include "l1a_seawifs.h"
#include "call1a_proto.h"
#include "getcal_proto.h"
#include "st_proto.h"

#define   LAC_PIXEL_NUM           1285
#define   GAC_PIXEL_NUM           248
#define   NREC_IN_BUF             10
#define   STBUFSIZ                5
#define   NOTDONE                 0
#define   FIRST_KNEE	          1
#define   MASK_HIGHLT1            16
#define   GENBUFSIZ               NREC_IN_BUF*sizeof(float)*40 /* size of inst_ana */

static int evalmask;

int16   syear, sday;       /* data start date                  */
int32   smsec;             /* data start time                  */
int16   eyear, eday;       /* data end date                    */
int32   emsec;             /* data end time                    */
int32   nscan;             /* number of scans                  */
int32   npix;              /* number pixels per scan           */
int32   spix;              /* start pixel number (from 0)      */
int32   dpix;              /* LAC pixel increment              */
float   dark_mean[8];	   /* Mean of dark restore counts      */
float   dark_std[8];	   /* Std dev  of dark restore counts  */

int16   tdi[BANDS_DIMS_1A] = {0,0,0,0,0,0,0,0};

int16           entry_year;
int16           entry_day;
int16           ref_year;
int16           ref_day;
int16           ref_minute;
float32         fp_temps[256][BANDS_DIMS_1A];
float32         scan_mod[2][1285];
float32         counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
float32         rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];

float64         t_const[BANDS_DIMS_1A];
float64         t_linear_1[BANDS_DIMS_1A];
float64         t_exponential_1[BANDS_DIMS_1A];
float64         t_linear_2[BANDS_DIMS_1A];
float64         t_exponential_2[BANDS_DIMS_1A];
float64         cal_offs[BANDS_DIMS_1A];
float64         inst_tcorr[BANDS_DIMS_1A];
float64         inst_tref[BANDS_DIMS_1A];
float64         fp_tcorr[BANDS_DIMS_1A];
float64         fp_tref[BANDS_DIMS_1A];
float64         ms1_const[BANDS_DIMS_1A];
float64         ms1_linear_1[BANDS_DIMS_1A];
float64         ms1_exponential_1[BANDS_DIMS_1A];
float64         ms1_linear_2[BANDS_DIMS_1A];
float64         ms1_exponential_2[BANDS_DIMS_1A];
float64         ms2_const[BANDS_DIMS_1A];
float64         ms2_linear_1[BANDS_DIMS_1A];
float64         ms2_exponential_1[BANDS_DIMS_1A];
float64         ms2_linear_2[BANDS_DIMS_1A];
float64         ms2_exponential_2[BANDS_DIMS_1A];

char            cal_path_tab[128];
char            dtype[8];

extern float32 cal_counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
extern float32 cal_rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];


float32 pcal_counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
float32 pcal_rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
short l2_flags_buffer[LAC_PIXEL_NUM];

int32           nsta;
int32           ninc;

int16   *l1a_data     = NULL;  /* raw radiances band-int-by-pix    */
int16   *l1a_back     = NULL;  /* raw radiances band-int-by-pix    */
float   *l1b_buffer   = NULL;  /* l1b radiances band-int-by-pix    */
int32   *msec;
int16   *side;
int16   *dark_rest;
float32 *tilt;

float32 ylat[LAC_PIXEL_NUM];
float32 xlon[LAC_PIXEL_NUM];
float32 solz[LAC_PIXEL_NUM];
float32 sola[LAC_PIXEL_NUM];
float32 senz[LAC_PIXEL_NUM];
float32 sena[LAC_PIXEL_NUM];


int16   stray_light = -1;     /* num pixels to apply sl correct.  */
float   Ltyp_frac   = 0.25;   /* fraction of B8 for sl threshold  */
int16   out_band    = 0;      /* if <= 0, no out of band correct. */


typedef struct {
  int16   gain[8];
  int16   tdi[8];
  int16   scan_temp[8];
  float32 inst_temp[40];
  float32 orb_vec[3];
  float32 sun_ref[3];
  float32 sen_mat[3*3];
  float32 scan_ell[6];
} inputBuffer;


int32 do_st = 1;


int32 get_l1a_rec(int32 sd_id, int32 recno, cal_mod_struc *cal_mod,
		  int16 *l1a_dum, float32 **l1b_data, int16 **l2_flags)
{
    /*                                                                 */
    /* local vars                                                      */
    /*                                                                 */
    int32   ipix;                  /* pixel number                     */
    int32   idet;                  /* detector number                  */
    int32   irec;
    int32   start[3]={0,0,0};
    int32   edges[3];

    int16   band_buf[LAC_PIXEL_NUM * 8];
    int32   scan_no;
    int16   gain8;
    int32   AS_pixels;
    float   Styp_frac = 0.9;
    int32   sl_scan;
    int16   *l1a_ptr;
    int16   *gain_ptr;

    static inputBuffer rdBuf[NREC_IN_BUF];
    static inputBuffer bkBuf[STBUFSIZ-2];
    byte   genBuf[GENBUFSIZ];

    static int32 crec = -1;
    static int16 recursive_flag = 0;
    static int16 n_read = 0;
    static int32 max_rec_in_rdbuf;
    static int32 offset;
    static int32 i,j;
    static int32 stray_light_scan_no;

    static float *st_l1b_data;
    static int16 *st_l2_flags;

    static int32 initial = 1;
    static byte  first = 1;


    static int16 hi_Lt[LAC_PIXEL_NUM];  /* whether radiance > knee */
    int	  pixvalue;/* l1a pixel value			*/
    int	  pixdark;/* dark_restore value for each pixel	*/
    float kneevalue; /* radiance at first knee		*/

    /*                                                                 */
    /* Read L1A data scan by scan, build L1B record, and write.        */

    /* If reading out-of-sequence, force restart from beginning to get stlight */
    /* baf, 05-oct-2001.                                                       */
    if (first || recno < crec) {
      max_rec_in_rdbuf = recno - 1;
      stray_light_scan_no = recno;
    }

    crec = recno;

    if (crec > nscan) return 0;


    if (crec > max_rec_in_rdbuf) {

      memcpy(&bkBuf, &rdBuf[NREC_IN_BUF - STBUFSIZ + 2], 
	     sizeof(inputBuffer) * (STBUFSIZ - 2));
      memcpy(l1a_back, &l1a_data[(NREC_IN_BUF - STBUFSIZ + 2)*npix*8], 
	     npix * 8 * (STBUFSIZ - 2) * sizeof(int16));

      n_read = (nscan - max_rec_in_rdbuf >= NREC_IN_BUF) ? 
	NREC_IN_BUF : nscan - max_rec_in_rdbuf;


      edges[0] = 1;
      start[1] = 0;
      edges[1] = npix;
      edges[2] = 8;
      for (irec=0;irec<n_read;irec++) {
	start[0] = recno-1 + irec;
	SDreaddata(SDselect(sd_id, SDnametoindex(sd_id,"l1a_data")), 
			    start, NULL, edges, (VOIDP) band_buf);
	for (ipix=0;ipix<npix;ipix++)
	  for (idet=0;idet<8;idet++) 
	    l1a_data[irec*(8*npix)+idet*npix+ipix] = 
	      band_buf[ipix*8+idet];
      }


      rdSDS(sd_id,"gain",recno-1,0,n_read,8,(VOIDP) genBuf);
      for (i=0;i<n_read;i++) 
	memcpy(&rdBuf[i].gain,&genBuf[i*8*2],8*2);
 
      rdSDS(sd_id,"tdi",recno-1,0,n_read,8,(VOIDP) genBuf);
      for (i=0;i<n_read;i++) 
	memcpy(&rdBuf[i].tdi,&genBuf[i*8*2],8*2);

      rdSDS(sd_id,"scan_temp",recno-1,0,n_read,8,(VOIDP) genBuf); 
      for (i=0;i<n_read;i++) 
	memcpy(&rdBuf[i].scan_temp,&genBuf[i*8*2],8*2);

      rdSDS(sd_id,"inst_ana",recno-1,0,n_read,40,(VOIDP) genBuf); 
      for (i=0;i<n_read;i++)
        memcpy(&rdBuf[i].inst_temp,&genBuf[i*40*4],40*4);

      rdSDS(sd_id,"orb_vec",recno-1,0,n_read,3,(VOIDP) genBuf); 
      for (i=0;i<n_read;i++) 
	memcpy(&rdBuf[i].orb_vec,&genBuf[i*3*4],3*4);

      rdSDS(sd_id,"sun_ref",recno-1,0,n_read,3,(VOIDP) genBuf); 
      for (i=0;i<n_read;i++) 
	memcpy(&rdBuf[i].sun_ref,&genBuf[i*3*4],3*4);

      rdSDS(sd_id,"scan_ell",recno-1,0,n_read,6,(VOIDP) genBuf); 
      for (i=0;i<n_read;i++) 
	memcpy(&rdBuf[i].scan_ell,&genBuf[i*6*4],6*4);


      start[0] = recno-1;
      edges[0] = n_read;
      start[1] = 0;
      edges[1] = 3;
      edges[2] = 3;
      SDreaddata(SDselect(sd_id, SDnametoindex(sd_id,"sen_mat")), 
			  start, NULL, edges, (VOIDP) genBuf);
      for (i=0;i<n_read;i++) 
	memcpy(&rdBuf[i].sen_mat,&genBuf[i*9*4],9*4);


      offset = max_rec_in_rdbuf + 1;
      max_rec_in_rdbuf += n_read;
    }

    if ((crec - offset) >= 0) {
      j = crec - offset;
      gain8 = rdBuf[j].gain[7];
      gain_ptr = rdBuf[j].gain;
      l1a_ptr = &l1a_data[j*8*npix];

      l1b_rad(syear,sday,smsec,msec[recno-1],
	   	    dtype, nsta, ninc, npix, 
		    dark_mean, rdBuf[j].gain, rdBuf[j].tdi, 
		    rdBuf[j].scan_temp, rdBuf[j].inst_temp, side[recno-1],
		    &l1a_data[j*8*npix], l1b_buffer, cal_mod);
/*
          calibrate_l1a(cal_path_tab,syear,sday,smsec,eday,msec[recno-1],
	   	    dtype, nsta, ninc, npix, 
		    dark_mean, rdBuf[j].gain, rdBuf[j].tdi, 
		    rdBuf[j].scan_temp, rdBuf[j].inst_temp, side[recno-1],
		    &l1a_data[j*8*npix], l1b_buffer, cal_mod);
*/

      if (!recursive_flag) {
	geonav_(rdBuf[j].orb_vec,rdBuf[j].sen_mat,rdBuf[j].scan_ell,
		rdBuf[j].sun_ref,&nsta,&ninc,&npix,
		ylat,xlon,solz,sola,senz,sena);

      }
    }	else {
      j = (STBUFSIZ - 2) + (crec - offset);
      gain8 = bkBuf[j].gain[7];
      gain_ptr = bkBuf[j].gain;
      l1a_ptr = &l1a_back[j*8*npix];

      if (!recursive_flag) { 
	geonav_(bkBuf[j].orb_vec,bkBuf[j].sen_mat,bkBuf[j].scan_ell,
		bkBuf[j].sun_ref,&nsta,&ninc,&npix,
		ylat,xlon,solz,sola,senz,sena);

      }
    }

    if (first) {
      memcpy(pcal_counts,cal_counts,sizeof(cal_counts));
      memcpy(pcal_rads,cal_rads,sizeof(cal_rads));
      first = 0;
    }

    if (!recursive_flag) {
      for (ipix=0; ipix < npix ; ipix++) {
	hi_Lt[ipix] = FALSE;
        /* HILT check limited to bands 7 & 8, BAF, 23 Jan 2002 */
        /* Use dark_mean rather than dark_rest, BAF, 7 July 2003 */
	for(idet=6;(idet<8)&&(hi_Lt[ipix] == FALSE);idet++) {
	  pixvalue = l1a_ptr[ipix + npix * idet];
	  /*
	  pixdark = dark_rest[8*(recno-1)+idet];
	  */
	  pixdark = dark_mean[idet];
	  kneevalue = pcal_counts[idet][gain_ptr[idet]][FIRST_KNEE];
	  hi_Lt[ipix] = ((pixvalue - pixdark) >= kneevalue);
	}
      }
    }


    scan_no = crec;
    AS_pixels = stray_light;

    if (do_st == 0) {

      for (ipix=0; ipix < npix ; ipix++) 
          l2_flags_buffer[ipix] = 0;

    } else if (!recursive_flag) {

      do {
	if (stray_light_scan_no != scan_no) {

	  recursive_flag = 1;

	  get_l1a_rec(sd_id, stray_light_scan_no, cal_mod, 
			       l1a_dum, &st_l1b_data, &st_l2_flags);
	}
      }
      while (stray_light_corr(&initial,Ltyp_frac,Styp_frac,
			      nscan,npix,stray_light_scan_no++,dtype,
			      gain8,&pcal_rads[0][0][0],l1b_buffer,&sl_scan,
			      l2_flags_buffer,&AS_pixels) == NOTDONE);
      recursive_flag = 0;
      crec = scan_no;

    } 



    for (ipix=0; ipix < npix ; ipix++) {
        if (hi_Lt[ipix]) 
            l2_flags_buffer[ipix] |= MASK_HIGHLT1;
    } 

    *l1b_data = l1b_buffer;
    *l2_flags = l2_flags_buffer;
      
    return (0);
}


int openl1a_seawifs(filehandle *file)
{

    /*                                                                 */
    /* input files                                                     */
    /*                                                                 */
    char    *cal_path = file->calfile;

    /*                                                                 */
    /* get_l1a_open interface                                          */
    /*                                                                 */
    int32 fileID;
    int32 status;


    /*                                                                 */
    /* set static global evaluation switch                             */
    /*                                                                 */
    evalmask = file->input->evalmask;

    /* open the file and get the file ID */
    fileID = SDstart(file->name, DFACC_RDONLY);

    if (fileID < 0) {
      fprintf(stderr,
	      "-E %s Line %d: Error opening %s for reading.\n",
	      __FILE__,__LINE__,file->name);
      return(1);
    }

    status = getHDFattr(fileID, "Start Year", "", (VOIDP) &syear);
    status = getHDFattr(fileID, "Start Day", "", (VOIDP) &sday);
    status = getHDFattr(fileID, "Start Millisec", "", (VOIDP) &smsec);
    status = getHDFattr(fileID, "End Year", "", (VOIDP) &eyear);
    status = getHDFattr(fileID, "End Day", "", (VOIDP) &eday);
    status = getHDFattr(fileID, "End Millisec", "", (VOIDP) &emsec);
    status = getHDFattr(fileID, "Number of Scan Lines", "", (VOIDP) &nscan);
    status = getHDFattr(fileID, "Data Type", "", (VOIDP) &dtype);

    status = getHDFattr(fileID, "Pixels per Scan Line", "", (VOIDP) &npix);
    status = getHDFattr(fileID, "LAC Pixel Start Number", "", (VOIDP) &nsta);
    status = getHDFattr(fileID, "LAC Pixel Subsampling", "", (VOIDP) &ninc);

    status = getHDFattr(fileID, "Orbit Node Longitude", "", (VOIDP)&file->orbit_node_lon);
    status = getHDFattr(fileID, "Orbit Number", "", (VOIDP)&file->orbit_number);
    status = getHDFattr(fileID, "Node Crossing Time", "", (VOIDP)&file->node_crossing_time[0]);

    if (file->input != NULL) {
        Ltyp_frac   = file->input->sl_frac;
        stray_light = file->input->sl_pixl;
    } 

    if (stray_light == 0)
        do_st = 0;
    else {
        do_st = 1;
        if (stray_light < 0) {
            if (strcmp(dtype,"GAC") == 0) 
                stray_light = 4;
            else
                stray_light = 3;
        }
    }


    /* get the current calibration model solution */
    if (cal_path == NULL) {
       fprintf(stderr,
              "-E %s Line %d: No calibration file specified.\n",
              __FILE__,__LINE__);
        exit(1);
    }

    read_caltable(cal_path);

/*
        status = get_cal(cal_path, syear, sday, eday, smsec, npix, nsta, ninc,
		     dtype, tdi, &entry_year, &entry_day, &ref_year, 
		     &ref_day, &ref_minute, fp_temps, scan_mod, t_const,
                     t_linear_1, t_exponential_1, t_linear_2, t_exponential_2,
		     cal_offs, inst_tcorr, inst_tref, fp_tcorr, fp_tref, 
                     ms1_const, ms1_linear_1, ms1_exponential_1, ms1_linear_2, 
                     ms1_exponential_2, ms2_const, ms2_linear_1, 
                     ms2_exponential_1, ms2_linear_2, ms2_exponential_2, 
                     counts, rads);
*/

    if (status < 0) {
      fprintf(stderr,
          "-E- %s line %d: Error applying calibration table \"%s\".\n",
          __FILE__,__LINE__,cal_path);
      exit(status);
    }


    /* call cdata.f to initialize global FORTRAN common block data	*/
    cdata_();


    file->npix     = npix;
    file->nscan    = nscan;
    file->sensorID = SEAWIFS;
    file->sd_id = fileID;
    if (strcmp(dtype,"GAC") == 0)
        strcpy(file->spatialResolution, "4.5 km");
    else
        strcpy(file->spatialResolution,"1.1 km");

    strcpy(cal_path_tab, cal_path);

    l1a_data = (int16 *) calloc(npix*8*NREC_IN_BUF, sizeof(int16));
    l1a_back = (int16 *) calloc(npix*8*(STBUFSIZ-2), sizeof(int16));
    l1b_buffer = (float *) calloc(npix*8, sizeof(float));

    msec = (int32 *) calloc(nscan, sizeof(int32));
    status = rdSDS(file->sd_id,"msec",0,0,0,0, (VOIDP) msec);
    side = (int16 *) calloc(nscan, sizeof(int16));
    status = rdSDS(file->sd_id,"side",0,0,0,0, (VOIDP) side);

    dark_rest = (int16 *) calloc(nscan*8, sizeof(int16));
    status = rdSDS(file->sd_id,"dark_rest",0,0,0,0, (VOIDP) dark_rest);
    dark_rest_stat(dark_rest, nscan, dark_mean, dark_std);
    free(dark_rest);

    tilt = (float32 *) calloc(nscan, sizeof(float32));
    status = rdSDS(file->sd_id,"tilt",0,0,0,0, (VOIDP) tilt);

    spix = nsta - 1;
    dpix = ninc;
    
    return (status);
}


int readl1a_seawifs(filehandle *file, int32 recnum, l1str *l1rec)
{

    /*                                                                 */
    /* get_l1a_rec interface                                           */
    /*                                                                 */
    static cal_mod_struc cal_mod;  /* cal modification structure       */
    static int16   *l2_flags;      /* radiance quality flags for L2    */
    static float32 *l1b_data;

    int32   i;
    int32   status = 0;
    
    int32  nwave = l1rec->nbands;
    int32 *bindx = (int32*)l1rec->bindx;

    /*                                                                 */
    /* local vars                                                      */
    /*                                                                 */
    int32   ipix;                  /* pixel number                     */
    int32    iw, ib;

    int16 *l1a_dum = NULL;

    static int32  prev_recnum = -1;

    static int16 first = 1;
    static int32 ntilts;
    static int16 tilt_flags[20];
    static int16 tilt_ranges[2*20];

    int32  nflag[8];


    /* If reading out-of-sequence, force restart from beginning to get stlight */
    /* baf, 05-oct-2001.                                                       */
    if (recnum < prev_recnum) {
        printf("Reading out-of-sequence %d %d\n",recnum,prev_recnum);
        prev_recnum = -1;
    }


    /* Get tilt info */
    /* ------------- */
    if (first) {
      status = rdSDS(file->sd_id,"ntilts",0,0,0,0,(VOIDP) &ntilts); 
      status = rdSDS(file->sd_id,"tilt_flags",0,0,0,0,(VOIDP) tilt_flags); 
      status = rdSDS(file->sd_id,"tilt_ranges",0,0,0,0,(VOIDP) tilt_ranges); 
      first = 0;
    }


    /* Check for bad or changing tilt */
    /* ------------------------------ */
    for (i=0; i<ntilts; i++) {
      if (tilt_ranges[2*i] <= recnum+1 && tilt_ranges[2*i+1] >= recnum)
	break;
//      if (tilt_flags[i] == 0 ||    /* nadir */
//	  tilt_flags[i] == 1 ||    /* fwd */
//	  tilt_flags[i] == 2 )     /* aft */
//        bad_tilt = 0;              /* tilt result is ok */
//      else if (tilt_flags[i] == 3) /* tilt is changing */
//        bad_tilt = 1;
//      else                         /* tilt is unknown */
//        bad_tilt = 2;
    }


    /* Get nav flag */
    /* ------------ */
    status = rdSDS(file->sd_id,"nflag",recnum,0,1,8,(VOIDP) &nflag); 


    /* Get l1a data */
    /* ------------ */
    for (i=prev_recnum+1; i<=recnum; i++)
      status = get_l1a_rec(file->sd_id, i+1, &cal_mod, l1a_dum, 
			   &l1b_data, &l2_flags);

    /*                                                              */
    /* Copy scan geolocation and view geometry                      */
    /*                                                              */
    memcpy(l1rec->lat , ylat, npix*sizeof(float));
    memcpy(l1rec->lon , xlon, npix*sizeof(float));
    memcpy(l1rec->solz, solz, npix*sizeof(float));
    memcpy(l1rec->sola, sola, npix*sizeof(float));
    memcpy(l1rec->senz, senz, npix*sizeof(float));
    memcpy(l1rec->sena, sena, npix*sizeof(float));

    /*                                                              */
    /* Copy L1B radiances, pixel interlaced by band. Add per-band   */
    /* view angles.                                                 */
    /*                                                              */
    for (ipix = 0; ipix < file->npix; ipix++) {
      l1rec->pixnum[ipix] = spix + ipix*dpix;
      for (iw=0; iw<nwave; iw++) {
        ib = bindx[iw];
	l1rec->Lt   [ipix*nwave+ib] = l1b_data[iw*npix+ipix];
      }
      l1rec->stlight[ipix]  = ((l2_flags[ipix] & STRAYLIGHT) > 0);
      l1rec->hilt   [ipix]  = ((l2_flags[ipix] & HILT)    > 0);
      l1rec->navwarn[ipix]  = (nflag[7] & 1) | (nflag[0] & 1);
    }


    /*                                                              */
    /* Set scan time in L1B output record                           */
    /*                                                              */
    *l1rec->year = syear;
    *l1rec->day = sday;
    *l1rec->msec = msec[recnum];
    if (msec[recnum] < smsec) {                 /* adjust for day rollover */
        *l1rec->day += 1;
        if (*l1rec->day > (365 + (*l1rec->year % 4 == 0))) {
             *l1rec->year += 1;
             *l1rec->day   = 1;
        }
    }

    l1rec->tilt  = tilt[recnum];
    l1rec->mside = side[recnum];

    prev_recnum = recnum;

    return(status);
}


int readl1a_lonlat_seawifs(filehandle *file, int32 recnum, l1str *l1rec)
{
    int32   start[3];
    int32   edges[3];

    inputBuffer rdBuf;

    if (recnum > nscan) return 0;

    rdSDS(file->sd_id,"orb_vec",recnum,0,1,3,(VOIDP) rdBuf.orb_vec);
    rdSDS(file->sd_id,"sun_ref",recnum,0,1,3,(VOIDP) rdBuf.sun_ref);
    rdSDS(file->sd_id,"scan_ell",recnum,0,1,6,(VOIDP) rdBuf.scan_ell);

    start[0] = recnum;
    edges[0] = 1;
    start[1] = 0;
    edges[1] = 3;
    start[2] = 0;
    edges[2] = 3;
    SDreaddata(SDselect(file->sd_id, SDnametoindex(file->sd_id, "sen_mat")), 
                        start, NULL, edges, (VOIDP) rdBuf.sen_mat);

    geonav_lonlat_(rdBuf.orb_vec, rdBuf.sen_mat, rdBuf.scan_ell,
                   rdBuf.sun_ref,&nsta,&ninc,&npix,
                   l1rec->lat,l1rec->lon);

    return(0);
}


int closel1a_seawifs(filehandle *file)
{
  free(l1a_data);
  free(l1a_back);
  free(l1b_buffer);
  free(msec);
  free(side);
  free(tilt);

  SDend(file->sd_id);

  return(0);
}



