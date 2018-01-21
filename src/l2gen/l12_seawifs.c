#include <strings.h>
#include "l1a.h"
#include "l1_octs_hdf.h"
#include "l12_proto.h"

#define   LAC_PIXEL_NUM           1285
#define   NREC_IN_BUF             10
#define   STBUFSIZ                5


static unsigned char genBuf[8192];

int32 seawifs_meta(filehandle *l1file, filehandle *outfile, 
		   int32 sscan, int32 escan, int32 dscan,
                   int32 spix, int32 dpix);


int32 l2_seawifs(filehandle *l1file, filehandle *l2file)
{
  PTB(seawifs_meta(l1file, l2file,
	       l1file->input->sline - 1, l1file->input->eline -1, 
	       l1file->input->dline, l1file->input->spixl - 1,
	       l1file->input->dpixl));

  return 0;
}




int32 l1b_seawifs(filehandle *l1file, filehandle *ofile, 
		  int32 sscan, int32 escan, int32 dline, 
                  int32 spixl, int32 epixl, int32 dpixl)
{
  int32 numBands = l1file->nbands;
  int32 i;
  int32 k;
  int32 n_read;
  int32 max_rec_in_rdbuf;
  int32 start[3] = {0,0,0};
  int32 edges[3];
  int32 stride[3] = {1,1,1};
  int32 sline = sscan+1;  
  int32 eline = escan+1;

  int32 nscan = (eline-sline)/dline + 1;

  int32 write_scan = 1;

  char *fld_name[] = {"s_satp","s_zerop","start_syn","stop_syn","dark_rest",
		      "sc_id","sc_ttag","sc_soh","inst_tlm","sc_ana","sc_dis",
		      "END"};
  int32 fld_2dim[] = {8,8,8,8,8,2,4,775,44,40,40};

  int32 dm[3];
  const char dm_name[3][80];

  if(seawifs_meta(l1file, ofile, sscan, escan, dline, spixl, dpixl)) {
      printf("-E- Can not write SeaWiFS metadata\n");
      exit(EXIT_FAILURE);
  }

  int32 l1a_sd_id = l1file->sd_id;

  idDS ds_id;
  ds_id.deflate = 0;
  ds_id.fid = ofile->sd_id;
  if(ofile->format == FMT_L1BNCDF)
      ds_id.fftype = DS_NCDF;
  else
      ds_id.fftype = DS_HDF;

  if ( ofile->format == FMT_L1BNCDF) {
    int dumdim;
    if ( nc_def_dim(ds_id.fid, "2", 2, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "4", 4, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "775", 775, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "44", 44, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "40", 40, &dumdim)
	 != NC_NOERR) exit(1);

    dm[0] = nscan;
    strcpy((char *) dm_name[0], "number_of_lines");
    dm[1] = numBands;
    strcpy((char *) dm_name[1], "number_of_reflective_bands");

  } else {

      dm[0] = nscan;
      strcpy((char *) dm_name[0], "Number of Scan Lines");
      dm[1] = numBands;
      strcpy((char *) dm_name[1], "band number");
  }

  PTB( createDS(ds_id, ofile->sensorID, "s_satp", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID, "s_zerop", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID, "start_syn", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID, "stop_syn", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID, "dark_rest", dm, dm_name));

  dm[1] = 2;
  strcpy((char *) dm_name[1], "2");
  PTB( createDS(ds_id, ofile->sensorID, "sc_id", dm, dm_name));

  dm[1] = 4;
  strcpy((char *) dm_name[1], "4");
  PTB( createDS(ds_id, ofile->sensorID, "sc_ttag", dm, dm_name));

  dm[1] = 775;
  strcpy((char *) dm_name[1], "775");
  PTB( createDS(ds_id, ofile->sensorID, "sc_soh", dm, dm_name));

  dm[1] = 44;
  strcpy((char *) dm_name[1], "44");
  PTB( createDS(ds_id, ofile->sensorID, "inst_tlm", dm, dm_name));

  dm[1] = 40;
  strcpy((char *) dm_name[1], "40");
  PTB( createDS(ds_id, ofile->sensorID, "sc_ana", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID, "sc_dis", dm, dm_name));

  /* Copy SDS data */
  max_rec_in_rdbuf = sline-1;

  for (i=sline; i<=eline; i++) {

    if (i >  max_rec_in_rdbuf) {

      n_read = (eline - max_rec_in_rdbuf) / dline;
      if ((eline - max_rec_in_rdbuf) % dline != 0) n_read++;
      n_read = (n_read >= NREC_IN_BUF) ? NREC_IN_BUF : n_read;

      k = 0;
      while (strcmp(fld_name[k],"END") != 0) {

	start[0] = i-1;
	stride[0] = dline;
	edges[0] = n_read;
	edges[1] = fld_2dim[k];
	SDreaddata(SDselect(l1a_sd_id, SDnametoindex(l1a_sd_id,fld_name[k])), 
		   start, stride, edges, (VOIDP) genBuf);
	
	start[0] = write_scan - 1;
	stride[0] = 1;
	PTB( writeDS(ds_id, fld_name[k], genBuf, 
		     start[0],start[1],start[2],edges[0],edges[1],edges[2]) );
	/*
	  SDwritedata(SDselect(sd_id, 
	  SDnametoindex(sd_id,fld_name[k])), 
	  start, stride, edges, (VOIDP) genBuf);
	  */
	k++;
	}
      
      max_rec_in_rdbuf += n_read * dline;
      write_scan += n_read;
    }  
  }

  return 0;
}



int32 seawifs_meta(filehandle *l1file, filehandle *ofile, 
		   int32 sscan, int32 escan, int32 dscan, 
		   int32 spix, int32 dpix)
{
  int32 numBands = l1file->nbands;
  int32 i;
  int32 k;
  int32 sd_id;
  int32 max_rec_in_rdbuf;
  int32 n_read;
  int32 start[3] = {0,0,0};
  int32 edges[3];
  int32 stride[3] = {1,1,1};

  int32 sline;
  int32 eline;
  int32 dline;
  int32 spixl;
  int32 dpixl;
  int32 nscan;

  int32 l1a_spixl;
  int32 l1a_dpixl;

  int32 write_scan = 1;

  int32 ntilts;
  int16 tilt_flags[20];
  int16 tilt_ranges[2*20];
  int16 *tf;

  char *fld_name[] = {"orb_vec","sun_ref","att_ang","scan_ell","nflag",
		      "sen_mat","END"};
  int32 fld_2dim[] = {3,3,3,6,8,3};

  char    buf1[32], buf2[32];
  int32_t year, day, hr, mn, sc;
  int16_t month, dom;

  extern float64         t_const[BANDS_DIMS_1A];
  extern float64         t_linear_1[BANDS_DIMS_1A];
  extern float64         t_exponential_1[BANDS_DIMS_1A];
  extern float64         t_linear_2[BANDS_DIMS_1A];
  extern float64         t_exponential_2[BANDS_DIMS_1A];
  extern float64         cal_offs[BANDS_DIMS_1A];
  extern float64         ms1_const[BANDS_DIMS_1A];
  extern float64         ms1_linear_1[BANDS_DIMS_1A];
  extern float64         ms1_exponential_1[BANDS_DIMS_1A];
  extern float64         ms1_linear_2[BANDS_DIMS_1A];
  extern float64         ms1_exponential_2[BANDS_DIMS_1A];
  extern float64         ms2_const[BANDS_DIMS_1A];
  extern float64         ms2_linear_1[BANDS_DIMS_1A];
  extern float64         ms2_exponential_1[BANDS_DIMS_1A];
  extern float64         ms2_linear_2[BANDS_DIMS_1A];
  extern float64         ms2_exponential_2[BANDS_DIMS_1A];

  extern int16           entry_year;
  extern int16           entry_day;
  extern int16           ref_year;
  extern int16           ref_day;
  extern int16           ref_minute;
  extern float32         counts[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];
  extern float32         rads[BANDS_DIMS_1A][GAINS_DIMS_1A][KNEES_DIMS_1A];

  int32 dm[3];
  const char dm_name[3][80];

  int32 l1a_sd_id = l1file->sd_id;

  sd_id = ofile->sd_id;
  nscan = (escan - sscan)/dscan + 1;

  sline = sscan + 1;
  eline = escan + 1;
  dline = dscan;
  
  idDS ds_id;
  ds_id.deflate = 0;
  ds_id.fid = sd_id;
  if(ofile->format == FMT_L2NCDF || ofile->format == FMT_L1BNCDF)
      ds_id.fftype = DS_NCDF;
  else
      ds_id.fftype = DS_HDF;

  if (ds_id.fftype == DS_NCDF) {
      int dumdim;
    if ( nc_def_dim(ds_id.fid, "vec", 3, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "ell", 6, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "nflag", 8, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "mtilt", 20, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "ltilt", 2, &dumdim)
	 != NC_NOERR) exit(1);
    if ( nc_def_dim(ds_id.fid, "unit", 1, &dumdim)
	 != NC_NOERR) exit(1);

    ds_id.fid = ofile->grp_id[4];
    strcpy((char *) dm_name[0], "number_of_lines");
  } else {
      strcpy((char *) dm_name[0], "Number of Scan Lines");
  }
  dm[0] = nscan;
 
  dm[1] = 3;
  strcpy((char *) dm_name[1], "vec");
  dm[2] = 3;
  strcpy((char *) dm_name[2], "vec");
  PTB( createDS(ds_id, ofile->sensorID, "orb_vec", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID, "sun_ref", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID ,"att_ang", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID ,"sen_mat", dm, dm_name));

  dm[1] = 6;
  strcpy((char *) dm_name[1], "ell");
  PTB( createDS(ds_id, ofile->sensorID,"scan_ell", dm, dm_name));

  dm[1] = 8;
  strcpy((char *) dm_name[1], "nflag");
  PTB( createDS(ds_id, ofile->sensorID, "nflag", dm, dm_name));

  if ( ofile->format == FMT_L2NCDF) ds_id.fid = ofile->grp_id[1];


  dm[0] = 20;
  strcpy((char *) dm_name[0], "mtilt");
  dm[1] = 2;
  strcpy((char *) dm_name[1], "ltilt");

  PTB( createDS(ds_id, ofile->sensorID, "tilt_flags", dm, dm_name));
  PTB( createDS(ds_id, ofile->sensorID, "tilt_ranges", dm, dm_name));

  dm[0] = 1;
  strcpy((char *) dm_name[0], "unit");
  PTB( createDS(ds_id, ofile->sensorID, "ntilts", dm, dm_name));


  /* Copy SDS data in "fld_name" array */
  /* --------------------------------- */
  max_rec_in_rdbuf = sline - 1;

  if ( ofile->format == FMT_L1HDF || ofile->format == FMT_L2HDF) {
    ds_id.fid = ofile->sd_id;
  } else if ( ofile->format == FMT_L2NCDF) {
    ds_id.fid = ofile->grp_id[4];
  }

  for (i=sline; i<=eline; i++) {

    if (i > max_rec_in_rdbuf) {

      n_read = (eline - max_rec_in_rdbuf) / dline;
      if ((eline - max_rec_in_rdbuf) % dline != 0) n_read++;
      n_read = (n_read >= NREC_IN_BUF) ? NREC_IN_BUF : n_read;

      k = 0;
      while (strcmp(fld_name[k],"END") != 0) {
	if (strcmp(fld_name[k],"sen_mat") == 0) edges[2] = 3;
	
	start[0] = i-1;
	stride[0] = dline;
	edges[0] = n_read;
	edges[1] = fld_2dim[k];
	PTB(SDreaddata(SDselect(l1a_sd_id, SDnametoindex(l1a_sd_id,fld_name[k])),
		   start, stride, edges, (VOIDP) genBuf));

	start[0] = write_scan - 1;
	stride[0] = 1;

	PTB( writeDS(ds_id, fld_name[k], genBuf, 
		     start[0],start[1],start[2],edges[0],edges[1],edges[2]) );

	  /*
	  status = SDwritedata(SDselect(sd_id, 
				       SDnametoindex(sd_id,fld_name[k])), 
			      start, stride, edges, (VOIDP) genBuf);
	  */
	k++;
      }

      max_rec_in_rdbuf += n_read * dline;
      write_scan += n_read;
    }
  }

  /* Write Nav Flag SDS */

  /* Write Tilt SDS */
  /* -------------- */
  rdSDS(l1a_sd_id,"ntilts",0,0,0,0,(VOIDP) &ntilts); 
  rdSDS(l1a_sd_id,"tilt_flags",0,0,0,0,(VOIDP) tilt_flags); 
  rdSDS(l1a_sd_id,"tilt_ranges",0,0,0,0,(VOIDP) tilt_ranges); 


  /* Handle subsetted and/or subsampled case */
  /* --------------------------------------- */
  if (sline != 1 || eline != l1file->nscan || dline != 1) {

    /* Build tilt flag (tf) array */
    tf = (short int *) calloc(l1file->nscan, sizeof(short int));
    for (k=0; k<ntilts; k++) 
      for (i=tilt_ranges[k*2]-1; i<tilt_ranges[k*2+1]; i++) 
	tf[i] = tilt_flags[k]; 

    /* Loop through tf array and build tilt_flags & tilt_ranges */
    k = 0;
    ntilts = 1;
    tilt_flags[0] = tf[sline-1];
    for (i=sline-1; i<eline; i+=dline) {
      k++;
      if (tf[i] != tilt_flags[ntilts-1]) {
	tilt_flags[ntilts] = tf[i];

	tilt_ranges[2*ntilts-1] = k-1;
	tilt_ranges[2*ntilts] = k;

	ntilts++;
      }
    }
    tilt_ranges[2*ntilts-1] = k;
    for (i=2*ntilts; i<40; i++) tilt_ranges[i] = 0;
    for (i=ntilts; i<20; i++) tilt_flags[i] = 0;

    free(tf);
  }

  if ( ofile->format == FMT_L1HDF || ofile->format == FMT_L2HDF) {
    ds_id.fid = ofile->sd_id;
  } else if ( ofile->format == FMT_L2NCDF) {
    ds_id.fid = ofile->grp_id[1];
  }

  PTB( writeDS(ds_id, "ntilts", &ntilts, 0, 0, 0, 1, 0, 0) );
  PTB( writeDS(ds_id, "tilt_flags", tilt_flags, 0, 0, 0, 20, 0, 0) );
  PTB( writeDS(ds_id, "tilt_ranges", tilt_ranges, 0, 0, 0, 20, 2, 0) );

  /* Write Global Attributes */
  ds_id.fid = ofile->sd_id;
  if ( ofile->format == FMT_L1HDF || ofile->format == FMT_L2HDF) {
    ds_id.sid = ds_id.fid;
  } else if ( ofile->format == FMT_L2NCDF) {
    ds_id.sid = NC_GLOBAL;
  }

  if (ds_id.fftype == DS_NCDF) {
      int32 tmp32;
      getHDFattr(l1a_sd_id, "Orbit Number", "", (VOIDP) &tmp32);
      PTB(SetI32GA(ds_id, "orbit_number", tmp32));

      bzero(genBuf, sizeof(genBuf));
      getHDFattr(l1a_sd_id, "Data Center", "", (VOIDP) &genBuf);
      PTB(SetChrGA( ds_id, "project", (const char*) genBuf) );

      bzero(genBuf, sizeof(genBuf));
      getHDFattr(l1a_sd_id, "Mission", "", (VOIDP) &genBuf);
      PTB(SetChrGA( ds_id, "platform", (const char*) genBuf) );

      float tmpFloat;
      getHDFattr(l1a_sd_id, "Orbit Node Longitude", "", (VOIDP) &tmpFloat);
      PTB(SetF32GA( ds_id, "equatorCrossingLongitude", tmpFloat) );

      bzero(genBuf, sizeof(genBuf));
      getHDFattr(l1a_sd_id, "Data Type", "", (VOIDP) &genBuf);
      PTB(SetChrGA( ds_id, "Data Type", (const char*) genBuf) );

      getHDFattr(l1a_sd_id, "LAC Pixel Start Number", "", (VOIDP) &l1a_spixl);
      getHDFattr(l1a_sd_id, "LAC Pixel Subsampling", "",  (VOIDP) &l1a_dpixl);

      spixl = (l1a_dpixl * spix) + l1a_spixl;
      dpixl = l1a_dpixl * dpix;

      PTB(SetI32GA( ds_id, "LAC Pixel Start Number", spixl) );
      PTB(SetI32GA( ds_id, "LAC Pixel Subsampling", dpixl) );

      bzero(genBuf, sizeof(genBuf));
      getHDFattr(l1a_sd_id, "Node Crossing Time", "", (VOIDP) &genBuf);
      strncpy( buf1, (const char *) genBuf, 4);
      buf1[4] = 0;
      year = atol(buf1);
      strncpy( buf1, (const char *) &genBuf[4], 3);
      buf1[3] = 0;
      day = atol(buf1);
      yd2md(year, day, &month, &dom);

      strncpy( buf2, (const char *) &genBuf[7], 2);
      buf2[2] = 0;
      hr = atol(buf2);
      strncpy( buf2, (const char *) &genBuf[9], 2);
      buf2[2] = 0;
      mn = atol(buf2);
      strncpy( buf2, (const char *) &genBuf[11], 2);
      buf2[2] = 0;
      sc = atol(buf2);
      sprintf(buf1, "%d-%02d-%02dT%02d:%02d:%02dZ", year, month, dom, hr, mn, sc);

      PTB(SetChrGA( ds_id, "equatorCrossingDateTime", (const char*) buf1) );

  } else {
      PTB(getHDFattr(l1a_sd_id, "Orbit Number", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Orbit Number", DFNT_INT32, 1, genBuf) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Data Center", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Data Center", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP) genBuf) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Station Name", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Station Name", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP)  genBuf) );

      PTB(getHDFattr(l1a_sd_id, "Station Latitude", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Station Latitude", DFNT_FLOAT32, 1, (VOIDP) genBuf) );

      PTB(getHDFattr(l1a_sd_id, "Station Longitude", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Station Longitude", DFNT_FLOAT32, 1, (VOIDP) genBuf) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Mission", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Mission", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP) genBuf) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Mission Characteristics", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Mission Characteristics", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP) genBuf) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Sensor", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Sensor", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP) genBuf) );

      PTB(getHDFattr(l1a_sd_id, "Orbit Node Longitude", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Orbit Node Longitude", DFNT_FLOAT32, 1, (VOIDP) genBuf) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Sensor Characteristics", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Sensor Characteristics", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP) genBuf) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Data Type", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Data Type", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP) genBuf) );

      PTB(getHDFattr(l1a_sd_id, "LAC Pixel Start Number", "", (VOIDP) &l1a_spixl));
      PTB(getHDFattr(l1a_sd_id, "LAC Pixel Subsampling", "",  (VOIDP) &l1a_dpixl));

      spixl = (l1a_dpixl * spix) + l1a_spixl;
      dpixl = l1a_dpixl * dpix;

      PTB(sd_setattr( sd_id, "LAC Pixel Start Number", DFNT_INT32, 1, (VOIDP) &spixl) );
      PTB(sd_setattr( sd_id, "LAC Pixel Subsampling", DFNT_INT32, 1,  (VOIDP) &dpixl) );

      bzero(genBuf, sizeof(genBuf));
      PTB(getHDFattr(l1a_sd_id, "Node Crossing Time", "", (VOIDP) &genBuf));
      PTB(sd_setattr( sd_id, "Node Crossing Time", DFNT_CHAR, strlen((const char *) genBuf)+1, (VOIDP) genBuf) );
  }

  return 0;
}



