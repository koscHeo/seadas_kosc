
/*-----------------------------------------------------------------------------
    Function: getl3b_misc.c

    Contents:
	get_fsize	-  returns file size (SCENE / DAILY / OTHER)
	rdvdata		-  reads requested information from the given file
	rdmeta		-  reads metadata from the given file
	getindex	-  returns the index for accessing data of the 
				requested sbin
 	out_metadata	-  writes metadata to the output parameters
	out_data	-  writes requested data to the output data parameter 
	rd_data		-  reads requested data 
	attach_slave	-  accesses a requested data parameter
	attach_vdata	-  accesses a requested set of data
 	alloc_parm_bufs -  allocates buffer space for the data
	free_parm_bufs  -  frees the buffer space
   	init		-  initializes buf_sbin and buf_ebin arrarys 

    Other relevant files: 
	seabin.h	-  various #defined constants for level 3 binned
				data, also #includes hdf.h
	seaproto.h	-  prototypes for public level 3 output/input functions
	getl3b.h	-  prototypes for low-layer level 3 input functions
	getl3b.c	-  a higher layer of level 3 output functions

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development
        Lakshmi Kumar    Hughes STX      04/26/94    Modified to return    
					 	     packed array of data	
        Lakshmi Kumar    Hughes STX      07/01/94    Removed get_bins function
						     The no. of bins_infile is 							     obtained by reading the 
						     attribute "Data Bins"
        Lakshmi Kumar    Hughes STX      06/18/96    Changed defn. of MAX to
                                                     MAXVAL inorder to remove
                                                     compile time warning

        Joel Gales       Futuretech      11/24/99    Modify to handle variable
                                                     products.
-----------------------------------------------------------------------------*/

#include "seabin.h"
#include "seaproto.h"
#include "getl3b.h"

#define MAXNAMELEN  255

PRIVATE float32 *p[MAX_IN][NPARMS];
PRIVATE int32   *binno_list[MAX_IN]; 
PRIVATE int16   *nobs_list[MAX_IN], *nscenes_list[MAX_IN]; 
PRIVATE int16   *timerec_list[MAX_IN];
PRIVATE int8    *selcat_list[MAX_IN];
PRIVATE int32   *flags_list[MAX_IN];
PRIVATE float32 *weights_list[MAX_IN];

extern int32 NUMROWS;


/*----------------------------------------------------------------------------
    Function: get_fsize

    Returns: intn (file size)
        The return code is FAIL (-1) if an error occurs, file size 
	(SCENE/DAILY/OTHER) otherwise.

    Description:
        The function get_fsize returns the file size. 
	1 - if Scene, 
        2 - if Daily, and 
 	3 - for Other. 

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      char *       prod_type  I      input product type  

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development
        Lakshmi Kumar    Hughes STX      04/18/94    Converts prod_type string
						     to lower case in order to
						     check file type

-----------------------------------------------------------------------------*/
intn
#ifdef PROTOTYPE
  get_fsize(char *prod_type)
#else
  get_fsize(prod_type)
  char *prod_type;
#endif
{
  int32  i;
  char lprod_type[MAXVAL];

  for (i = 0; i < (int)strlen(prod_type); i++)   /* convert all chars to */
      lprod_type[i] = tolower(prod_type[i]);     /*   lower case         */
  lprod_type[i] = 0;

  if ((strcmp(lprod_type, "scene")) == 0)
     return SCENE;
  
  if ((strcmp(lprod_type, "daily")) == 0 || (strcmp(lprod_type, "day")) == 0)
     return DAILY;

  return OTHER;
}


/*-----------------------------------------------------------------------------
    Function: out_data
  
    Returns: intn (Status) 
      	Status returns as the number of data containing bins read or as
		-ve value if any error occurs.  

    Description:
        The function out_data outputs the bin numbers, number of observations, 
	time record, number of scenes, statistical weights, selection category,
	sums and sum squares of each bin for all the twelve parameters for 
	the requested nrec bins starting with sbin.

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      fst_call     I      flag to indicate if it is a first call
                                        or not for the given product file
      int32      fid          I      HDF file ID
      int32      sdfid        I      HDF ID 
      int32      prod_ID      I      ID number of the given product
      int32      parm_opt     I      option flag array of parameters to be read
      int32      sbin         I      start/first bin number to retrieve data
      int32      nrec         I      number of bins to retrieve
      int32      binlist_key  I      BinList access ID
      int32 *    begin        I      list of first bin in file for each row
      int32 *    extent       I      no. of bins written to file for each row
      int32 *    p_vskeys     O      parameter vdatas access IDs
      int32 *    last_bin     O      contains next bin no. to be processed
      int32 *    binno        O      bin numbers
      int16 *    nobs         O      number of observations of each bin
      int16 *    time_rec     O      time record for each bin
      int16 *    nscenes      O      no. of scenes contributing to each bin
      float32 *  weights      O      statistical weight for each bin
      int8  *    sel_cat      O      selection category for each binno bin
      int16 *    flags_set    O      flag field for each binno bin
      float32 *  l3b_data     O      sums and sumsquares of all the parameters

    Notes:  Status will be less than nrec if fewer than nrec data containing
 		bins were found starting with bin sbin and will be zero if no
		such bins were found.

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development
        Lakshmi Kumar    Hughes STX      04/26/94    Modified to return packed
					 	     array of data  
	Lakshmi Kumar 	 Hughes STX	 07/01/94    bins-infile value is
						     obtained by reading attr
						     "Data Bins"
	Lakshmi Kumar	 Hughes STX 	 02/13/95    added code to read sel_cat
						     data from BinList vdata
						     (ref. to I/O specs v4.2)
	Lakshmi Kumar	 HITC		 05/18/95    added code to read flags_
						     set from BinList vdata
						     (ref. to I/O specs V4.3)

        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure in
                                                     call to rd_data.
------------------------------------------------------------------------------*/

intn
#ifdef PROTOTYPE
  out_data(int32 fst_call, int32 fid, int32 sdfid, int32 prod_ID, l3b_prod *parm_opt, 
         int32 sbin, int32 nrec, int32 binlist_key, int32 *begin, int32 *extent, 
	 int32 *p_vskeys, int32 *last_bin, int32 *binno, int16 *nobs, 
	 int16 *time_rec, int16 *nscenes, float32 *weights, int8 *sel_cat,
	 int32 *flags_set, float32 *l3b_data)
#else
  out_data(fst_call, fid, sdfid, prod_ID, parm_opt, sbin, nrec, binlist_key, begin, 
		extent, p_vskeys, last_bin, binno, nobs, time_rec, nscenes, 
		weights, sel_cat, flags_set, l3b_data)
  int32 fst_call, fid, sdfid, prod_ID, sbin, nrec, binlist_key, *p_vskeys; 
  int32 *begin, *extent, *binno, *last_bin, *flags_set;
  int16 *nobs, *time_rec, *nscenes;
  int8  *sel_cat;
  float32 *weights, *l3b_data;
  l3b_prod *parm_opt;
#endif
{
  int32 	i, j, parm, k, done = 1, pr_extent = 0, loc_sbin = sbin;
  int32 	l3b_index, n_databins = 0, bins_infile, element_index, rd_flg;
  int32         sum_extent = 0, attrnum; 
  static int32 	start[MAX_IN], nelts[MAX_IN], eof[MAX_IN];

   attrnum = SDfindattr(sdfid, DATABINS);
   if ((SDreadattr(sdfid, attrnum, &bins_infile)) < 0)
        return FAIL;

   if ((sbin < 1 && *last_bin > 0) || sbin > 0) {
      if (sbin < 1)		/* set sbin to bin following the last bin */ 
         sbin = *last_bin;	/* retrieved by the last call             */
      pr_extent = 0;
      for(i = 0; i < NUMROWS && begin[i] < sbin; i++) {
        if (begin[i] != 0) {
           sum_extent = sum_extent + pr_extent;
           pr_extent = extent[i];
         }
       }
    }
   else { /* sbin < 1 and it is first call */
      for (i = 0; i < NUMROWS && begin[i] == 0; i++) ;
      sbin = begin[i];		/* set sbin to first bin in the product */
      sum_extent = 0;
    }
   loc_sbin = sbin;
   element_index = sum_extent;

   if (fst_call == 1) {
      rd_flg = 1;
      eof[prod_ID] = 0;
      if ((alloc_parm_bufs(prod_ID, BUFSZ)) == FAIL) 
         return FAIL; 
    }
   else {                        /* check if requested bin data in buffer */
     if (sbin >= binno_list[prod_ID][0] && 
		sbin<=binno_list[prod_ID][nelts[prod_ID]-1]){
        rd_flg = 0;
        element_index = get_index(sbin, prod_ID, nelts[prod_ID]);
      }
     else 
        rd_flg = 1;
    }
   parm = 0;
   done = 0;
   l3b_index = 0;
   j = 0;
   while (!done) {
      if(rd_flg) {
         eof[prod_ID] = 0;
         start[prod_ID] = element_index;
         nelts[prod_ID] = BUFSZ;
         if ((nelts[prod_ID] + start[prod_ID]) > bins_infile) {
            nelts[prod_ID] = bins_infile - start[prod_ID];
            eof[prod_ID] = 1;
          }
         if (nelts[prod_ID] > 0) { 
            if((rd_data(prod_ID, binlist_key, p_vskeys, start[prod_ID], 
		nelts[prod_ID], parm_opt))<0)
                return FAIL;
          }
         element_index = 0;
       }

      for(i = element_index; i<nelts[prod_ID] && 
		binno_list[prod_ID][i] < loc_sbin; i++) ;  

      for (k = i; i < nelts[prod_ID] && j < nrec; i++, j++) {
         n_databins++;
         binno[j] = binno_list[prod_ID][i];
         nobs[j] =  nobs_list[prod_ID][i];
         nscenes[j] = nscenes_list[prod_ID][i];
         time_rec[j] = timerec_list[prod_ID][i];
         weights[j] = weights_list[prod_ID][i];
	 sel_cat[j] = selcat_list[prod_ID][i];
	 flags_set[j] = flags_list[prod_ID][i];
         k = i * 2;
         for(parm = 0; parm < NPARMS; parm++) { 
            l3b_data[l3b_index++] = p[prod_ID][parm][k];
            l3b_data[l3b_index++] = p[prod_ID][parm][k+1];
          }
         loc_sbin = binno[j] + 1;
       }

      if (j < nrec && eof[prod_ID] != 1) {
         element_index = start[prod_ID] + nelts[prod_ID];
         done = 0;
         rd_flg = 1;
       }
      else
         done = 1;
    }
         
  *last_bin = binno[n_databins-1] + 1;

  return n_databins;
}

/*----------------------------------------------------------------------------
    Function: rd_data

    Returns: intn (Status)

    Description:
        The function rd_data reads the requested data from the appropriate 
	vdata

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      prod_ID      I      ID number of the given product
      int32      binlist_key  I      BinList vdata access ID  
      int32 *    p_vskey      I      parameter vdatas access IDs
      int32      start        I      start element number to read  
      int32      nelts        I      number of elements to read     
      int32 *    parm_opt     I      option flag array of parameters to be read


    Notes: The static variable *p[prod_ID][parm] is a global variable which
		points to a buffer to read the data in.

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development

        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure.
                                                     Replace parm_opt[i]
                                                     with parm_opt->code[i].
------------------------------------------------------------------------------*/

intn
#ifdef PROTOTYPE
  rd_data(int32 prod_ID, int32 binlist_key, int32 *p_vskey, int32 start, 
		int32 nelts, l3b_prod *parm_opt)
#else
  rd_data(prod_ID, binlist_key, p_vskey, start, nelts, parm_opt)
    int32 prod_ID, binlist_key, *p_vskey, start, nelts;
    l3b_prod *parm_opt;
#endif
{
  int32  i, ret;
  unsigned char *databuf;

   if((rdvdata(binlist_key, "bin_num", start, nelts,
                        (unsigned char *) binno_list[prod_ID])) < 0)
       return FAIL;

   if((rdvdata(binlist_key, "nobs", start, nelts, 
		(unsigned char *) nobs_list[prod_ID])) < 0)
	return FAIL;
   
   if((rdvdata(binlist_key, "nscenes",start,nelts,
		(unsigned char *)nscenes_list[prod_ID])) < 0)
	return FAIL;

   if((rdvdata(binlist_key,"time_rec",start,nelts,
		(unsigned char *)timerec_list[prod_ID])) < 0)
 	return FAIL;	

   if((rdvdata(binlist_key,"weights",start,nelts,
		(unsigned char *)weights_list[prod_ID])) < 0)
	return FAIL;
   
   if((rdvdata(binlist_key,"sel_cat", start,nelts,
		(unsigned char *)selcat_list[prod_ID])) < 0)
	return FAIL;
 
   if((rdvdata(binlist_key,"flags_set", start,nelts,
		(unsigned char *)flags_list[prod_ID])) < 0)
	return FAIL;
 
  for(i = 0; i < NPARMS; i++) {
    if (parm_opt->code[i] == 1)
    {
      if ((ret = VSseek(p_vskey[i], start)) < 0) {
        printf("Error: VSseek returned %d\n", ret);
        HEprint(stderr, 0);
        printf("***Error stack complete.\n");
        HEclear();
        return FAIL;
      }

      databuf = (unsigned char *)p[prod_ID][i];
      if ((ret = VSread(p_vskey[i], databuf, nelts, FULL_INTERLACE)) < 0) {
        printf("Error: VSread returned %d on parameter %d\n", ret, i);
        HEprint(stderr, 0);
        printf("***Error stack complete.\n");
        HEclear();
        return FAIL;
      }
    }
   }

 return SUCCEED; 
}

/*----------------------------------------------------------------------------
    Function: attach_slave

    Returns: intn (Status)

    Description:
        The function attach_slave attaches to the requested parameter vdata
	and returns the status.

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      fid          I      HDF file ID
      char *     sname        I      parameter name 
      int32      parm         I      parameter number
      int32 *    p_vskeys     O      parameter vdatas access IDs

    Notes: 

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development

------------------------------------------------------------------------------*/

intn
#ifdef PROTOTYPE
  attach_slave(int32 fid, char *sname, int32 parm, int32 *p_vskeys) 
#else
  attach_slave(fid, sname, parm, p_vskeys) 
    int32 fid, *p_vskeys, parm;
    char *sname;
#endif
{
   int32 vsid; 
   int32 nelements, interlace,  vsize;
   char fields[255],  name[255];

   if ((vsid = VSfind(fid, sname)) < 0)
        return FAIL;
   if ((p_vskeys[parm] = VSattach(fid, vsid, "r")) < 0)
        return FAIL;
   if ((VSinquire(p_vskeys[parm], &nelements, &interlace, fields,
                &vsize, name)) < 0)
        return FAIL;
   if ((VSsetfields(p_vskeys[parm], fields)) < 0)
        return FAIL;
   return SUCCEED; 
} 


/*-----------------------------------------------------------------------------
    Function: alloc_parm_bufs

    Returns: intn (Status)

    Description:
        Allocates buffers.  The parameter buffer size = nbins * 2 (sums and
	sum squares) * size of float.  The buffer pointers are global 
	variables to this file.

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      prod_ID      I      ID number of the product
      int32      nbins        I      number of bins 

    Notes: 

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development

------------------------------------------------------------------------------*/

intn
#ifdef PROTOTYPE
  alloc_parm_bufs(int32 prod_ID, int32 nbins)
#else
  alloc_parm_bufs(prod_ID, nbins)
  int32 prod_ID, nbins;
#endif
{
  int32 i;
  
  
   if((binno_list[prod_ID] = (int32 *) malloc(nbins * sizeof(int32)))==NULL)
          return FAIL;

   if ((nobs_list[prod_ID] = (int16 *) malloc(nbins * sizeof(int16))) == NULL)
         return FAIL;

   if((nscenes_list[prod_ID] = (int16 *) malloc(nbins * sizeof(int16))) == NULL)
        return FAIL;

   if((timerec_list[prod_ID] = (int16 *) malloc(nbins * sizeof(int16))) == NULL)
        return FAIL;

   if((weights_list[prod_ID] = (float32 *)malloc(nbins*sizeof(float32)))==NULL)
        return FAIL;

   if((selcat_list[prod_ID] = (int8 *)malloc(nbins*sizeof(int8)))==NULL)
        return FAIL;

   if((flags_list[prod_ID] = (int32 *)malloc(nbins*sizeof(int32)))==NULL)
        return FAIL;

  for (i = 0; i < NPARMS; i++) {
     if((p[prod_ID][i]=(float32 *)malloc(2 * nbins * sizeof(float32))) == NULL)
        return FAIL;

   }


  return SUCCEED;  
}

/*-----------------------------------------------------------------------------
    Function: free_parm_bufs

    Returns: intn (Status)

    Description:
        Frees data buffers.  

    Arguments: (in calling order)
      Type       Name        I/O     Description
      ----       ----        ---     -----------
      int32      prod_ID      I      ID number of the product

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar    Hughes STX      11/22/93    Original development

------------------------------------------------------------------------------*/

intn 
#ifdef PROTOTYPE
  free_parm_bufs(int32 prod_ID)
#else
  free_parm_bufs(prod_ID)
  int32 prod_ID;
#endif
{

  int32 i;

  free(binno_list[prod_ID]);
  free(nobs_list[prod_ID]);
  free(nscenes_list[prod_ID]);
  free(timerec_list[prod_ID]);
  free(weights_list[prod_ID]);
  free(selcat_list[prod_ID]);
  free(flags_list[prod_ID]);

  for(i = 0; i < NPARMS; i++)
    if (p[prod_ID][i] != NULL)
       free(p[prod_ID][i]);

  return SUCCEED;
} 



intn
#ifdef PROTOTYPE
  get_index(int32 sbin, int32 prod_ID, int32 bufsz)
#else
  get_index(sbin, prod_ID, bufsz)
  int32 sbin, prod_ID, bufsz;
#endif
{

  int32 i, done = 0, bufst, bufend, bufmid, element_index;
 
  bufst = 0;
  bufend = bufsz-1; 
  if (bufsz > 50) {
     while(!done) {
       bufmid = (bufend - bufst) /2;
       bufmid = bufst + bufmid; 
       if (sbin > binno_list[prod_ID][bufmid]) {
          bufst = bufmid;
          if((bufend - bufst) <= 50) {
             done = 1;
             element_index = bufst;
           }
        }
       else {
         if (sbin < binno_list[prod_ID][bufmid]) {
            bufend = bufmid;
            if((bufend - bufst) <= 50) {
               done = 1;
               element_index = bufst;
             }
          }
         else {
           done = 1;
           element_index = bufmid;
          }
        }
      }
   }
  if (bufsz <= 50)
     element_index = 0;

  for(i = element_index; i <= 50 && (binno_list[prod_ID][i] < sbin); i++) ;  

  return i;
}   
       
