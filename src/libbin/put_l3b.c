/*------------------------------------------------------------------------------
    File:  put_l3b.c

    Contents:
	put_l3b_open   - initializes file and data structures; writes out
			 most metadata values
	put_l3b_record - writes out any number of data "records" to the file
	put_l3b_close  - writes out index for grid; writes out some final
			 metadata values; performs clean-up chores

    Notes:
	o These routines must be called in order as listed above.
	o put_l3b_record may be called multiple times.
	o Files created by these functions follow the guidelines set forth
	  in EOSDIS' "Format System Implementation Guidelines," section 5.6.
	  The SEAGrid scheme is used.

    Other relevant files:
	seabin.h    - various #defined constants for level 3 binned data;
		      also #includes hdf.h and seaproto.h
	seaproto.h  - prototypes for public level 3 output functions
	seaprotoi.h - prototypes for low-layer level 3 output functions
	l3b_misc.c  - a lower layer of level 3 output functions

    NOTE:  The minimum file sizes of daily, scene, and other files are defined
	   in seabin.h file.

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
	Lakshmi Kumar    Hughes STX      02/10/94    Made changes in 
						     put_l3b_close to add some
						     more global attributes
        Doug Ilg         Hughes STX      06/30/94    Updated to reflect I/O
						     Specs v3.2 changes 
	Lakshmi Kumar	 Hughes STX	 09/01/94    Removed defn. of TESTCODE
						     and the main function
						     which was defined before
						     for creating test files.
						     Made changes to get soft-
						     ware ID info from Makefile
	Lakshmi Kumar	 Hughes STX 	 11/17/94    Removed station from 
						     put_l3b_close call 
						     (I/O Specs v4.0) 
	Lakshmi Kumar	 Hughes STX	 03/01/95    Added few global attrs &
						     a field "sel_cat" in 
						     BinList vdata(I/O specs
						     V4.2, product specs V2.6)
						     Fixed a bug and changed
						     slvid var name to slvidp
						     Changed epsilon name to
						     eps_78
	Lakshmi Kumar	 HITC		 05/18/95    Added flags_set field to
						     BinList vdata and start_
						     num fld to BinIndex vdata
						     Changed datatype of attr. 
						     "L2 Flag Usage" 
	Lakshmi Kumar    HSTX	 	 09/25/95    Added "First Orbit" &
						     corner coordinates as 
						     global attributes (V4.4 
						     I/O & 2.8 product specs)
	Lakshmi Kumar	 HSTX		 10/27/95    Changed "First Orbit" to
						     "Orbit" & added "Start 
						     Orbit" & "End Orbit" attrs.
	Lakshmi Kumar	 HSTX		 12/21/95    Added global attribute
						     fstcall to indicate put_
						     l3b_record is called for
						     the first time for that
						     file id
	Lakshmi Kumar    HSTX		02/29/96     oops...reinitialization of
						     lastrec variable was done
						     before it was used to 
						     findout nlat	
        Lakshmi Kumar    Hughes STX      06/18/96    Changed defn. of MAX to
                                                     MAXVAL inorder to remove
                                                     compile time warning
	Lakshmi Kumar	 Hughes STX	 07/25/97    As a request from DAAC,
						     this code has been updated
						     for creating main file of
						     multifiles with the ext.
						     ".main"
        Joel Gales       Futuretech      11/24/99    Modify to handle variable
                                                     products.
        Joel Gales       Futuretech      03/21/01    Increase start_num array
                                                     size by one.
----------------------------------------------------------------------------*/
#include <math.h>
#include "seabin.h"
#include "seaprotoi.h"

/*  Values for "continuity" - Used to ensure  **
**  proper ordering of function calls         */
#define CLEAR   0
#define INIT    1
#define RECORDS 2

#define max(a,b)               (a<b ? b : a)
#define min(a,b)               (a>b ? b : a)

/*  Private Global Variables  */
PRIVATE int32 sdfid[MAX_OUT], fid[MAX_OUT], geomid[MAX_OUT], ndxid[MAX_OUT];
PRIVATE int32 gridid[MAX_OUT], mstrid[MAX_OUT], *slvidp[MAX_OUT];
PRIVATE int32 continuity[MAX_OUT];
PRIVATE int32 lastrec[MAX_OUT];
PRIVATE int32 fstcall[MAX_OUT];
PRIVATE int32 row[MAX_OUT], rowcount[MAX_OUT];

PRIVATE int32 *begin[MAX_OUT];
PRIVATE int32 *extent[MAX_OUT];
PRIVATE int32 *lastdatabin[MAX_OUT];

PRIVATE int32 start_num[2160*16+1];
PRIVATE int32 lastbin[2160*16];
PRIVATE int32 *numbin;

extern  int32 NUMROWS;
extern  int32 TOTBINS;

/*, numbin[NUMROWS+1];*/
/*PRIVATE float64 latbin[NUMROWS+1];*/



/*  Private function  */
PRIVATE void init(void);


/*------------------------------------------------------------------------------
    Function: put_l3b_open

    Returns: int32 (ID number of product file)
	The return code is prod_ID, the ID number for the product file
	being written, or FAIL (-1) if an error occurs.  The value of
	prod_ID must be sent to successive calls to put_l3b_record and
	put_l3b_close to identify the product file to use.

    Description:
        The function put_l3b_open creates a level 3 binned HDF file with
	the file name given by l3b_path.  Level 3 binned file metadata
	will also be written.  See "SEAWIFS HDF INPUT/OUTPUT SUBPROGRAM
	INTERFACES" and "SeaWiFS Level 3 Binned Data CDL" (both by Fred
	Patt) for details.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      char * 	l3b_path    I     path for level 1A file
      int32	fsize	    I	  product type (SCENE, DAILY,...)
      char *	prod_type   I	  binning period descr. ("Scene", "Daily",...)
      char *   	ptime       I     data processing time (ASCII)
      char *	proc_con    I	  input to processing routines
      int32 *   parm_opt    I     option array of parameters to be written

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
        Doug Ilg         Hughes STX      06/30/94    Modified some of the 
						     attribute names(product
							specs)	
	Lakshmi Kumar	 Hughes STX	 09/01/94    Software ID is provided 
					 	     through the Makefile
	Lakshmi Kumar	 Hughes STX	 11/09/94    Added code to remove dir
						     path from the product name
 	Lakshmi Kumar	 Hughes STX	 09/25/95    Added argument "orbit"
						     Added code to write orbit
						     as global attribute
        Lakshmi Kumar    Hughes STX      07/25/97    As a request from DAAC,
                                                     this code has been updated
                                                     for creating main file of
                                                     multifiles with the ext.
                                                     ".main"
        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure.
                                                     Replace parm_opt[i] with
                                                     parm_opt->code[i].
        Joel Gales       Futuretech      03/29/01    Write some attributes
                                                     from metadata values rather
                                                     than from hard-coded values.
------------------------------------------------------------------------------*/
int32
put_l3b_open(char *l3b_path, char *replaces, int32 fsize, char *prod_type, 
	     char *ptime, int32 orbit, int32 start_orb, int32 end_orb,
	     char *proc_con, char *soft_name, char *soft_ver, 
	     char *input_parms, l3b_prod *parm_opt, meta_l3bType *meta_l3b)
{

    int32 prod_ID;
    int32 i, c;
    static int first = 1;
    /*
    char *softwareName = "BIN"; 			       	
    char *softwareVer = "5.1";					
    */
    char *loc_l3b_path;						/* LK */
    char loc_replace[MAXVAL];					/* LK */
    char mainfile[2048];					/* LK */
    char loc_prod_type[30];


    /*  init fids to "unused"  */
    if(first == 1) {
	init();
	first = 0;
    }

    /*  get a prod_ID  */
    for (i=0; i<MAX_OUT; i++)
	if (fid[i] == -1) break;
    prod_ID = i;

    /*  verify prod_ID  */
    if (prod_ID >= MAX_OUT) {
	fprintf(stderr, "*** Error:  attempt to open too many files.\n");
	return FAIL;
    }

    /*  check for proper calling order  */
    if (continuity[prod_ID] != CLEAR) {
	fprintf(stderr, "*** Error:  put_l3b_open() called out of order.\n");
	return FAIL;
    }

    /*  open file; fire up interfaces  */
    /* copy the file path to local variable mainfile and add '.main' to it. */
	
    strcpy(mainfile, l3b_path);

    for (c = 0; c < strlen(prod_type); c++)
	loc_prod_type[c] = tolower(prod_type[c]);
    loc_prod_type[c] = (char) 0;

    if ((strcmp(loc_prod_type, "scene")) != 0)
        strcat(mainfile, ".main");

    if ((sdfid[prod_ID] = SDstart(mainfile, DFACC_CREATE)) < 0) {
	printf("\nput_l3b_open: Error creating output file - %s\n", mainfile);
        return FAIL;
     }

    SDsetblocksize(sdfid[prod_ID], 4096*10);
    printf("setting sd blocksize\n");

    if ((fid[prod_ID] = Hopen(mainfile, DFACC_RDWR, 0)) < 0) {
	printf("\nput_l3b_open: Error opening output file - %s\n", mainfile);
        return FAIL;
     }
    Vstart(fid[prod_ID]);

#ifdef DEBUG
printf("put_l3b_open/fid = %d\n", fid[prod_ID]);
#endif

    begin[prod_ID] = calloc(NUMROWS+1, sizeof(int32));
    extent[prod_ID] = calloc(NUMROWS+1, sizeof(int32));
    lastdatabin[prod_ID] = calloc(NUMROWS+1, sizeof(int32));

    /*  initialize index  */
    for (i=0; i<NUMROWS; i++) {
	begin[prod_ID][i] = 0;
	extent[prod_ID][i] = 0;
    }
    rowcount[prod_ID] = 0;
    row[prod_ID] = 0;

    /*  write out global attributes  */

    /*  Remove direcotry path from product name */		/* LK */
    if ((loc_l3b_path = strrchr(mainfile, '/')) != NULL)	/* LK */
        loc_l3b_path++;						/* LK */
    else							/* LK */
        loc_l3b_path = mainfile;				/* LK */

    SDsetattr(sdfid[prod_ID], "Product Name", DFNT_CHAR, 
	      strlen(loc_l3b_path) + 1, (VOIDP)loc_l3b_path); 	/* LK */


    sprintf(meta_l3b->title, "%s%s", meta_l3b->sensor_name, " Level-3 Binned Data"); 
    SDsetattr(sdfid[prod_ID], "title", DFNT_CHAR, 
	      strlen(meta_l3b->sensor_name) + strlen(" Level-3 Binned Data") + 1, 
	      meta_l3b->title); 

    SDsetattr(sdfid[prod_ID], SENSOR_NAME, DFNT_CHAR, strlen(meta_l3b->sensor_name) + 1, 
	      meta_l3b->sensor_name);


    SDsetattr(sdfid[prod_ID], "project", DFNT_CHAR, 
	      strlen(DCENTER_VAL) + 1, DCENTER_VAL);



    if (meta_l3b->mission != 0x0) {
      SDsetattr(sdfid[prod_ID], "mission", DFNT_CHAR, strlen(meta_l3b->mission) + 1, 
		(VOIDP)meta_l3b->mission);
    }

    SDsetattr(sdfid[prod_ID], "Product Type", DFNT_CHAR,
	      strlen(prod_type) + 1, (VOIDP)prod_type);

    if (replaces == NULL || (strcmp(replaces, "")) == 0) {       /* LK */ 	
       strcpy(loc_replace, "ORIGINAL");
       replaces = loc_replace;
     }
    SDsetattr(sdfid[prod_ID], "Replacement Flag", DFNT_CHAR,
	      strlen(replaces) + 1, (VOIDP)replaces);	        /* LK */
    SDsetattr(sdfid[prod_ID], "Software Name", DFNT_CHAR,
              strlen(soft_name) + 1, soft_name);
    SDsetattr(sdfid[prod_ID], "Software Version", DFNT_CHAR,
              strlen(soft_ver) + 1, soft_ver);
								/* LK */
    SDsetattr(sdfid[prod_ID], "Orbit", DFNT_INT32, 1, (VOIDP)&orbit);

    SDsetattr(sdfid[prod_ID], "Start Orbit", DFNT_INT32, 1, (VOIDP)&start_orb);

    SDsetattr(sdfid[prod_ID], "End Orbit", DFNT_INT32, 1, (VOIDP)&end_orb); 

    SDsetattr(sdfid[prod_ID], "Processing Time", DFNT_CHAR, strlen(ptime) + 1,
              (VOIDP)ptime);
    SDsetattr(sdfid[prod_ID], "Processing Control", DFNT_CHAR,
              strlen(proc_con) + 1, (VOIDP)proc_con);

    SDsetattr(sdfid[prod_ID], "Input Parameters", DFNT_CHAR,
              strlen(input_parms) + 1, (VOIDP)input_parms);

    /*  create grid structure (Vgroup)  */
    gridid[prod_ID] = setupgrid(prod_ID, fid[prod_ID]);

    /*  write Geometry Vdata and add to grid Vgroup  */
    geomid[prod_ID] = writegeom(fid[prod_ID], (int32) 2*NUMROWS);
    Vaddtagref(gridid[prod_ID], DFTAG_VH, geomid[prod_ID]);
/*
    Vaddtagref(gridid[prod_ID], DFTAG_VH, VSQueryref(geomid[prod_ID]));
*/
    /*  init DataMaster Vdata and add to grid Vgroup  */
    mstrid[prod_ID] = setupmaster(fid[prod_ID], fsize);
    Vaddtagref(gridid[prod_ID], DFTAG_VH, VSQueryref(mstrid[prod_ID]));

    /*  init DataSlave Vdatas and add to grid Vgroup  */
    slvidp[prod_ID] = setupslaves(fid[prod_ID], fsize, l3b_path, parm_opt);
    for (i=0; i<NPARMS; i++)
      if (parm_opt->code[i] == 1)
        Vaddtagref(gridid[prod_ID], DFTAG_VH, VSQueryref(slvidp[prod_ID][i]));

    continuity[prod_ID] = INIT;
    return prod_ID;
}


/*------------------------------------------------------------------------------
    Function: put_l3b_record

    Returns: intn (status)
	The return code is FAIL (-1) if an error occurs, SUCCEED (0)
	otherwise.

    Description:
        The function put_l3b_record stores nrec "records" of binned data
	into the HDF file.  The data are buffered for improved performance.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      int32        prod_ID   I       id of file to write to
      int32        nrec      I       number or bins to write
      int32 *      binno     I       bin number for each bin
      int16 *      nobs      I       # of observations contributing to each bin
      int16 *      time_rec  I       time record for each bin
      int16 *      nscenes   I       # of scene contributing to each bin
      float32 *    weights   I	     statistical weight for each bin
      int8  *  	   sel_cat   I       selection category for eh binno bin
      float32 *    l3b_data  I	     3-D array of sums and sum-squares
      int32 *      parm_opt  I       option array of parameters to be written 
    
    Notes:
	o Except for prod_ID, nrec, and l3b_data, all arguments are pointers
	  to one-dimensional arrays of values (one per bin).
	o l3b_data must be stored in such a way as to make it equivalent to
	  the [row-major] C-language declaration
	                  float32 l3b_data[12][nrec][2];
	  Where 12 is the number of arguments in the data set, nrec is the
	  number of records being written by this call, and 2 is the number
	  of statistics stored for each argument (i.e., sum and sum-squared).

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
        Doug Ilg         Hughes STX      06/30/94    Changed time_rec from
						     4byte int to 2byte int
	Lakshmi Kumar    Hughes STX	 02/13/95    Added argument sel_cat
						     (ref. I/O specs v4.2)
	Lakshmi Kumar	 Hughes STX	 09/25/95    Added code to save last
						     databins in an array, wh
						     will be used in figuring 
						     out easternmost longitude
	Lakshmi Kumar	 Hughes STX	 12/29/95    Added fstcall argument to
						     buffbins
        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure
                                                     in call to buffbins.
------------------------------------------------------------------------------*/
intn
put_l3b_record(int32 prod_ID, int32 nrec, int32 *binno, int16 *nobs,
		int16 *time_rec, int16 *nscenes, float32 *weights, 
		int8 *sel_cat, int32 *flags_set, float32 *l3b_data, 
                l3b_prod *parm_opt)
{
    register int i;


    /*  verify prod_ID  */
    if (prod_ID >= MAX_OUT) {
	fprintf(stderr, "*** Error:  bad file id sent to put_l3b_record().\n");
	return FAIL;
    }
    if (fid[prod_ID] == -1) {
	fprintf(stderr, "*** Error:  bad file id sent to put_l3b_record().\n");
	return FAIL;
    }

    /*  check for proper calling order  */
    if (continuity[prod_ID] == CLEAR) {
	fprintf(stderr, "*** Error:  put_l3b_record() called out of order.\n");
	return FAIL;  /* exit(1) instead? */
    }

    /*  update index  */
    for (i=0; i<nrec; i++) {
	if ((binno[i] > TOTBINS) || (binno[i] < 1)) {
	    fprintf(stderr, "*** Warning:  undefined bin number.  ");
	    fprintf(stderr, "Bin %d skipped on prod_ID %d\n", binno[i],prod_ID);
	}
        if (binno[i] < lastrec[prod_ID]) {
	   fprintf(stderr, "*** Error:  Bin numbers not in ascending order.\n");
	   fprintf(stderr, "*** calling bin#%d after bin#%d\n", binno[i],
			lastrec[prod_ID]);
	   fprintf(stderr, "*** Warning: BinIndex may get corrupted. \n");
         }
		
	if (binno[i] > lastbin[row[prod_ID]]) { /* next row */
	    extent[prod_ID][row[prod_ID]] = rowcount[prod_ID];
	    while (binno[i] > lastbin[row[prod_ID]])
	        row[prod_ID]++;
	    rowcount[prod_ID] = 1;
            lastdatabin[prod_ID][row[prod_ID]] = binno[i-1];	/* LK */ 
            if (i == 0)
		lastdatabin[prod_ID][row[prod_ID]] = lastrec[prod_ID];
	    begin[prod_ID][row[prod_ID]] = binno[i];
	} else { /* add bin to rowcount */
	    rowcount[prod_ID]++;
	}
    }
    lastrec[prod_ID] = binno[i-1];				/* LK */

    buffbins(fstcall, prod_ID, mstrid[prod_ID], slvidp[prod_ID], nrec, binno, 
	     nobs, time_rec, nscenes, weights, sel_cat, flags_set, l3b_data,
             parm_opt);
    continuity[prod_ID] = RECORDS;
    return SUCCEED;
}


/*------------------------------------------------------------------------------
    Function: put_l3b_close

    Returns: intn (status)
	The return code is FAIL (-1) if an error occurs, SUCCEED (0)
	otherwise.

    Description:
        The function put_l3b_close writes out some final metadata, shuts
	down the HDF interfaces, and closes the HDF file.

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      int32        prod_ID   I       id of file to close
      int16        bin_syear I       start year for binning period
      int16        bin_sday  I       start day-of-year for binning period
      int16        bin_eyear I       end year for binning period
      int16        bin_eday  I       end day-of-year for binning period
      int16        syear     I       data start time(year)
      int16        sday      I       data start time(day)
      int32        smsec     I       data start time(msec)
      int16        eyear     I       data end time(year)
      int16        eday      I       data end time(day)
      int32        emsec     I       data end time(msec)
      char *	   infiles   I	     list of input files
      char *	   flag_names I      algorithm names set for bits in l2_flags
      char *       flag_use  I	     algorithm names for bits in l2_flags       
      uint8 *      eng_q_use I       values for exclusion during binning
      char *	   proc_log  I       output of processing routines
      int32 *      parm_opt  I       option array of parameters to be written 
    
    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
	Lakshmi Kumar    Hughes STX      02/10/94    added bin_syear,bin_sday,
						     bin_eyear and bin_eday
        Doug Ilg         Hughes STX      06/30/94    Added "data bins" attribute
	Lakshmi Kumar	 Hughes STX	 11/17/94    Removed station from the
						     call (I/O Specs v4.0)
	Lakshmi Kumar	 Hughes STX	 02/21/95    Added global attributes
						     flag_use and eng_q_use
						     (ref. v4.2 I/O specs)
						     Fixed a bug.
	Lakshmi Kumar	 Hughes STX	 05/09/95    Changed datatype of global
						     attribute "flag_use"
						     (ref. V4.3 I/O specs)
	Lakshmi Kumar	 HITC		 05/25/95    Added global attribute 
						     flag_names (V4.3 I/O specs)
	Lakshmi Kumar	 Hughes STX	 09/25/95    Added a function call to
						     get_coords to get nlat,
						     slat, elon & wlon
						     Added code to write north-
						     ern/southernmost lats &
						     eastern/westernmost lons
						     as global attributes
	Lakshmi Kumar	 Hughes STX	 12/21/95    Reinitialize lastrec &
						     fstcall of this freed 
						     product id
	Lakshmi Kumar	 Hughes STX	 02/29/96    reinitialized lastrec 
						     before its use.  Moved it
						     after the call to 
						     get_coords

        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure in
                                                     call to closedata.

        Joel Gales       Futuretech      03/14/00    Write "Units" attribute.
        Joel Gales       Futuretech      11/09/01    Fix "End Time" roll-over
	                                             problem
------------------------------------------------------------------------------*/
intn
put_l3b_close(int32 prod_ID, int16 bin_syear, int16 bin_sday, int16 bin_eyear, 
		int16 bin_eday, int16 syear, int16 sday, int32 smsec,
	      	int16 eyear, int16 eday, int32 emsec, char *infiles,
	      	char *flag_names, char *flag_use, uint8 *eng_q_use, 
		l3b_prod *parm_opt)
{
    int32 total;
    float pct;
    char string[20], *loc_infiles;				/* LK */
    div_t quot1, quot2, quot3;
    float32  stlat = 37.9272;
    float32  stlon = -75.4753;
    float32  slat, nlat, elon, wlon;				/* LK */


    printf("in close\n");

    /*  verify prod_ID  */
    if (prod_ID >= MAX_OUT) {
	fprintf(stderr, "*** Error:  bad file id sent to put_l3b_close().\n");
	return FAIL;
    }
    if (fid[prod_ID] == -1) {
	fprintf(stderr, "*** Error:  bad file id sent to put_l3b_close().\n");
	return FAIL;
    }

    /*  check for proper calling order  */
    if (continuity[prod_ID] == INIT) {
	fprintf(stderr, "*** Warning: no data records written to file.\n");
    } else {
        if (continuity[prod_ID] != RECORDS) {
	    fprintf(stderr,
	            "*** Error:  put_l3b_close() called out of order.\n");
	    return FAIL;
        }
    }

    /*  create index and add to grid Vgroup  */
    ndxid[prod_ID] = setupindex(fid[prod_ID]);
    Vaddtagref(gridid[prod_ID], DFTAG_VH, VSQueryref(ndxid[prod_ID]));
    /*  don't forget to count the last [partial] row of bins  */
    extent[prod_ID][row[prod_ID]] = rowcount[prod_ID];
    writeindex(ndxid[prod_ID], start_num, begin[prod_ID], extent[prod_ID], 
		numbin);

    /*  close all open structures  */
    closeindex(ndxid[prod_ID]);
    closedata(prod_ID, mstrid[prod_ID], slvidp[prod_ID], parm_opt);
    total = closegrid(prod_ID, gridid[prod_ID]);


    /*  Remove direcotry path from input file name */           /* LK */
    if ((loc_infiles = strrchr(infiles, '/')) != NULL)          /* LK */
        loc_infiles++;                                          /* LK */
    else                                                        /* LK */
        loc_infiles = infiles;
    SDsetattr(sdfid[prod_ID], "Input Files", DFNT_CHAR, 
		strlen(loc_infiles) + 1, (VOIDP)loc_infiles);

    SDsetattr(sdfid[prod_ID], "L2 Flag Names", DFNT_CHAR, 
		strlen(flag_names) + 1, (VOIDP)flag_names);

    /*
    SDsetattr(sdfid[prod_ID], "L2 Flag Usage", DFNT_CHAR, 
		strlen(flag_use) + 1, (VOIDP)flag_use);

    SDsetattr(sdfid[prod_ID], "L2 Engineering Quality Usage", DFNT_UINT8, 4, 
		(VOIDP)eng_q_use);
		*/

    quot1 = div(smsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);
    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", syear, sday,
            quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid[prod_ID], "time_coverage_start", DFNT_CHAR, strlen(string) + 1,
              (VOIDP)string);
    quot1 = div(emsec, MSECHOUR);
    quot2 = div(quot1.rem, MSECMIN);
    quot3 = div(quot2.rem, MSECSEC);


    /* Fix roll-over problem */
    /* --------------------- */
    if (quot1.quot >= 24) {
      quot1.quot -= 24;
      eday++;

      if (((eyear % 4) == 0) && (eday > 366)) {
	eday = 1;
	eyear++;
      }

      if (((eyear % 4) != 0) && (eday > 365)) {
	eday = 1;
	eyear++;
      }
    }


    sprintf(string, "%4.4d%3.3d%2.2d%2.2d%2.2d%3.3d", eyear, eday,
            quot1.quot, quot2.quot, quot3.quot, quot3.rem);
    SDsetattr(sdfid[prod_ID], "time_coverage_end", DFNT_CHAR, strlen(string) + 1,
              (VOIDP)string);

    SDsetattr(sdfid[prod_ID], "geospatial_lat_units", DFNT_CHAR, 
		strlen(LATUNITS_VAL) + 1, LATUNITS_VAL);
    SDsetattr(sdfid[prod_ID], "geospatial_lon_units", DFNT_CHAR,
              strlen(LONUNITS_VAL) + 1, LONUNITS_VAL);
								/* LK */
    if ((get_coords(begin[prod_ID], lastdatabin[prod_ID], lastrec[prod_ID],
 			&nlat, &slat, &elon, &wlon)) < 0) {
        fprintf(stderr, "*** Error:  get_coords failed for prod_ID %d.\n",
			prod_ID);
	return FAIL;
     }
								/* LK */
    SDsetattr(sdfid[prod_ID], "northernmost_latitude", DFNT_FLOAT32, 1,  /*LK*/
				(VOIDP)&nlat);		        	 /*LK*/
    SDsetattr(sdfid[prod_ID], "southernmost_latitude", DFNT_FLOAT32, 1,  /*LK*/
				(VOIDP)&slat);		      		 /*LK*/
    SDsetattr(sdfid[prod_ID], "easternmost_longitude", DFNT_FLOAT32, 1,  /*LK*/ 
				(VOIDP)&elon);         		         /*LK*/
    SDsetattr(sdfid[prod_ID], "westernmost_longitude", DFNT_FLOAT32, 1,  /*LK*/
				(VOIDP)&wlon);                 	 	 /*LK*/
    SDsetattr(sdfid[prod_ID], "geospatial_lat_max", DFNT_FLOAT32, 1,  /*LK*/
				(VOIDP)&nlat);		        	 /*LK*/
    SDsetattr(sdfid[prod_ID], "geospatial_lat_min", DFNT_FLOAT32, 1,  /*LK*/
				(VOIDP)&slat);		      		 /*LK*/
    SDsetattr(sdfid[prod_ID], "geospatial_lon_max", DFNT_FLOAT32, 1,  /*LK*/ 
				(VOIDP)&elon);         		         /*LK*/
    SDsetattr(sdfid[prod_ID], "geospatial_lon_min", DFNT_FLOAT32, 1,  /*LK*/
				(VOIDP)&wlon);                 	 	 /*LK*/

    SDsetattr(sdfid[prod_ID], "Data Bins", DFNT_INT32, 1, (VOIDP)&total);
    pct = ((float)total / (float)TOTBINS)*100.0;
    SDsetattr(sdfid[prod_ID], "Percent Data Bins", DFNT_FLOAT32, 1,
		(VOIDP)&pct);

    SDsetattr(sdfid[prod_ID], "Units", DFNT_CHAR, 
	      strlen(parm_opt->l3b_units)+1,(VOIDP) parm_opt->l3b_units);


    /*  reset fstcall and lastrec for this product id */
    fstcall[prod_ID] = 1;
    lastrec[prod_ID] = 0;

    /*  close down interfaces and close file  */
    SDend(sdfid[prod_ID]);
    Vend(fid[prod_ID]);
    Hclose(fid[prod_ID]);

    /*  clear ids  */
    fid[prod_ID] = -1;
    sdfid[prod_ID] = -1;
    gridid[prod_ID] = -1;
    ndxid[prod_ID] = -1;
    mstrid[prod_ID] = -1;
    free(slvidp[prod_ID]);					/* LK */

    free(begin[prod_ID]);
    free(extent[prod_ID]);
    free(lastdatabin[prod_ID]);

/*
    The following for loop is replaced by the above line.   This was causing
	the creation of some external files smaller than they should have 
	been.  Changed 'slvid' name to 'slvidp' to indicate it is an array 
	of ptrs.  NCSA helped in resolving this problem on 2/19/95

    for (i=0; i<NPARMS; i++)
	free(slvid[i]);
*/
    continuity[prod_ID] = CLEAR;
    return SUCCEED;
}


/*------------------------------------------------------------------------------
    Function: init

    Returns: <none>

    Description:
        The function init performs some important one-time-only
	initialization of some global data structures.  This function
	is intentionally declared *INACCESSIBLE* to functions outside
	of this file.

    Arguments: <none>
    
    Notes:
	o This function is declared "PRIVATE" (="static") because only
	  functions in this file need to access it.  It is not intended
	  to be called by processing programs.
	o The function appears in this file (as opposed to l3b_misc.c)
	  because it intializes some global variables that are only
	  visible from this file.

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

------------------------------------------------------------------------------*/
PRIVATE void
init(void)
{
    register int i;
    /*    int32 basebin[NUMROWS+1];*/
    extern int32_t *NUMBIN, *BASEBIN;
    /*    extern float *LATBIN;*/

    /*  invalidate all fids  */
    for (i=0; i<MAX_OUT; i++) {
	fid[i] = -1;
   	lastrec[i] = 0;
	fstcall[i] = 1;
    }

    /*  init index arrays  */
    /*  this code is adapted into C from FORTRAN code originally  **

    **  developed by Jim Brown (University of Miami)		  */
    lastbin[0] = 0;

    /*
    basebin[1] = 1;
    latbin[1] = (90.0 / NUMROWS) - 90.0;
    numbin[1] = (int)(cos(latbin[1]*RADCONV) * (2 * NUMROWS) + 0.5);
    for (i=2; i<=NUMROWS; i++) {
        latbin[i] = ((float64)i - 0.5) * (180.0 / NUMROWS) - 90.0;
        numbin[i] = (int)(cos(latbin[i]*RADCONV) * (2 * NUMROWS) + 0.5);
        basebin[i] = basebin[i-1] + numbin[i-1];
    }
    */

    numbin = (int32 *) (NUMBIN-1);

    for (i=1; i<NUMROWS; i++)
	lastbin[i] = BASEBIN[i] - 1;

    lastbin[NUMROWS] = TOTBINS;

    for (i = 1; i < NUMROWS+1; i++)			/* LK */
       start_num[i] = BASEBIN[i-1];	 	/* save the beg bin nos. */
     
	
}


/*------------------------------------------------------------------------------
    Function: get_coords

    Returns: intn (status)
        The return code is FAIL (-1) if an error occurs, SUCCEED (0)
        otherwise.

    Description:
        The function get_coords calls init_bins to Initialize arrays to 
	account for binned data - used in locating an observation and calls 
	binxf_bin2ll to Convert bin # to lat, lon.  It returns the 
	southernmost, northermost latitudes and easternmost and westernmost 
	longitudes to the calling routines through the output arguments. 

    Arguments: (in calling order)
      Type         Name      I/O     Description
      ----         ----      ---     -----------
      int32     *begin        I      starting data bin numbers of each row
      int32	*lastdatabin  I      ending data bin numbers of each row  
      int32     *lastrec      I      ending data bin number of last data row
      float32   *nlat	      O      northernmost latitude
      float32   *slat	      O      southernmost latitude
      float32   *elon	      O      easternmost longitude
      float32   *wlon	      O      westernmost longitude

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Lakshmi Kumar 	 Hughes STX      09/25/95    Original development
        Joel Gales 	 Futuretech      01/05/00    Initialize elon to -180
-----------------------------------------------------------------------------*/
int32
  get_coords(int32 *begin, int32 *lastdatabin, int32 lastrec, float32 *nlat, 
	     float32 *slat, float32 *elon, float32 *wlon)
{

   int32 fst_one = 1, i;
   float32 lat, lon;

   /*
   if ((tot_bins = init_bins(NUMROWS, SEAMLON)) < 0){
      fprintf(stderr, "*** Error:  init_bins returned %d.\n", tot_bins);
      return FAIL;
    }
    tot_bins = TOTBINS;*/

   for (i = 0; i < NUMROWS; i++) {
      if (begin[i] > 0) {
	bin2ll(begin[i], &lat, &lon);
	/*	binxf_bin2ll(begin[i], &lat, &lon, &dellat, &dellon);*/
         if (fst_one) {
            fst_one = 0;
            *slat = lat;
            *wlon = lon;
          }
         *wlon = min(*wlon, lon); 
       }
    }
 
   *elon = -180.0;
   for (i = 0; i < NUMROWS; i++) {
      if (lastdatabin[i] > 0) {
	bin2ll(lastdatabin[i], &lat, &lon);
	/*	binxf_bin2ll(lastdatabin[i], &lat, &lon, &dellat, &dellon);*/
	*elon = max(*elon, lon);
       }
    }
   bin2ll(lastrec, &lat, &lon);
   /*   binxf_bin2ll(lastrec, &lat, &lon, &dellat, &dellon);*/
   *nlat = lat;

   return SUCCEED;
}
