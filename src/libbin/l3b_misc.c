/*------------------------------------------------------------------------------
    File:  l3b_misc.c

    Contents:
	setupgrid   - create and initialize grid Vgroup
	closegrid   - detach from grid Vgroup
	writegeom   - define and write out Geometry-class Vdata
	setupindex  - create and initialize Index-class Vdata
	writeindex  - write out accumulated index data
	closeindex  - detach from Index-class Vdata
	setupmaster - create and initialize DataMaster-class Vdata
	setupslaves - create and initialize DataSlave-class Vdatas
	closedata   - detach from DataMaster- and DataSlave-class Vdatas
	buffbins    - buffer a set of bins
	flushbuff   - flush the bin buffers
	writeattr   - write a generic attribute (not used)

    Notes:
        o These routines have been designed for use by a higher level of
	  SeaWiFS I/O routines.  They are not intended to be called directly
	  by processing applications.
        o Files created by these functions follow the guidelines set forth
          in EOSDIS' "Format System Implementation Guidelines," section 5.6.
          The SEAGrid scheme is used.

    Other relevant files:
        seabin.h    - various #defined constants for level 3 binned data;
                      also #includes hdf.h
        seaprotoi.h - prototypes for this low-layer of level 3 output
	              functions
        put_l3b.c   - a higher layer of level 3 output functions

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
        Lakshmi Kumar    Hughes STX      11/08/95    Changed vdata name of
                                                        eps_68 to eps_78
	Lakshmi Kumar	 Hughes STX	 21/12/95    Fixed a bug in buffbins
	Lakshmi Kumar	 Hughes STX	 01/19/96    Last fix introduced a new
						     bug, and is fixed
	Lakshmi Kumar	 Hughes STX	 05/31/96    Added Hendaccess call to
						     close any open ids before
						     attempting to close a file
        Lakshmi Kumar    Hughes STX      07/30/96    Now passing external file
                                                     ref. no. to HXcreate call
                                                     rather than its id .

        Joel Gales       Futuretech      11/24/99    Modify to handle variable
                                                     products.
------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "seabin.h"
#include "seaprotoi.h"

/*  bin counters  */
intn count[MAX_OUT];
int32 total[MAX_OUT];

/*  static buffers  */
static uchar8 mstrbuff[MAX_OUT][NBINS*MSTRSIZE];
static uchar8 slvbuff[MAX_OUT][NPARMS][NBINS*SLVSIZE];

extern int32 NUMROWS;
extern int32 TOTBINS;

/*------------------------------------------------------------------------------
    Function: setupgrid

    Returns: int32 (ID of grid Vgroup)

    Description:
	The function setupgrid creates the grid Vgroup and gives it a name
	and class.  Bin counters are also initialized.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     prod_ID     I     ID number of product
      int32     fid         I     HDF file ID

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

------------------------------------------------------------------------------*/
int32
#ifdef PROTOTYPE
setupgrid(int32 prod_ID, int32 fid)
#else
setupgrid(prod_ID, fid)
int32 prod_ID, fid;
#endif
{
    int32 gridid;

    /*  create; set name & class  */
    gridid = Vattach(fid, -1, "w");
    Vsetclass(gridid, GRPCLASS);
    Vsetname(gridid, GRPNAME);

    /* initialize bin counters  */
    count[prod_ID] = 0;
    total[prod_ID] = 0;

    return gridid;
}

/*------------------------------------------------------------------------------
    Function: closegrid

    Returns: int32 (number of bins written to grid)

    Description:
	Detaches from the grid Vgroup and returns the number of bins written
	to the grid.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     prod_ID     I     ID number of product
      int32     gridid      I     HDF Vgroup ID of grid

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

------------------------------------------------------------------------------*/
int32
#ifdef PROTOTYPE
closegrid(int32 prod_ID, int32 gridid)
#else
closegrid(prod_ID, gridid)
int32 prod_ID, gridid;
#endif
{
    Vdetach(gridid);
    return total[prod_ID];
}


/*------------------------------------------------------------------------------
    Function: writegeom

    Returns: int32 (ID number of geometry Vdata)

    Description:
	The function writegeom creates and fills out a Geometry-class
	Vdata and returns its Vdata ID.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     HDF file ID

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

------------------------------------------------------------------------------*/
int32
writegeom(int32 fid, int32 nbins)
{
    int32 geomid, geom_ref, i32;
    float64 f64;
    uchar8 dummy[GEOMSIZE];

    /*  create; define fields; set name & class  */
    geomid = VSattach(fid, -1, "w");
    VSfdefine(geomid, "registration", DFNT_INT32, 1);
    VSfdefine(geomid, "straddle", DFNT_INT32, 1);
    VSfdefine(geomid, "bins", DFNT_INT32, 1);
    VSfdefine(geomid, "radius", DFNT_FLOAT64, 1);
    VSfdefine(geomid, "max_north", DFNT_FLOAT64, 1);
    VSfdefine(geomid, "max_south", DFNT_FLOAT64, 1);
    VSfdefine(geomid, "seam_lon", DFNT_FLOAT64, 1);
    VSsetclass(geomid, GEOMCLASS);
    VSsetname(geomid, GEOMNAME);
    VSsetfields(geomid,
	"registration,straddle,bins,radius,max_north,max_south,seam_lon");

    /*  copy data into buffer  */
    i32 = CENTER;
    memcpy(&dummy[0],&i32,sizeof(int32));
    i32 = NO;
    memcpy(&dummy[4],&i32,sizeof(int32));
    i32 = nbins;
    memcpy(&dummy[8],&i32,sizeof(int32));
    f64 = RADIUS;
    memcpy(&dummy[12],&f64,sizeof(float64));
    f64 = MAXNORTH;
    memcpy(&dummy[20],&f64,sizeof(float64));
    f64 = MAXSOUTH;
    memcpy(&dummy[28],&f64,sizeof(float64));
    f64 = SEAMLON;
    memcpy(&dummy[36],&f64,sizeof(float64));

    /*  write out buffer and detach  */
    VSwrite(geomid, dummy, 1, FULL_INTERLACE);
    geom_ref = VSQueryref(geomid);
    VSdetach(geomid);

    return geom_ref;
}


/*------------------------------------------------------------------------------
    Function: setupindex

    Returns: int32 (ID number of Index Vdata)

    Description:
	The function setupindex creates an Index-class Vdata and returns
	its Vdata ID number.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     HDF file ID

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
	Lakshmi Kumar	 HITC		 05/17/95    Added start_num field to
							bin_index 
------------------------------------------------------------------------------*/
int32
#ifdef PROTOTYPE
setupindex(int32 fid)
#else
setupindex(fid)
int32 fid;
#endif
{
    int32 ndxid;

    /*  create; define fields; set name & class  */
    ndxid = VSattach(fid, -1, "w");
    if (ndxid == FAIL) {
        printf("VSattach failed for index.\n***Here's the error stack:\n");
        HEprint(stderr, 0);
        printf("***Error stack complete.\n");
        HEclear();
    }
    VSfdefine(ndxid, "row_num", DFNT_INT32, 1);
    VSfdefine(ndxid, "vsize", DFNT_FLOAT64, 1);
    VSfdefine(ndxid, "hsize", DFNT_FLOAT64, 1);
    VSfdefine(ndxid, "start_num", DFNT_INT32, 1); 		/* LK */
    VSfdefine(ndxid, "begin", DFNT_INT32, 1);
    VSfdefine(ndxid, "extent", DFNT_INT32, 1);
    VSfdefine(ndxid, "max", DFNT_INT32, 1);
    VSsetclass(ndxid, NDXCLASS);
    VSsetname(ndxid, NDXNAME);
    VSsetfields(ndxid, "row_num,vsize,hsize,start_num,begin,extent,max");/*LK*/

    return ndxid;
}


/*------------------------------------------------------------------------------
    Function: writeindex

    Returns: intn (status)
	If all goes well, SUCCEED (0) is returned.  Otherwise FAIL (-1)
	is returned.

    Description:
	The function writeindex writes accumulated index information to
	an Index-class Vdata previously created by a call to setupindex.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     HDF Vdata ID
      int32 *   start_num   I     1-D array - # of first bin in each row
      int32 *   begin       I     1-D array - # of first data bin in each row
      int32 *   extent      I     1-D array - # of bins in each row
      int32 *   maxbin      I     1-D array - max # of bins in each row

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
        Lakshmi Kumar    HITC            05/17/95    Added start_num field to
							bin_index vdata
------------------------------------------------------------------------------*/
intn
#ifdef PROTOTYPE
writeindex(int32 ndxid, int32 *start_num, int32 *begin, int32 *extent, 
		int32 *maxbin)
#else
writeindex(ndxid, start_num, begin, extent, maxbin)
int32 ndxid, *start_num, *begin, *extent, *maxbin;
#endif
{
    uchar8 *buffer, *b;
    float64 vs, hs;
    int32 row;

    /*  initialize pointers  */
    b = buffer = (uchar8 *)HDgetspace(NUMROWS * NDXSIZE);

    /*  copy index info into buffer  */
    for (row=0; row<NUMROWS; row++) {
        memcpy(b,&row,sizeof(int32));
        b+=sizeof(int32);
        vs = (float64)180.0 / (float64)(NUMROWS / 1.0);
        memcpy(b,&vs,sizeof(float64));
        b+=sizeof(float64);
        hs = (float64)360.0 / (float64)(maxbin[row+1]);
        memcpy(b,&hs,sizeof(float64));
        b+=sizeof(float64);
        memcpy(b,&start_num[row+1],sizeof(int32));
        b+=sizeof(int32);
        memcpy(b,&begin[row+1],sizeof(int32));
        b+=sizeof(int32);
        memcpy(b,&extent[row+1],sizeof(int32));
        b+=sizeof(int32);
        memcpy(b,&maxbin[row+1],sizeof(int32));
        b+=sizeof(int32);

    }

    /*  write out and free buffer  */
    VSwrite(ndxid, buffer, NUMROWS, FULL_INTERLACE);
    HDfreespace(buffer);

    return SUCCEED;
}

/*------------------------------------------------------------------------------
    Function: closeindex

    Returns: intn (status)
	If all goes well, SUCCEED (0) is returned.  Otherwise FAIL (-1)
	is returned.

    Description:
	The function closeindex detaches from the Index-class Vdata.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     ndxid       I     HDF Vdata ID

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

------------------------------------------------------------------------------*/
intn
#ifdef PROTOTYPE
closeindex(int32 ndxid)
#else
closeindex(ndxid)
int32 ndxid;
#endif
{
    VSdetach(ndxid);
    return SUCCEED;
}


/*------------------------------------------------------------------------------
    Function: setupmaster

    Returns: int32 (ID number of DataMaster Vdata)

    Description:
	The function setupmaster creates the DataMaster-class Vdata, writes
	out a number of dummy records determined by the output product size,
	and converts the Vdata to an appendable (linked- block) element.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     HDF file ID
      int32     fsize       I     type number for output product

    Notes:
	o Writing out dummy records is necessary at this point in order to
	  force a reasonable size for the initial data block of the data
	  element.  This should have the result of reducing overall linked-
	  block overhead both in file space and in future access time.

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
	Lakshmi Kumar	 Hughes STX	 02/21/95    Added sel_cat field
						     (ref: V4.2 I/O & V2.6
						      product specs)
------------------------------------------------------------------------------*/
int32
#ifdef PROTOTYPE
setupmaster(int32 fid, int32 fsize)
#else
setupmaster(fid, fsize)
int32 fid, fsize;
#endif
{
    int32 mstrid, ref, stat;
    uchar8 *buffer;

    /*  create; define fields; set name & class  */
    mstrid = VSattach(fid, -1, "w");
    if (mstrid == FAIL) {
        printf("VSattach failed for master.\n***Here's the error stack:\n");
        HEprint(stderr, 0);
        printf("***Error stack complete.\n");
        HEclear();
    }
    VSfdefine(mstrid, "bin_num",  DFNT_INT32,   1);
    VSfdefine(mstrid, "nobs",     DFNT_INT16,   1);
    VSfdefine(mstrid, "nscenes",  DFNT_INT16,   1);
    VSfdefine(mstrid, "time_rec", DFNT_INT16,   1);
    VSfdefine(mstrid, "weights",  DFNT_FLOAT32, 1);
    VSfdefine(mstrid, "sel_cat",  DFNT_INT8,    1);	
    VSfdefine(mstrid, "flags_set", DFNT_INT32,  1);
    VSsetclass(mstrid, MSTRCLASS);
    VSsetname(mstrid, MSTRNAME);
    VSsetfields(mstrid, 
	"bin_num,nobs,nscenes,time_rec,weights,sel_cat,flags_set");

    /*  create and write out zeroed dummy buffer based on fsize  */
    switch(fsize) {
    case SCENE:
        buffer = (uchar8 *)calloc(SCNINIT, MSTRSIZE);
        if (buffer == NULL) {
            fprintf(stderr, "ERROR: Out of memory in setupmaster()\n");
            exit(1);
        }
        stat = VSwrite(mstrid, buffer, SCNINIT, FULL_INTERLACE);
        break;
    case DAILY:
        buffer = (uchar8 *)calloc(DAYINIT, MSTRSIZE);
        if (buffer == NULL) {
            fprintf(stderr, "ERROR: Out of memory in setupmaster()\n");
            exit(1);
        }
        stat = VSwrite(mstrid, buffer, DAYINIT, FULL_INTERLACE);
        break;
    case OTHER:
        buffer = (uchar8 *)calloc(BIGINIT, MSTRSIZE);
        if (buffer == NULL) {
            fprintf(stderr, "ERROR: Out of memory in setupmaster()\n");
            exit(1);
        }
        stat = VSwrite(mstrid, buffer, BIGINIT, FULL_INTERLACE);
        break;
    default:
        fprintf(stderr, "Bad fsize value\n");
    }


    if (stat == FAIL) {
        fprintf(stderr, "Failure writing 'DataMaster' Vdata.\n");
        HEprint(stderr, 0);
        exit(1);
    }

    /*  force update to file  */
    ref = VSQueryref(mstrid);
    VSdetach(mstrid);
    /*  re-attach and make appendable with reasonble block size  */
    mstrid = VSattach(fid, ref, "w");
    VSappendable(mstrid, (BLKSIZE * MSTRSIZE));
    /*  seek to beginning to overwrite dummy records  */
    VSseek(mstrid, 0);

    free(buffer);

    return mstrid;
}


/*------------------------------------------------------------------------------
    Function: setupslaves

    Returns: int32 * (array of ID numbers of DataSlave Vdatas)

    Description:
	The function setupslaves creates the DataSlave-class Vdatas, writes
	out a number of dummy data records determined by the output product
	file size, and converts the Vdatas to appendable (linked-block)
	elements if product is a SCENE file, or external elements otherwise.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     HDF file ID
      int32     fsize       I     type number for output product
      char *    l3b_path    I     name of master file
      int32 *   parm_opt    I     option array of parameters to be written 

    Notes:
	o Only one dummy record is needed for each DataSlave-class Vdata
	  in a non-SCENE product because solitary external elements are
	  inherently appendable without using linked-blocks.
	o The names of the files that hold any external Vdatas are
	  synthesized from the name of the master file by adding the
	  suffix ".x##", where the ## represents a two-digit decimal
	  number on the range [00, 01, 02,..., NPARM].

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
        Doug Ilg         Hughes STX      06/17/94    Added code to put main
						     file name at biginning
						     of external element files
    	Lakshmi Kumar	 Hughes STX	 11/09/94    Fixed memory leak 
						     Added code to remove
					  	     direcotry path from file
						     name	
	Lakshmi Kumar    Hughes STX	 11/08/95    Changed vdata name of 
							eps_68 to eps_78
	Lakshmi Kumar	 Hughes STX	 05/31/96    Added Hendaccess call to
						     close any open ids before
						     attempting to close a file

        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure.
                                                     Get fieldnames from this
                                                     structure.
                                                     Replace parm_opt[i] with
                                                     parm_opt->code[i].

        Joel Gales       Futuretech      03/14/00    Increase fname & prodname
                                                     to 256 characters.
------------------------------------------------------------------------------*/
int32 *
#ifdef PROTOTYPE
setupslaves(int32 fid, int32 fsize, char *l3b_path, l3b_prod *parm_opt)
#else
setupslaves(fid, fsize, l3b_path, parm_opt)
int32 fid, fsize;
char *l3b_path;
l3b_prod *parm_opt;
#endif
{
    int32 *slvid, aid;
    int32 stat, i, slvref;
    uchar8 *buffer;
    char fname[256];
    char mname[512];
    char prodname[256];
    char *loc_l3b_path;

    /*  remove direcotry path from file name */ 		    /* LK */
    if ((loc_l3b_path = strrchr(l3b_path, '/')) != NULL)	    /* LK */
        loc_l3b_path++;						    /* LK */
    else							    /* LK */
        loc_l3b_path = l3b_path;				    /* LK */

    /*  make room for slave IDs and open slave Vdatas  */
    slvid = (int32 *)malloc(NPARMS*sizeof(int32));
    for (i=0; i<NPARMS; i++) {
      if (parm_opt->code[i] == 1)
      {
        slvid[i] = VSattach(fid, -1, "w");
        if (slvid[i] == FAIL) {
            printf("VSattach failed for slaves.\n***Error stack:\n");
            HEprint(stderr, 0);
            printf("***Error stack complete.\n");
            HEclear();
        }
      }
    }


    /* Define output fields from product names stored in parm_opt (JMG) */
    /* ---------------------------------------------------------------- */
    for (i=0; i<NPARMS; i++) {
      if (parm_opt->code[i] == 1) {
	strcpy(prodname, parm_opt->prodname[i]);
	strcat(prodname, "_sum");
	VSfdefine(slvid[i], prodname,    DFNT_FLOAT32, 1);

	strcpy(prodname, parm_opt->prodname[i]);
	strcat(prodname, "_sum_sq");
	VSfdefine(slvid[i], prodname,    DFNT_FLOAT32, 1);

	VSsetclass(slvid[i], SLVCLASS);
	VSsetname(slvid[i], parm_opt->prodname[i]);

	strcpy(prodname, parm_opt->prodname[i]);
	strcat(prodname, "_sum,");
	strcat(prodname, parm_opt->prodname[i]);
	strcat(prodname, "_sum_sq");
	VSsetfields(slvid[i], prodname);
      }
    }


    /*  define fields; set name & class for each slave  (old code) */
    /*
    if (parm_opt->code[0] == 1) {
    VSfdefine(slvid[0], "nLw_412_sum",    DFNT_FLOAT32, 1);
    VSfdefine(slvid[0], "nLw_412_sum_sq", DFNT_FLOAT32, 1);
    VSsetclass(slvid[0], SLVCLASS);
    VSsetname(slvid[0], SNAME0);
    VSsetfields(slvid[0], "nLw_412_sum,nLw_412_sum_sq");
    }
    .
    .
    .
    */


    /*  make dummy buffer depending on fsize  */
    if (fsize == SCENE) {
	buffer = (uchar8 *)calloc(SCNINIT, SLVSIZE);
    } else {
	buffer = (uchar8 *)calloc(1, SLVSIZE);
    }

    /*  write out dummy record(s) for each slave  */
    for (i=0; i<NPARMS; i++) {

      if (parm_opt->code[i] == 1) {

        if (fsize == SCENE) {
            stat = VSwrite(slvid[i], buffer, SCNINIT, FULL_INTERLACE);
        } else {
            stat = VSwrite(slvid[i], buffer, 1, FULL_INTERLACE);
	}

        if (stat == FAIL) {
            fprintf(stderr, "Failure writing 'DataSlave' Vdata %d.\n", i);
            HEprint(stderr, 0);
            exit(1);
        }

        slvref = VSQueryref(slvid[i]);
	/*  force write to file  */
        VSdetach(slvid[i]);
	strncpy(mname, loc_l3b_path, 512);			    /* LK */
/*
	strncpy(mname, l3b_path, 512);   
*/

	/*  make appendable or external, as appropriate  */
        slvid[i] = VSattach(fid, slvref, "w");                      /* LK */
	if (fsize == SCENE) {
            VSappendable(slvid[i], (BLKSIZE * SLVSIZE));
	} else {
	    static FILE *sfile;

	    /*  synthesize file name  */
            strcpy(fname, l3b_path);
	    sprintf(&fname[strlen(fname)], ".x%02d", i);
	    /*            sprintf(&fname[strlen(fname)], ".%s", parm_opt->prodname[i]);*/

	    if (strlen(fname) >= 256) {
	      printf("\"fname\" greater than 256 characters\n");
	      exit(1);
	    }

	    /*  open file and write name of main file  **
	    **  in first 512 bytes                     */
	    sfile = fopen(fname, "w");
	    stat = fwrite(mname, 512, 1, sfile);
	    fclose(sfile);

	    /*  convert to external element  */
	    aid = HXcreate(fid, DFTAG_VS, (uint16)slvref, fname, 512, 0);
	    Hendaccess(aid);
	}
        VSdetach(slvid[i]);

	/*  re-attach and seek to beginning to overwrite dummy record(s)  */
        slvid[i] = VSattach(fid, slvref, "w");                      /* LK */
/*
        slvid[i] = VSattach(fid, slvid[i], "w");
*/
        VSseek(slvid[i], 0);
      }
    }
    free(buffer);						    /* LK */
    return slvid;
}


/*------------------------------------------------------------------------------
    Function: closedata

    Returns: intn (status)
	If all goes well, SUCCEED (0) is returned.  Otherwise FAIL (-1)
	is returned.

    Description:
	The function closedata calls flushbuff to flush all data to the file
	and detaches from DataMaster- and DataSlave-class Vdatas.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     prod_ID     I     ID number of product
      int32     mstrid      I     ID number of master Vdata
      int32 *   slvid       I     array of IDs of slave Vdatas
      int32 *   parm_opt    I     option array of parameters to be written 

    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure.
                                                     Replace parm_opt[i] with
                                                     parm_opt->code[i].
------------------------------------------------------------------------------*/
intn
#ifdef PROTOTYPE
closedata(int32 prod_ID, int32 mstrid, int32 *slvid, l3b_prod *parm_opt)
#else
closedata(prod_ID, mstrid, slvid, parm_opt)
int32 prod_ID, mstrid, *slvid;
l3b_prod *parm_opt;
#endif
{
    intn i;

    /*  flush bins left in buffer  */
    flushbuff(prod_ID, mstrid, slvid, parm_opt);

    /*  detach from all data Vdatas  */
    VSdetach(mstrid);

    for (i=0; i<NPARMS; i++) {
      if (parm_opt->code[i] == 1)
        VSdetach(slvid[i]);
    }

    return SUCCEED;
}


/*------------------------------------------------------------------------------
    Function: buffbins

    Returns: intn (status)
	If all goes well, SUCCEED (0) is returned.  Otherwise FAIL (-1)
	is returned.

    Description:
	The function buffbins maintains one data buffer for each DataMaster-
	and DataSlave-class Vdata.  When the buffers fill up, they are
	written out by a call to flushbuff.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     prod_ID     I     ID number of product
      int32     mstrid      I     ID number of master Vdata
      int32 *   slvid       I     array of ID numbers of slave Vdatas
      int32     nrec        I     number of bins being buffered
      int32 *   binno       I     bin number for each bin
      int16 *   nobs        I     # of observations contributing to each bin
      int32 *   ttag        I     time record for each bin
      int16 *   nscenes     I     # of scenes contributing for each bin
      float32 * weights     I     statistical weight for each bin
      int8    * sel_cat	    I     selection category for each binno list
      int16   * flags_set   I     flags corresponding to parent l2_flags
      float32 * l3b_data    I     3-D array of sums and sum-squares
      int32 *   parm_opt    I     option array of parameters to be written 


    Notes:
	o All pointer arguments except slvid, and l3b_data are pointers to
	  one-dimensional arrays of values (one per bin).
        o l3b_data must be stored in such a way as to make it equivalent to
          the [row-major] C-language declaration
                          float32 l3b_data[12][nrec][2];
          Where 12 is the number of parameters in the data set, nrec is the
          number of records being written by this call, and 2 is the number
          of statistics stored for each parameter (i.e., sum and sum-squared).
	o The size of the buffers (in bins) is equal to the number of bins
	  in the equatorial row of the grid (NBINS from seaproto.h).

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development
	Doug Ilg	 Hughes STX	 12/09/93    Changed dimensions of
						     input buffer to:
						     nrec x NPARMS x 2
						     Out-of-range bins
						     get ignored on writes
	Lakshmi Kumar	 Hughes STX	 02/13/95    updated to include sel_cat
						     field in vdata BinList
						     (ref. I/O specs v4.2)
 	Lakshmi Kumar	 HITC		 05/18/95    Added flags_set field to 
						     BinList vdata (ref. V4.3
						     I/O Specs.)
	Lakshmi Kumar	 HSTX		 21/12/95    Fixed a bug that was 
						     causing buffering problem
						     It reinitializes buf ptrs
						     if it is fst call for that
						     product ID
	Lakshmi Kumar	 HSTX		 01/19/96    Prev fix caused a bug and
						     if fixed by setting 
						     fstcall[prod_ID] = 0 in 
						     1st if condition. 

        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure.
                                                     Replace parm_opt[i] with
                                                     parm_opt->code[i].
------------------------------------------------------------------------------*/
intn
#ifdef PROTOTYPE
buffbins(int32 *fstcall, int32 prod_ID, int32 mstrid, int32 *slvid, 
	 int32 nrec, int32 *binno, int16 *nobs, int16 *ttag, int16 *nscenes, 
	 float32 *weights, int8 *sel_cat, int32 *flags_set, float32 *l3b_data,
         l3b_prod *parm_opt)
#else
buffbins(fstcall, prod_ID, mstrid, slvid, nrec, binno, nobs, ttag, 
	 nscenes, weights, sel_cat, flags_set, l3b_data, parm_opt)
int32 *fstcall, prod_ID, mstrid, *slvid, nrec, *binno, *flags_set;
int16 *nobs, *nscenes, *ttag;
int8  *sel_cat;
float32 *weights, *l3b_data;
l3b_prod *parm_opt;
#endif
{

    static intn first = 1;

    static uchar8 *bn[MAX_OUT], *ns[MAX_OUT], *nsc[MAX_OUT],
    		  *tt[MAX_OUT], *sw[MAX_OUT], *sc[MAX_OUT], 		/*LK*/
		  *fs[MAX_OUT];						/*LK*/
    static float32 *x0[MAX_OUT][NPARMS], *d0[MAX_OUT][NPARMS];

    static int32 slvincr = NPARMS * 2;
    static int32 mstrincr = MSTRSIZE;

    register intn i, j;

    if (first) {
	/*  initialize all buffer and input pointers  */
	for (i=0; i<MAX_OUT; i++) {
            bn[i]  = (uchar8 *) &mstrbuff[i][0];
            ns[i]  = bn[i]  + sizeof(*binno);
            nsc[i] = ns[i]  + sizeof(*nobs);
            tt[i]  = nsc[i] + sizeof(*nscenes);
            sw[i]  = tt[i]  + sizeof(*ttag);
	    sc[i]  = sw[i]  + sizeof(*weights); 		/* LK */
	    fs[i]  = sc[i]  + sizeof(*sel_cat);			/* LK */
	    for (j=0; j<NPARMS; j++) {
              if (parm_opt->code[j] == 1)
              {
	        x0[i][j] = (float32 *) &slvbuff[i][j][0];
                d0[i][j] = (float32 *) &l3b_data[j*2];
              }
	    }
	}
        first = 0;
        fstcall[prod_ID] = 0;
     }  
    else {
	/*  initialize just the input pointers  */
	for (j=0; j<NPARMS; j++) {
          if (parm_opt->code[j] == 1)
            d0[prod_ID][j] = (float32 *) &l3b_data[j*2];
         }
     }
    if (fstcall[prod_ID]) {				   	/*LK*/
	/* initialize buffer and input pointers for this product ID */  /*LK*/
           bn[prod_ID]  = (uchar8 *) mstrbuff[prod_ID];			/*LK*/
           ns[prod_ID]  = bn[prod_ID]  + sizeof(*binno);		/*LK*/
           nsc[prod_ID] = ns[prod_ID]  + sizeof(*nobs);			/*LK*/
           tt[prod_ID]  = nsc[prod_ID] + sizeof(*nscenes);		/*LK*/	
           sw[prod_ID]  = tt[prod_ID]  + sizeof(*ttag);			/*LK*/
           sc[prod_ID]  = sw[prod_ID]  + sizeof(*weights);         	/*LK*/
           fs[prod_ID]  = sc[prod_ID]  + sizeof(*sel_cat);         	/*LK*/
           for (j=0; j<NPARMS; j++) {				        /*LK*/
             if (parm_opt->code[j] == 1)                                      /*GF*/
             {
               x0[prod_ID][j] = (float32 *) &slvbuff[prod_ID][j][0];	/*LK*/	
	       d0[prod_ID][j] = (float32 *) &l3b_data[j*2];		/*LK*/
             }                                                          /*GF*/
            }								/*LK*/	
           fstcall[prod_ID] = 0;					/*LK*/
      }									/*LK*/

    /*  main buffering loop  */
    for (i=0; i<nrec; i++) {
	if ((binno[i] > 0) && (binno[i] < TOTBINS)) {
	    /*  copy DataMaster data into buffer  */
            memcpy(bn[prod_ID],  (char *)&binno[i],   sizeof(*binno));
            memcpy(ns[prod_ID],  (char *)&nobs[i],    sizeof(*nobs));
            memcpy(nsc[prod_ID], (char *)&nscenes[i], sizeof(*nscenes));
            memcpy(tt[prod_ID],  (char *)&ttag[i],    sizeof(*ttag));
            memcpy(sw[prod_ID],  (char *)&weights[i], sizeof(*weights)); 
	    memcpy(sc[prod_ID],  (char *)&sel_cat[i], sizeof(*sel_cat)); /*LK*/
	    memcpy(fs[prod_ID],  (char *)&flags_set[i], sizeof(*flags_set));

	    /*  copy DataSlave data into buffers  */
	    for (j=0; j<NPARMS; j++) {
              if (parm_opt->code[j] == 1)
                memcpy(x0[prod_ID][j],  (char *)d0[prod_ID][j],
	    	       2 * sizeof(*d0[prod_ID][j]));
	    }

	    /*  increment bin counters  */
            count[prod_ID]++;
	    total[prod_ID]++;

	    /*  buffers full?  */
            if (count[prod_ID] == NBINS) {
                /*  flush buffers and reset pointers/count  */
                bn[prod_ID]  = (uchar8 *) mstrbuff[prod_ID];
                ns[prod_ID]  = bn[prod_ID]  + sizeof(*binno);
                nsc[prod_ID] = ns[prod_ID]  + sizeof(*nobs);
                tt[prod_ID]  = nsc[prod_ID] + sizeof(*nscenes);
                sw[prod_ID]  = tt[prod_ID]  + sizeof(*ttag);
		sc[prod_ID]  = sw[prod_ID]  + sizeof(*weights);		/*LK*/
 		fs[prod_ID]  = sc[prod_ID]  + sizeof(*sel_cat);		/*LK*/
	        for (j=0; j<NPARMS; j++) {
                  if (parm_opt->code[j] == 1)
	            x0[prod_ID][j] = (float32 *) &slvbuff[prod_ID][j][0];
	        }
                flushbuff(prod_ID, mstrid, slvid, parm_opt);
                count[prod_ID] = 0;
            } else {
                /*  increment output buffer pointers  */
                bn[prod_ID]  += mstrincr;
                ns[prod_ID]  += mstrincr;
                nsc[prod_ID] += mstrincr;
                tt[prod_ID]  += mstrincr;
                sw[prod_ID]  += mstrincr;
		sc[prod_ID]  += mstrincr;			/* LK  */	
		fs[prod_ID]  += mstrincr;			/* LK  */
	        for (j=0; j<NPARMS; j++) {
                  if (parm_opt->code[j] == 1)
                  {
	            x0[prod_ID][j]++; x0[prod_ID][j]++;
                  }
	        }
            }
	}
        /*  increment input buffer pointers  */
	for (j=0; j<NPARMS; j++) {
          if (parm_opt->code[j] == 1)
	    d0[prod_ID][j] += slvincr;
	}
    }
  return SUCCEED;						/*  LK  */
}


/*------------------------------------------------------------------------------
    Function: flushbuff

    Returns: intn (status)
	If all goes well, SUCCEED (0) is returned.  Otherwise FAIL (-1)
	is returned.

    Description:
	The function flushbuff writes the contents of all data buffers for
	the given prod_ID to the appropriate Vdatas in the appropriate file.

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     prod_ID     I     ID number of product
      int32     mstrid      I     ID number of master Vdata
      int32 *   slvid       I     array of ID numbers of slave Vdatas
      int32 *   parm_opt    I     option array of parameters to be written 


    Notes:

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

        Joel Gales       Futuretech      11/24/99    Change parm_opt from
                                                     int32 to structure.
                                                     Replace parm_opt[i] with
                                                     parm_opt->code[i].
------------------------------------------------------------------------------*/
intn
#ifdef PROTOTYPE
flushbuff(int32 prod_ID, int32 mstrid, int32 *slvid, l3b_prod *parm_opt)
#else
flushbuff(prod_ID, mstrid, slvid, parm_opt)
int32 prod_ID, mstrid, *slvid;
l3b_prod *parm_opt;
#endif
{
    register intn i;
    int stat;

    if (count[prod_ID] > 0) {
        VSwrite(mstrid, mstrbuff[prod_ID], count[prod_ID], FULL_INTERLACE);

        for (i=0; i<NPARMS; i++) {
          if (parm_opt->code[i] == 1) {
            stat = VSwrite(slvid[i], &slvbuff[prod_ID][i][0], count[prod_ID],
	                   FULL_INTERLACE);

	    if (stat == FAIL) {
	      fprintf(stderr, "Failure writing 'DataSlave' Vdata: %d.\n",i);
	      HEprint(stderr, 0);
	      exit(1);
	    }

	  }
        }
    }

    return SUCCEED;
}

/*------------------------------------------------------------------------------
    Function: writeattr

    Returns: int32 (ID number of attribute Vdata)

    Description:
	The function writeattr writes a local attribute to a grid and returns
	its Vdata ID. 

    Arguments: (in calling order)
      Type      Name        I/O   Description
      ----      ----        ---   -----------
      int32     fid         I     HDF file ID
      char *    name        I     attribute name
      int32     nt          I     HDF number type constant
      intn	count       I     number of objects of type nt to be written
      void *    data        I     array of data to be written

    Notes:
	o This function provides general attribute capability similar to
	  that provided by the HDF SDS interface.
	o The attributes written by this function are identical in structure
	  to those written by the HDF SDS interface.
	o This function is not used by any functions in put_l3b.c.  It is
	  included here for completeness.

    Modification history:
          Programmer     Organization      Date      Description of change
        --------------   ------------    --------    ---------------------
        Doug Ilg         Hughes STX      10/14/93    Original development

------------------------------------------------------------------------------*/
int32
#ifdef PROTOTYPE
writeattr(int32 fid, char *name, int32 nt, intn count, VOIDP data)
#else
writeattr(fid, name, nt, count, data)
int32 fid;
char  *name;
int32 nt;
intn  count;
VOIDP data;
#endif
{
    int32 order, size, status;
/* don't forget arg checks */

    /*  chars take a size of 1 and an order    **
    **  equal to the length of the string.     **
    **  Other types take an order of 1 and     **
    **  a size equal to the number of values.  **
    **  See the SDS interface for details.     */
    if (nt == DFNT_CHAR) {
        order = count;
        size = 1;
    } else {
        order = 1;
        size = count;
    }

    status = VHstoredatam(fid, "Values", (uchar8 *) data, size, nt, name,
   			  ATTR, order);

    return status;
}
