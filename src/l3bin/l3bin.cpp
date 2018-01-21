#include <stdio.h>
#include <math.h>
#include <time.h>
#include <libgen.h>
#include <sys/types.h>

#include "netcdf.h"

#include "seabin.h"
#include "getl3b.h"
#include "l3bin_input.h"
#include <timeutils.h>
#include <genutils.h>
#include "hdf_bin.h"
#include "sensorInfo.h"


#define MAXNFILES 128

#define BYTE    unsigned char

#define BINCHECK -1

#define VERSION "4.20"

static instr input;

extern "C" int l3bin_input(int argc, char **argv, instr *input);

/*
  Revision 4.20 09/03/15
  Add support for time_rec
  J. Gales

  Revision 4.10 08/03/15
  Add support for L2C SMAP binfile input
  J. Gales

  Revision 4.03 04/28/15
  In the case of reduce_fac != 1 for quality field allow
  "i" and "n_write" to differ by 2 or less due to adjacent
  map rows of integralized sinesodial projection being of different lengths. 
  J. Gales

  Revision 4.02 08/01/14
  Check for both "Mission" and "mission" for HDF5 (Aquarius) test
  J. Gales

  Revision 4.01 03/19/14
  Add support for netCDF4 compression
  Modify isHDF5/isCDF4 code to handle new netcdf4 library
  Check for Aquarius mission attribute if possible
  J. Gales

  Revision 4.00 01/30/14
  Add support for NETCDF4 I/O
  J. Gales

  Revision 3.56 11/06/13
  Write Start/End Time metadata to HDF5 (Aquarius) files
  Fixes bad "Start/End Day" metadata for monthly binfiles.
  J. Gales

  Revision 3.55 10/23/13
  Read input Start/End Time metadata for HDF5 (Aquarius) files
  Fixes bad "End Day" metadata for 7-day/monthly binfiles.
  J. Gales

  Revision 3.54 08/05/13
  Initial n_data_prod to 0 in hdf5_bin constructor (in bin_io.cpp)
  J. Gales

  Revision 3.53 06/27/13
  Write "Sensor Name", "Units", "Software Version" and "Input Parameters"
  metadata for NETCDF4 files
  J. Gales

  Revision 3.52 06/27/13
  Write "instrumentInformation", "Product Type", "L2 Flag Names",
  "Period ...", "Start/End Time" and "Start/End Orbit" metadata 
  for NETCDF4 files
  J. Gales

  Revision 3.51 06/22/13
  Add "Binning Scheme" metadata
  J. Gales

  Revision 3.50 06/14/13
  Add support for NETCDF4
  J. Gales

  Revision 3.42 07/19/12
  Fix bad qual in GBR SST bug
  Clean up qual read code
  J. Gales

  Revision 3.41 04/10/12
  Fix in_qual_buf 255 fill bug
  J. Gales

  Revision 3.40 03/04/12
  Add support for products with quality fields for reduce_fac > 1
  J. Gales

  Revision 3.36 08/05/11
  Don't perform quality write check for granules without quality field
  J. Gales

  Revision 3.35 07/29/11
  Fix problem with incorrect quality for longitudinal subsampling
  J. Gales

  Revision 3.30 07/12/11
  Use inheritance and virtual functions
  J. Gales

  Revision 3.28 06/27/11
  Remove Data_Center and Station metadata
  Add Map Resolution attribute
  J. Gales

  Revision 3.27 09/03/10
  Fix readQual bug in bin_io.cpp
  J. Gales

  Revision 3.26 07/26/10
  Fix writeQual bug in bin_io.cpp
  J. Gales

  Revision 3.25 02/16/10
  Fix problems with multiple data fields for HDF5 bin files
  J. Gales

  Revision 3.20 08/11/09
  Add pversion parameter
  Remove repimg
  J. Gales

  Revision 2.20 11/16/08
  Add support for unit weights and median binning
  J. Gales

  Revision 2.18 09/10/08
  Add checks for lon/lat limit values
  J. Gales

  Revision 2.18 10/11/07
  Add support VERBOSE parameter
  J. Gales

  Revision 2.17 02/16/07
  Add support for lon/lat subsetting
  J. Gales

  Revision 2.16 11/17/06
  Put open, close in l3bin_util.c
  J. Gales

  Revision 2.15 09/14/06
  Fix "bins" count in SEAGrid vdata
  J. Gales

  Revision 2.14 08/28/06
  Change vdata_id[0][2+nprod[ifile]] to vdata_id[ifile][2+nprod[ifile]]
  in "Read quality" section 
  J. Gales
*/

void usage (const char *progname)
{
  printf("%s %s (%s %s)\n",progname,VERSION,__DATE__,__TIME__);

  printf("\nUsage: %s in=input-file out=output-file out_parm=prodlist\n",
	 progname);
  printf("            [reduce_fac=reduce_fac] [noext=noext]\n");
  printf("    input-file  = listfile of input binfiles\n");
  printf("    output-file = output bin filename\n");
  printf("    out_parm    = data products list\n");
  printf("    parfile     = parameter filename\n");
  printf("    reduce_fac  = scale reduction factor (power of 2)\n");
  printf("    loneast     = Easternmost longitude (default=+180)\n");
  printf("    lonwest     = Westernmost longitude (default=-180)\n");
  printf("    latnorth    = Northernmost latitude (default=+90)\n");
  printf("    latsouth    = Southernmost latitude (default=-90)\n");
  printf("    noext       = set to 1 to suppress generation of external files\n");
  printf("                  [default=0]\n");
  printf("    verbose     = Allow more verbose screen messages [default=0]\n");
  printf("    pversion    = Processing Version [default=Unspecified]\n");
  printf("    oformat     = output format: 1 (HDF4), 2 (netCDF4), 3 (HDF5)\n");
  printf("                  default is the same format as input file(s)\n");
  printf("    deflate     = apply internal compression for netCDF output\n");

  exit(0);
}



int main (int argc, char **argv)
{
  intn  i;

  int32 irow;
  int32 kbin;
  int32 iprod;

  int32 ifile, nfiles;

  int32 nread;
  int32 offset;
  int32 offset_out;
  int32 offmin;
  int32 offmax;
  int32 bin_num;
  int32 bin_num_out;
  int32 row_write;
  int32 n_write_total=0;
  int32 n_write=0;
  int32 nprod[MAXNFILES];
  int32 ncols;
  int32 ncols_out;

  int32 nrows;

  int32 reduce_fac;

  float32 wgt;
  float32 *sort_array[MAXNVDATA];

  char buf[FILENAME_MAX];

  float32 *in_sum_buf[MAXNVDATA-2];
  float32 *out_sum_buf[MAXNVDATA-2];
  uint8 *in_qual_buf[MAXNFILES+1], *uint8_buf;

  int status;
  
  char     ptime[17];
  char     proc_con[2048];

  float32 f32;

  float32 minlon;
  float32 maxlon;
  float32 minlat;
  float32 maxlat;

  float32 lat;

  time_t tnow;
  struct tm *tmnow;

  FILE *fp;

  void insertion_sort(float a[], int length);

  setlinebuf(stdout);

  if (argc == 1) usage("l3bin");

  printf("%s %s (%s %s)\n","L3BIN",VERSION,__DATE__,__TIME__);

  if (l3bin_input(argc, argv, &input) != 0) {
    /*  usage(argv[0]);*/
    exit(1);
  }

  get_time(ptime);


  strcpy(proc_con, argv[0]);
  for (i=1; i < argc; i++) {
    strcat(proc_con, " ");
    strcat(proc_con, argv[i]);
  }

  reduce_fac = input.reduce_fac;
  if (reduce_fac != 1 &&
      reduce_fac != 2 &&
      reduce_fac != 4 &&
      reduce_fac != 8 &&
      reduce_fac != 16) {
    printf("Reduction factor must be power of 2 less than 16\n");
    exit(-1);
  }

  if (input.loneast <= input.lonwest) {
    printf("loneast: %f must be greater than lonwest: %f.\n",
	   input.loneast, input.lonwest);
    exit(-1);
  }

  if (input.latnorth <= input.latsouth) {
    printf("latnorth: %f must be greater than latsouth: %f.\n",
	   input.latnorth, input.latsouth);
    exit(-1);
  }

  /* Get lon/lat limits */
  minlon = input.lonwest;
  maxlon = input.loneast;
  minlat = input.latsouth;
  maxlat = input.latnorth;


  /* Determine number of input files */
  /* ------------------------------- */
   nfiles = 0;

   bool isHDF4=false;
   bool isHDF5=false;
   bool isCDF4=false;

   Hdf::hdf_bin *input_binfile[MAXNFILES];

   /* Single HDF input */
   /* ---------------- */
   if (Hishdf(input.infile) == TRUE || H5Fis_hdf5(input.infile) == TRUE) {
    printf("Single HDF input\n");
    nfiles = 1;

    if (Hishdf(input.infile) == TRUE) {
      isHDF4 = true;
      input_binfile[0] = new Hdf::hdf4_bin;
    }

    if (H5Fis_hdf5(input.infile) == TRUE) {
      int ncid;
      char nam_buf[256];
      status = nc_open( input.infile, NC_NOWRITE, &ncid);
      if ( status != NC_NOERR) {
	isHDF5 = true;
	input_binfile[0] = new Hdf::hdf5_bin;
      } else {
	status = nc_get_att( ncid, NC_GLOBAL, "Mission", nam_buf);
        if ( status != NC_NOERR)
          status = nc_get_att( ncid, NC_GLOBAL, "mission", nam_buf);

        if ( (strcmp( nam_buf, "SAC-D Aquarius") == 0) ||
             (strcmp( nam_buf, "SMAP") == 0)) {
              nc_close( ncid);
              isHDF5 = true;
              input_binfile[0] = new Hdf::hdf5_bin;	
	} else {
	  isCDF4 = true;
	  input_binfile[0] = new Hdf::cdf4_bin;
	}
      }
    }

    input_binfile[0]->open( input.infile);
    nprod[0] = input_binfile[0]->nprod();
   } else {
       fp = fopen(input.infile, "r");
       if (fp == NULL) {
           printf("Input listing file: \"%s\" not found.\n", input.infile);
           return -1;
       }
       while(fgets(buf, 256, fp) != NULL) nfiles++;
       fclose(fp);
       printf("%d input files\n", nfiles);


       /* Open L3 input files */
       /* ------------------- */
       fp = fopen(input.infile, "r");
       for (ifile=0; ifile<nfiles; ifile++) {
           fgets(buf, 256, fp);
           buf[strlen(buf)-1] = 0;

           if ( ifile == 0) {
               if (Hishdf(buf) == TRUE) isHDF4 = true;

               if (H5Fis_hdf5(buf) == TRUE) {
                   int ncid;
                   char nam_buf[256];
                   status = nc_open( buf, NC_NOWRITE, &ncid);
                   if ( status != NC_NOERR) {
                       isHDF5 = true;
                   } else {
                       status = nc_get_att( ncid, NC_GLOBAL, 
                                            "Mission", nam_buf);
                       if ( status != NC_NOERR)
                           status = nc_get_att( ncid, NC_GLOBAL, 
                                                "mission", nam_buf);

                       if ( (strcmp( nam_buf, "SAC-D Aquarius") == 0) ||
                            (strcmp( nam_buf, "SMAP") == 0)) {
                           nc_close( ncid);
                           isHDF5 = true;
                       } else {
                           isCDF4 = true;
                       }
                       nc_close( ncid);
                   }
               }
           }
           if ( isHDF4) input_binfile[ifile] = new Hdf::hdf4_bin;
           if ( isHDF5) input_binfile[ifile] = new Hdf::hdf5_bin;
           if ( isCDF4) input_binfile[ifile] = new Hdf::cdf4_bin;

           printf("%d %s\n", ifile, buf);
           input_binfile[ifile]->open( buf);
           nprod[ifile] = input_binfile[ifile]->nprod();

           //printf("open status: %d\n", status);

       } /* ifile loop */

       fclose(fp);
   }


  nrows = input_binfile[0]->nrows;

  /* Generate output product list from 1st input file if DEFAULT */
  /* ----------------------------------------------------------- */
  if (strcmp(input.out_parm, "DEFAULT") == 0) {

    // Read input file product list
    input_binfile[0]->query( input.out_parm);

    //    printf("out_parm: %s\n", &input.out_parm[1]);
  } else {
    strcpy(buf, &input.out_parm[1]);
    buf[strlen(buf)-1] = 0;
    for (i=0; i < (intn)strlen(buf); i++)
      if (buf[i] == ':') buf[i] = ',';
    strcpy(input.out_parm, buf);
  }

  /* Determine active products */
  /* ------------------------- */
  for (ifile=0; ifile<nfiles; ifile++) {
    int status = input_binfile[ifile]->read( input.out_parm);
    if ( status == -1) {
      printf("Not all output products found in input file %d (1-based)\n",
	     ifile);
      exit(-1);
    }
  }


  /* Create output file */
  /* ------------------ */
  Hdf::hdf_bin *output_binfile;

  if (getFileFormatName( input.oformat) == NULL){
      if (isHDF4) strcpy(input.oformat, "HDF4");
      if (isHDF5) strcpy(input.oformat, "HDF5");
      if (isCDF4) strcpy(input.oformat, "netCDF4");
  }

  strcpy(buf, input.ofile);

  if ( strcmp( input.oformat, "HDF4") == 0) {
    output_binfile = new Hdf::hdf4_bin;
    output_binfile->hasNoext = false;
    if (input.noext == 0)
        strcat(buf, ".main");
    else
        output_binfile->hasNoext = true;
  }
  if ( strcmp( input.oformat, "HDF5") == 0) output_binfile = new Hdf::hdf5_bin;
  if ( strcmp( input.oformat, "netCDF4") == 0) output_binfile = new Hdf::cdf4_bin;

  output_binfile->create( buf, nrows/reduce_fac);
  if ( isCDF4) output_binfile->deflate = input.deflate;


  /* Allocate I/O buffers */
  /* -------------------- */
  ncols = 2 * nrows;
  ncols_out = 2 * nrows/reduce_fac;

  for (iprod=0; iprod<nprod[0]; iprod++) {
    if ( input_binfile[0]->active_data_prod[iprod] == true) {
      in_sum_buf[iprod] = (float32 *) calloc(ncols, 2*sizeof(float32));
      out_sum_buf[iprod] = (float32 *) calloc(ncols_out, 2*sizeof(float32));
    }
  } /* iprod loop */


  /* Allocate quality buffer */
  /* ----------------------- */
  for (ifile=0; ifile<nfiles; ifile++) {
    in_qual_buf[ifile] = (uint8 *) calloc(ncols, sizeof(uint8));
  }
  in_qual_buf[nfiles] = (uint8 *) calloc(ncols, sizeof(uint8));

  uint8_buf = (uint8 *) calloc(ncols, sizeof(uint8));


  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
  /* For each scan ... (Main Loop) */
  /* ----------------------------- */
  for (irow=0; irow<nrows; irow++) {

    if ((irow % 500) == 0) {
      time(&tnow);
      tmnow = localtime(&tnow);
      printf("irow:%6d of %8d %s", irow, nrows, asctime(tmnow));
    }

    int32 max_out_kbin;
    if ( (irow % reduce_fac) == 0)
      max_out_kbin = -1;

    // Get basebin and numbin for this input row
    int basebin = input_binfile[0]->get_basebin( irow);
    int numbin = input_binfile[0]->get_numbin( irow);

    double ratio = 1.0;
    if ( reduce_fac > 1)
      ratio = 1.0 * input_binfile[0]->get_numbin(irow) / 
	output_binfile->get_numbin(irow/reduce_fac);


    // If median allocate and initialize storage
    if ( input.median) {
      for (iprod=0; iprod<nprod[0]; iprod++) {
	if ( input_binfile[0]->active_data_prod[iprod] == true) {
	  sort_array[iprod] = 
	    (float32 *) calloc(nfiles*numbin, sizeof(float32));
	  for (int32 i=0; i<nfiles*numbin; i++) sort_array[iprod][i] = -999;
	}
      }
    }


    /* Clear output binlist, sum and quality buffers */
    /* --------------------------------------------- */
    if ((irow % reduce_fac) == 0) {
      row_write = 0;

      output_binfile->clear_binlist();
      for (iprod=0; iprod<nprod[0]; iprod++) {
	if (input_binfile[0]->active_data_prod[iprod] == true) {
	  memset(&out_sum_buf[iprod][0], 0, ncols_out*2*sizeof(float32));
	}
      } /* iprod loop */

      memset(in_qual_buf[nfiles], 255, ncols);
      for (ifile=0; ifile<nfiles; ifile++) {
	memset(in_qual_buf[ifile], 255, ncols);
      }
    }

    /* Get bin & qual info */
    /* ------------------- */
    for (ifile=0; ifile<nfiles; ifile++) {
      input_binfile[ifile]->readBinIndex( irow);

      int ext = input_binfile[ifile]->get_ext();

      /* Read BinList */
      /* ------------ */
      // Get current binlist pointer
      int32 list_ptr = input_binfile[ifile]->get_list_ptr();

      if ( ext > 0) {
	nread = input_binfile[ifile]->readBinList( ext);

	if (nread == -1) {
	  printf("Unable to read bin numbers...: %d\n", ext);
	}
      }

      /* Read quality if vdata exists */
      /* Note: in_qual_buf is "uncompressed" along row */
      if ( input_binfile[ifile]->has_qual()) {
	if ( (irow % reduce_fac) == 0) {
	  int32 ext_qual = ext;
	  for (int32 i=0; i<reduce_fac; i++) {
	    switch (i) {
	    case 0:
	      input_binfile[ifile]->readQual( uint8_buf, ext_qual, list_ptr);
	      break;
	    default:
	      input_binfile[ifile]->readBinIndex( irow+i);
	      ext_qual = input_binfile[ifile]->get_ext();
	      if ( ext_qual != 0) {
		nread = input_binfile[ifile]->readBinList( ext_qual);
                input_binfile[ifile]->readQual( uint8_buf, ext_qual);
	      }
	    } // switch

	    // Update in_qual_buf for (irow+i)
	    if ( ext_qual != 0) {
              int basebin = input_binfile[0]->get_basebin( irow+i);

	      double ratio = 1.0 * input_binfile[0]->get_numbin(irow+i) / 
		output_binfile->get_numbin(irow/reduce_fac);

	      for (kbin=0; kbin<ext_qual; kbin++) {
		bin_num = input_binfile[ifile]->get_bin_num( kbin);
		offset = bin_num - basebin;
                if ( offset < 0) {
                  cout << "bin_num - basebin is negative for ifile: " << ifile;
                  cout << " irow: " << irow << " kbin: " << kbin << endl;
                  cout << "bin_num: " << bin_num << "  basebin: " << basebin 
                       << endl;
                  exit(1);
                }

		int32 j = reduce_fac * (int32) ((offset/ratio) + 0.0);
		if ( (j < ncols) && (in_qual_buf[ifile][j] > uint8_buf[kbin])) {
		  in_qual_buf[ifile][j] = uint8_buf[kbin];
		}
	      }
	    } // ext_qual != 0
	  } // i (reduce_fac) loop

	  // Reset to current row if reduce_fac != 1
          if ( reduce_fac != 1) {
	    input_binfile[ifile]->readBinIndex( irow);
	    nread = input_binfile[ifile]->readBinList( ext, list_ptr);
	  }

	} // if ( (irow % reduce_fac) == 0)
      } // if ( input_binfile[ifile]->has_qual())
    } /* ifile loop */


    /* Find best quality */
    /* ----------------- */
    for (kbin=0; kbin<ncols; kbin++) {
      for (ifile=0; ifile<nfiles; ifile++) {
	int32 j = reduce_fac * (int32) ((kbin/ratio) + 0.0);
	if ( j < ncols) {
	  if (in_qual_buf[ifile][j] < in_qual_buf[nfiles][j])
	    in_qual_buf[nfiles][j] = in_qual_buf[ifile][j];
	}
      }
    }
    
    /* For each file ... */
    /* ----------------- */
    for (ifile=0; ifile<nfiles; ifile++) {

      int beg = input_binfile[ifile]->get_beg();
      int ext = input_binfile[ifile]->get_ext();

      /* ========== If row has data ... ========== */
      /* ----------------------------------------- */
      if ( beg != 0) {
	/*	printf("row has data: %d\n", irow);*/

	/* Determine lon kbin limits */
	/* ------------------------- */
	offmin = 
	  (int32) ((minlon + 180) * (numbin / 360.0) + 0.5);
	offmax = 
	  (int32) ((maxlon + 180) * (numbin / 360.0) + 0.5);


	/* Get data values (sum, sum_sq) for each filled bin in row */
	/* -------------------------------------------------------- */
	int nbins_to_read = ext;
	for (iprod=0; iprod<nprod[ifile]; iprod++) {
	  if (input_binfile[0]->active_data_prod[iprod] == true) {

	    input_binfile[ifile]->readSums( &in_sum_buf[iprod][0], 
					    nbins_to_read, iprod);
	  }
	} /* iprod loop */
	if ( isHDF5 || isCDF4) 
	  input_binfile[ifile]->setDataPtr( nbins_to_read);
	//	row_write = 1;


	/* Skip row if not between minlat & maxlat */
	lat = ((irow + 0.5) / nrows) * 180.0 - 90.0;
	if (lat < minlat || lat > maxlat) {
	  // row_write = 0;
	  continue;
	}

	/* Fill output buffers with input bin data */
	/* --------------------------------------- */
	for (kbin=0; kbin<ext; kbin++) {

	  /* Store bin number */
	  /* ---------------- */
	  bin_num = input_binfile[ifile]->get_bin_num( kbin);
	  offset = bin_num - basebin;
          if ( offset < 0) {
            cout << "bin_num - basebin is negative for ifile: " << ifile;
            cout << " irow: " << irow << " kbin: " << kbin << endl;
            cout << "bin_num: " << bin_num << "  basebin: " << basebin << endl;
            exit(1);
          }

	  /* If bin outside lon range then skip */
	  /* ---------------------------------- */
	  if (offset < offmin || offset > offmax)
	    continue;

	  float weights = input_binfile[ifile]->get_weights( kbin);
	  float time_rec = input_binfile[ifile]->get_time_rec( kbin);

	  /* Skip if not good enough */
	  /* ----------------------- */
	  int32 j = reduce_fac * (int32) ((offset/ratio) + 0.0);
	  if ( j < ncols) {
	    if (in_qual_buf[ifile][j] > in_qual_buf[nfiles][j]) 
	      continue;
	  }

	  // Assign output offset & bin number
	  if (reduce_fac > 1) {
	    offset_out = (int32) ((offset/ratio) + 0.0);
	    //	    if ( offset_out == ncols_out)
	    // offset_out = ncols_out - 1;
	    bin_num_out = 
	      offset_out + output_binfile->get_basebin(irow/reduce_fac);
	  } else {
	    offset_out = offset;
	    bin_num_out = bin_num;
	  }

	  if (offset_out >= ncols_out) {
	    printf("Bad write to BINLIST: %d %d %d %d\n",
		   ifile, irow, ncols_out, offset_out);
	    exit(1);
	  }

	  if ( offset_out > max_out_kbin)
	    max_out_kbin = offset_out;

	  output_binfile->set_bin_num( offset_out, bin_num_out);
	  row_write = 1;

	  /* Sum & store number of observations,nscenes */
	  /* ------------------------------------------ */
	  int nobs = input_binfile[ifile]->get_nobs( kbin);
	  output_binfile->inc_nobs( offset_out, nobs);
	  int nscenes = input_binfile[ifile]->get_nscenes( kbin);
	  output_binfile->inc_nscenes( offset_out, nscenes);

	  /* Sum & store weights */
	  /* ------------------- */
	  if ( input.unit_wgt || input.median) {
	    output_binfile->set_weights( offset_out, 1);
	  } else {
	    output_binfile->inc_weights( offset_out, weights);
	  }

          output_binfile->inc_time_rec( offset_out, time_rec);

	  /* Product loop */
	  /* ------------ */
	  for (iprod=0; iprod<nprod[ifile]; iprod++) {
	    if (input_binfile[ifile]->active_data_prod[iprod] == true) {

	      if ( input.unit_wgt) {
		wgt = weights;
		f32 = in_sum_buf[iprod][2*kbin];
		out_sum_buf[iprod][2*offset_out] += f32/wgt;
		out_sum_buf[iprod][2*offset_out+1] += (f32/wgt)*(f32/wgt);
	      } else if ( input.median) {
		wgt = weights;
		f32 = in_sum_buf[iprod][2*kbin];
		sort_array[iprod][nfiles*offset+ifile] = f32/wgt;
	      } else {
		/* Add new sum to accumulated sum & sum2 */
		/* ------------------------------------- */
		f32 = in_sum_buf[iprod][2*kbin];
		out_sum_buf[iprod][2*offset_out] += f32;
		f32 = in_sum_buf[iprod][2*kbin+1];
		out_sum_buf[iprod][2*offset_out+1] += f32;
	      } /* input.unit_wgt */

	    } /* active_prod */

	  } /* product loop */
	  
	} /* kbin loop */

      } /* ========== if row has data ========== */

    } /* ifile loop */ 


    // Compute median
    if ( row_write && input.median) {
      for (kbin=0; kbin<numbin; kbin++) {
	for (iprod=0; iprod<nprod[0]; iprod++) {
	  if ( input_binfile[0]->active_data_prod[iprod] == true) {
	    float32 *sort_buf = (float32 *) calloc( nfiles, sizeof(float32));
	    int nsort = 0;
	    for (ifile=0; ifile<nfiles; ifile++) {
	      f32 = sort_array[iprod][nfiles*kbin+ifile];
	      if ( f32 != -999)
		sort_buf[nsort++] = sort_array[iprod][nfiles*kbin+ifile];
	    }
	    // Call insertion sort
	    if ( nsort > 0) {
	      insertion_sort(sort_buf, nsort);
	      out_sum_buf[iprod][2*(kbin/reduce_fac)] = sort_buf[nsort/2];
	      out_sum_buf[iprod][2*(kbin/reduce_fac)+1] = 1.0;
	      free( sort_buf);
	    }
	  }
	}
      }
    }

    /* Write output vdatas */
    /* ------------------- */
    if (row_write && (irow % reduce_fac) == reduce_fac-1) {
     
      n_write = 0;
      for (kbin=0; kbin<=max_out_kbin; kbin++) {

	int32 bin_num = output_binfile->get_bin_num( kbin);
	if (bin_num != 0) {

	  /* Loop over data products */
	  /* ----------------------- */
	  for (iprod=0; iprod<nprod[0]; iprod++) {
	    if (input_binfile[0]->active_data_prod[iprod] == true) {

	      /* Remove "blank" bin records */
	      /* -------------------------- */
	      if (n_write != kbin)
		memcpy(&out_sum_buf[iprod][2*n_write], 
		       &out_sum_buf[iprod][2*kbin], 8);
	    }
	  } /* iprod loop */

	  /* Remove "blank" bin records */
	  /* -------------------------- */
	  if (n_write != kbin) 
	    output_binfile->copy_binlist( kbin, n_write);

	  n_write++;
	  n_write_total++;

	} /* bin_num != 0 */
      } /* kbin loop */

      /* Write BinList & Data Products */
      /* ----------------------------- */
      if ( n_write > 0) {
	output_binfile->writeBinList( n_write);

	for (iprod=0; iprod<nprod[0]; iprod++) {
	  if (input_binfile[0]->active_data_prod[iprod] == true) {

	    input_binfile[0]->get_prodname( iprod, buf);
	    output_binfile->writeSums( &out_sum_buf[iprod][0], n_write, buf);
	  }
	}
	if ( strcmp(input.oformat,"HDF5") ==0 || 
             strcmp(input.oformat,"netCDF4") ==0 ) 
          output_binfile->incNumRec( n_write);
//	if ( isHDF5 || isCDF4) output_binfile->incNumRec( n_write);
      }

      /* Write to Quality Vdata */
      /* ---------------------- */
      i = 0;
      int32 nbin;
      if ( irow < nrows/2)
	nbin = input_binfile[0]->get_numbin(irow);
      else
	nbin = input_binfile[0]->get_numbin(irow-(reduce_fac/2));

      offmin = 
	(int32) ((minlon + 180) * (nbin / 360.0) + 0.5);
      offmax = 
	(int32) ((maxlon + 180) * (nbin / 360.0) + 0.5);

      for (kbin=0; kbin<ncols; kbin++) {
	if (kbin < offmin || kbin > offmax) continue;
	if (in_qual_buf[nfiles][kbin] != 255) {
	  in_qual_buf[nfiles][i++] = in_qual_buf[nfiles][kbin];
	}
      }

      if ( input_binfile[0]->has_qual()) {
	if ( (i-n_write) > 2) {
	  cout << "Problem with Quality write: irow: " << irow << " i: " <<
	    i << " n_write: " << n_write << " " << max_out_kbin << endl;
          exit(1);
	}

	output_binfile->writeQual( in_qual_buf[nfiles], n_write);
      } // has_qual() == true
    } /* row_write = 1 */


    // if median free storage arrays
    if ( input.median) {
      for (iprod=0; iprod<nprod[0]; iprod++) {
	if ( input_binfile[0]->active_data_prod[iprod] == true) {
	  free(sort_array[iprod]);
	}
      }
    }

  } /* irow loop (Main loop) */
  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

  for (ifile=0; ifile<nfiles; ifile++) {  
    free(in_qual_buf[ifile]);
  }
  free(uint8_buf);

  for (iprod=0; iprod<nprod[0]; iprod++) {
    if (input_binfile[0]->active_data_prod[iprod] == true) {
	free(in_sum_buf[iprod]);
	free(out_sum_buf[iprod]);
      }
  } /* iprod loop */

  // Copy metadata from input to output binfile
  output_binfile->copymeta( nfiles, &input_binfile[0]);
  int sensorID = sensorName2SensorId(output_binfile->meta_l3b.sensor_name);
  output_binfile->meta_l3b.sensorID = sensorID;

  // put in a fix to the missing Mission
  if(strlen(output_binfile->meta_l3b.mission) == 0) {
      if(sensorID != -1) {
          strcpy(output_binfile->meta_l3b.mission, platformName[sensorID]);
      }
  }

  strcpy(buf, input.ofile);
  strcpy( output_binfile->meta_l3b.product_name, buf);
  strcpy( output_binfile->meta_l3b.pversion, input.pversion);
  strcpy( output_binfile->meta_l3b.soft_name, "L3BIN");
  strcpy( output_binfile->meta_l3b.soft_ver, VERSION);
  if ( isCDF4 == 1) {
      // 1994-11-05T13:15:30Z
      strcpy(output_binfile->meta_l3b.ptime, unix2isodate(tnow, 'G'));
  } else {
      // yyyydddhhmmssmmm
      strcpy(output_binfile->meta_l3b.ptime, ydhmsf(tnow, 'L'));
  }
  strcpy( output_binfile->meta_l3b.ptime, ptime);
  strcpy( output_binfile->meta_l3b.proc_con, proc_con);
  strcpy( output_binfile->meta_l3b.input_parms, input.parms);

  char buf2[5000];
  if (Hishdf(input.infile) == TRUE || H5Fis_hdf5(input.infile) == TRUE) {
      strcpy( output_binfile->meta_l3b.infiles, basename(input.infile));

  } else {

    fp = fopen(input.infile, "r");
    buf2[0] = 0;
    for (ifile=0; ifile<nfiles; ifile++) {
      fgets(buf, 256, fp);
      buf[strlen(buf)-1] = 0;

      strcat(buf2, basename(buf));
      strcat(buf2, ",");
    } /* ifile loop */
    fclose(fp);
    buf2[strlen(buf2)-1] = 0;
    strcpy( output_binfile->meta_l3b.infiles, buf2);
  }

  for (ifile=0; ifile<nfiles; ifile++) input_binfile[ifile]->close();

  output_binfile->close();

  return 0;
}


/* Copyright (c) 2009 the authors listed at the following URL, and/or
the authors of referenced articles or incorporated external code:
http://en.literateprograms.org/Insertion_sort_(C)?action=history&offset=20081205204844

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Retrieved from: http://en.literateprograms.org/Insertion_sort_(C)?oldid=15530
*/

/* Sort an array of integers */
void insertion_sort(float a[], int length)
{
  int i;
  for (i=0; i < length; i++)
  {
      /* Insert a[i] into the sorted sublist */
    int j;
    float v = a[i];

      for (j = i - 1; j >= 0; j--)
      {
          if (a[j] <= v) break;
          a[j + 1] = a[j];
      }
      a[j + 1] = v;

  }
}

