/* ========================================================================
 * MSl1b2browse - converts any L1B file to a ppm file or browse format
 * 
 * Synopsis:
 *
 *   MSl1b2browse input-l1b-filename output-ppm-filename
 *
 * Description:
 * 
 * Modification history:
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Bryan A. Franz   GSC             Jul 1998     Original development
 *   Joel Gales	      Futuretech      Jan 2001     Added 8 and 24-bit hdf
 *   Norman Kuring    NASA/GSFC       Feb 2005     New scaling for true color
 *   Sean Bailey      Futuretech      Nov 2007     Made scaling consistent with
 *	                                               l1mapgen, use a sensor 
 *                                                 dependent defaults, and
 *                                                 generally clean up the code
 *   Joel Gales       Futuretech      Oct 2012     Define fftype to make 
 *                                                 compatible with new HDF/NCDF
 *                                                 libhdfutils
 *
 * ======================================================================== */

#include "l12_proto.h"
#include "input_struc.h"

#include "magick.h"
#include "image.h"
#include <png.h>

#define INT32   int32_t 
#define FLOAT32 float
#define BYTE    unsigned char

#define HDF8 0
#define HDF24 1
#define PPM 2
#define FLATBINARY 3
#define PNG 4

#define MALLOC(ptr,typ,num) {						\
  (ptr) = (typ *)malloc((num) * sizeof(typ));				\
  if((ptr) == NULL){							\
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",	\
    __FILE__,__LINE__);							\
    exit(EXIT_FAILURE);							\
  }									\
}

BYTE logscale (float val, float min, float max);
BYTE linscale (float val, float min, float max);
float toa_reflect(l1str * l1rec, int32_t ip, int32_t ib);

Image *initializeImage (int32 width, int32 height)
{
    Image *image;

    /* Allocate memory for an Image structure a la ImageMagick. */
    MALLOC (image, Image, 1);
    /* Initialize the Image structure.
     *  See AllocateImage() in image.c of the ImageMagick source code. */
    image->file = (FILE *) NULL;
    image->status = False;
    *image->filename = '\0';
    image->pipe = False;
    (void) strcpy (image->magick, "MIFF");
    image->comments = (char *) NULL;
    image->label = (char *) NULL;
    image->id = UndefinedId;
    image->class = DirectClass;
    image->matte = False;
    image->compression = RunlengthEncodedCompression;
    image->columns = width;
    image->rows = height;
    image->colors = 0;
    image->scene = 0;
    image->montage = (char *) NULL;
    image->directory = (char *) NULL;
    image->colormap = (ColorPacket *) NULL;
    image->signature = (char *) NULL;
    MALLOC (image->pixels, RunlengthPacket, width * height);
    image->packet = (RunlengthPacket *) NULL;
    image->packets = width * height;
    image->packet_size = 0;
    image->packed_pixels = (unsigned char *) NULL;
    image->orphan = False;
    image->previous = (Image *) NULL;
    image->next = (Image *) NULL;
    return (image);
}


void dfr8_addimage (char *fname,
                    unsigned char *raster,
                    int32 width,
                    int32 height, unsigned char *palette, char *label)
{
    uint16 ref;

    if (DFR8setpalette (palette) == FAIL) {
        fprintf (stderr,
                 "-E- %s line %d: DFR8setpalette() failed for file, %s.\n",
                 __FILE__, __LINE__, fname);
        exit (EXIT_FAILURE);
    }
    if (DFR8addimage (fname, raster, width, height, COMP_NONE) == FAIL) {
        fprintf (stderr,
                 "-E- %s line %d: DFR8addimage() failed for file, %s.\n",
                 __FILE__, __LINE__, fname);
        exit (EXIT_FAILURE);
    }
    ref = DFR8lastref ();
    if (ref == 0) {
        fprintf (stderr,
                 "-E- %s line %d: DFR8lastref() failed for file, %s.\n",
                 __FILE__, __LINE__, fname);
        exit (EXIT_FAILURE);
    }
    if (DFANputlabel (fname, DFTAG_RIG, ref, "brs_data")) {
        fprintf (stderr,
                 "-E- %s line %d: DFANputlabel() failed for file, %s.\n",
                 __FILE__, __LINE__, fname);
        exit (EXIT_FAILURE);
    }
}

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main (int argc, char *argv[])
{

    int32_t iscan = 0;             /* input scan number                  */
    int32_t spix = 0;              /* start pixel for subscene process   */
    int32_t epix = -1;             /* end pixel for subscene process     */
    int32_t dpix = 1;              /* pixel increment for sub-sampling   */
    int32_t sscan = 0;             /* start scan for subscene process    */
    int32_t escan = -1;            /* end scan for subscene process      */
    int32_t dscan = 1;             /* scan subsampling increment         */
    int32_t ip;                    /* input pixel number                 */
    int32_t iw;                    /* wavelength  number                 */
    int32_t sfactor = 10;
    double esdist;
    static float *Fobar;

    float  rhoc = 0.0;             /* cirrus reflectance                 */
    float  twc  = 0.8;             /* water vapor transmittance above cirrus */

    int32_t npix = 0;              /* Number of output pixels per scan   */
    int32_t nscan = 0;             /* Number of output scans             */

    BYTE rgb[3];
    int32_t r, g, b;
    int32_t *wave;
    int want_bin = 0;
    int want_ppm = 0;
    int want_hdf = 0;
    int want_8bit = 0;
    int want_24bit = 0;
    int want_atmocor = 0;
    int want_linscale = 0;
    int want_png =0;
    
    float sr_r, sr_g, sr_b;
    float min = 0.01;
    float max = 0.9;

    l1str l1rec;                /* generic level-1b scan structure      */
    filehandle l1file;          /* input file handle                    */
    FILE *outfp = NULL;
    int32_t hdr[5];
    int32_t length;

    char outfile[FILENAME_MAX] = "";
    char rflag[FILENAME_MAX] = "ORIGINAL";

    instr input;                /* input parameters structure         */
    int32 dims[8];
    int32 rank;
//    int32 HDFfid_w;
    int32 sd_id_w;
    int32 sds_id_rgb;
    int32 sds_id_lon;
    int32 sds_id_lat;
    int32 start_w[3] = { 0, 0, 0 };
    int32 count_w[3] = { 1, 1, 3 };
    float32 f32;
    char strbuf[80];
    char strbuf2[80];

    Image *img;
    RunlengthPacket *pix;
    ColorPacket *clr_ptr;

    unsigned char palette[768], *lut;
    unsigned char *raster, *ras_ptr;
    int32 numoutpix;
    int32_t nbands;

    int i, j, k;

    int16 *lonBuffer = NULL;
    int16 *latBuffer = NULL;
    BYTE *lineBuffer = NULL;

    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
        
    
    /* hold all of the command line options */
    clo_optionList_t* list;

    if (argc == 1) {
        l2gen_usage("l1brsgen");
        return 0;
    }

    // see if help on command line
    for(i=0; i<argc; i++) {
        if((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "-help") == 0)) {
            l2gen_usage("l1brsgen");
            return 1;
        }
    }

    cdata_();
    filehandle_init (&l1file);
    msl12_input_init(&input);

     // create an option list
    list = clo_createList();

    /* initialize the option list with descriptions and default values */
    l2gen_init_options(list, "l1brsgen");

    // change the default values for this program
    clo_setString(list, "atmocor", "off", "default");
    clo_setString(list, "proc_land", "1", "default");
    clo_setString(list, "sl_pixl", "0", "default");

    // Parse input parameters
    if (msl12_option_input(argc, argv, list, "l1brsgen", &input, &l1file) != 0) {
        printf ("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit (FATAL_ERROR);
    }
    strcpy (outfile, input.ofile[0]);

    want_linscale = input.stype;
    if (input.atmocor == 1)
        want_atmocor = 1;


    /*                                                                  */
    /* Open input file and get sensor and scan information from handle  */
    /*                                                                  */
    if (openl1 (&l1file) != 0) {
        printf ("-E- %s: Error opening %s for reading.\n", argv[0],
                l1file.name);
        exit (1);
    }

    r = bindex_get (input.rgb[0]);
    g = bindex_get (input.rgb[1]);
    b = bindex_get (input.rgb[2]);
    if (r < 0 || g < 0 || b < 0) {
        printf ("-E- Invalid RGB set for this sensor\n");
        if (r < 0)
            printf ("  Invalid R: %d\n", input.rgb[0]);
        if (g < 0)
            printf ("  Invalid G: %d\n", input.rgb[1]);
        if (b < 0)
            printf ("  Invalid B: %d\n", input.rgb[2]);
        exit (FATAL_ERROR);
    }

    min = input.datamin;
    max = input.datamax;


    /*                                                                  */
    /* Allocate memory for L1 scan data                                 */
    /*                                                                  */
    nbands = input.nbands;
    if (alloc_l1 (l1file.npix, nbands, NBANDSIR, l1file.n_inprods, &l1rec) ==
        0) {
        printf ("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit (FATAL_ERROR);
    }

    /* Set the end pixel if it was not set by command argument	        */
    if (input.epixl < 1 || input.epixl > l1file.npix)
        input.epixl = l1file.npix;
    if (input.eline < 1 || input.eline > l1file.nscan)
        input.eline = l1file.nscan;
    if (input.spixl < 1)
        input.spixl = 1;
    if (input.sline < 1)
        input.sline = 1;
    sfactor = MAX(input.subsamp, 1);
    spix  = MAX(input.spixl - 1, 0);
    epix  = MIN(input.epixl - 1, l1file.npix-1);
    sscan = MAX(input.sline - 1, 0);
    escan = MIN(input.eline - 1, l1file.nscan-1);
    dpix  = sfactor;
    dscan = sfactor;

    npix  = (epix - spix) / sfactor + 1;
    nscan = (escan - sscan) / sfactor + 1;

    l1file.spix = spix;
    l1file.epix = epix;
    l1file.dpix = dpix;

    input.sline = sscan + 1;
    input.eline = escan + 1;
    input.dline = dscan;
    input.spixl = spix + 1;
    input.epixl = epix + 1;
    input.dpixl = dpix;

    printf ("Using r,g,b = %d, %d, %d\n", input.rgb[0], input.rgb[1],
            input.rgb[2]);

    if(strcmp(input.oformat, "HDF4") == 0) {
        if(strcmp(input.oformat_depth, "8bit") == 0)
            want_8bit = 1;
        else
            want_24bit = 1;
    } else if(strcmp(input.oformat, "BIN") == 0)
        want_bin = 1;
    else if(strcmp(input.oformat, "PNG") == 0)
        want_png = 1;
    else if(strcmp(input.oformat, "PPM") == 0)
        want_ppm = 1;

    if (want_ppm + want_24bit + want_8bit + want_bin + want_png > 1) {
        printf ("%s: Error: More than one output format has been chosen.\n",argv[0]);
        exit (1);
    }
    if (want_ppm + want_24bit + want_8bit + want_bin + want_png < 1) {
        want_8bit = 1;
    }
    if (want_8bit || want_24bit) {
        want_hdf = 1;
    }

    if (want_ppm) {

        if ((outfp = fopen (outfile, "w")) == NULL) {
            printf ("%s: Error: Unable to open %s for writing.\n",
                    argv[0], outfile);
            exit (FATAL_ERROR);
        }
        /* Write output ppm file */
        fprintf (outfp, "P6\n");
        fprintf (outfp, "%d\n", npix);
        fprintf (outfp, "%d\n", nscan);
        fprintf (outfp, "255\n");

    } else if (want_hdf) {

//        HDFfid_w = Hopen (outfile, DFACC_CREATE, 0);
//        status = Vstart (HDFfid_w);
//        sd_id_w = SDstart (outfile, DFACC_RDWR);
        sd_id_w = SDstart (outfile, DFACC_CREATE);

        rank = 2;
        dims[0] = nscan;
        dims[1] = npix;
        sds_id_lon = SDcreate (sd_id_w, "longitude", DFNT_INT16, rank, dims);

        SDsetdimname (SDgetdimid (sds_id_lon, 0), "Number of Scans");
        SDsetdimname (SDgetdimid (sds_id_lon, 1), "Pixels per Scan");

        f32 = 1. / 180;
        SDsetattr (sds_id_lon, "slope", DFNT_FLOAT32, 1, &f32);
        f32 = 0.0;
        SDsetattr (sds_id_lon, "intercept", DFNT_FLOAT32, 1, &f32);

        rank = 2;
        dims[0] = nscan;
        dims[1] = npix;
        sds_id_lat = SDcreate (sd_id_w, "latitude", DFNT_INT16, rank, dims);

        SDsetdimname (SDgetdimid (sds_id_lat, 0), "Number of Scans");
        SDsetdimname (SDgetdimid (sds_id_lat, 1), "Pixels per Scan");

        f32 = 1. / 360;
        SDsetattr (sds_id_lat, "slope", DFNT_FLOAT32, 1, &f32);
        f32 = 0.0;
        SDsetattr (sds_id_lat, "intercept", DFNT_FLOAT32, 1, &f32);

        lonBuffer = (int16*) malloc(npix * 2);
        latBuffer = (int16*) malloc(npix * 2);

        count_w[1] = npix;

        if (want_24bit) {
            rank = 3;
            dims[0] = nscan;
            dims[1] = npix;
            dims[2] = 3;
            sds_id_rgb = SDcreate (sd_id_w, "rgb", DFNT_UINT8, rank, dims);

            SDsetdimname (SDgetdimid (sds_id_rgb, 0), "Number of Scans");
            SDsetdimname (SDgetdimid (sds_id_rgb, 1), "Pixels per Scan");
            SDsetdimname (SDgetdimid (sds_id_rgb, 2), "Number of Colors");
            
            lineBuffer = (unsigned char*) malloc(npix * 3);
        }


        if (want_8bit) {
            /* Initialize an Image structure to be used in color quantization. */
            img = initializeImage (npix, nscan);
            pix = img->pixels;
        }

    } else  if (want_bin) {

        if ((outfp = fopen (outfile, "w")) == NULL) {
            printf ("%s: Error: Unable to open %s for writing.\n",
                    argv[0], outfile);
            exit (FATAL_ERROR);
        }

        /* Write file header */
        length = 12 + npix * 2 * 4 + npix * 3;
        hdr[0] = length;
        hdr[1] = npix;
        hdr[2] = nscan;
        hdr[3] = l1file.sensorID;

        fwrite (hdr, sizeof (hdr), 1, outfp);
        fseek (outfp, length, SEEK_SET);

    } else if (want_png) {

        if ((outfp = fopen (outfile, "w")) == NULL) {
            printf ("%s: Error: Unable to open %s for writing.\n",
                    argv[0], outfile);
            exit (FATAL_ERROR);
        }

        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                NULL, NULL, NULL);
        if (!png_ptr) {
            fprintf(stderr, "%s: Error: Unable to create PNG write structure.\n",argv[0]);
            exit(1);
        }

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "%s: Error: Unable to create PNG info structure.\n",argv[0]);
            exit(1);
        }
        if (setjmp(png_jmpbuf(png_ptr))) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "%s: Error: Unable to call PNG setjmp().\n",argv[0]);
            exit(1);
        }
        png_init_io(png_ptr, outfp);
        png_set_IHDR(png_ptr, info_ptr, npix, nscan,
                8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(png_ptr, info_ptr);

        lineBuffer = (unsigned char*) malloc(npix * 3);
        if (!lineBuffer) {
            png_destroy_write_struct(&png_ptr, (png_infopp) NULL);
            fprintf(stderr, "%s: Error: Unable to allocate line buffer.\n",argv[0]);
            exit(1);
        }
        
        
        
    } else {
        printf ("%s: Error: No output format specified.\n", argv[0]);
        exit (FATAL_ERROR);
        
    }
    
    /*                                                          */
    /* Read file scan by scan, scale radiances,  and write.     */
    /*                                                          */
    rdsensorinfo (l1file.sensorID, 0, "Fobar", (void **) &Fobar);
    for (iscan = sscan; iscan <= escan; iscan += dscan) {

        if (readl1 (&l1file, iscan, &l1rec) != 0) {
            fprintf (stderr,
                     "-E- %s Line %d: error reading %s at scan %d.\n",
                     __FILE__, __LINE__, l1file.name, iscan);
            exit (FATAL_ERROR);
        }


        if (want_bin) {
            fwrite (l1rec.year, 4, 1, outfp);
            fwrite (l1rec.day, 4, 1, outfp);
            fwrite (l1rec.msec, 4, 1, outfp);

            for (ip = 0; ip < npix; ip++) {
                fwrite (&l1rec.lon[ip], 4, 1, outfp);
            }
            for (ip = 0; ip < npix; ip++) {
                fwrite (&l1rec.lat[ip], 4, 1, outfp);
            }

        } else if (want_hdf) {
            start_w[0] = (iscan - sscan) / sfactor;
            for (ip = 0; ip < npix; ip++) {
                if (l1rec.lon[ip] >= 0)
                    lonBuffer[ip] = (int16) ((l1rec.lon[ip] * 180) + 0.5);
                if (l1rec.lon[ip] < 0)
                    lonBuffer[ip] = (int16) ((l1rec.lon[ip] * 180) - 0.5);
 
                if (l1rec.lat[ip] >= 0)
                    latBuffer[ip] = (int16) ((l1rec.lat[ip] * 360) + 0.5);
                if (l1rec.lat[ip] < 0)
                    latBuffer[ip] = (int16) ((l1rec.lat[ip] * 360) - 0.5);
            }
            
            SDwritedata (sds_id_lon, start_w, NULL, count_w,
                                  (VOIDP) lonBuffer);

            SDwritedata (sds_id_lat, start_w, NULL, count_w,
                                  (VOIDP) latBuffer);

        }

        if (want_atmocor) {
            if (loadl1 (&l1file, &input, &l1rec) != 0) {
                fprintf (stderr,
                         "-E- %s Line %d: error loading %s at scan %d.\n",
                         __FILE__, __LINE__, l1file.name, iscan);
                exit (FATAL_ERROR);
            }
        } else {
            /* Get correction for Earth-Sun distance and apply to Fo  */
            esdist = esdist_ (l1rec.year, l1rec.day, l1rec.msec);
            l1rec.fsol = pow (1.0 / esdist, 2);
            for (iw = 0; iw < l1rec.nbands; iw++) {
                l1rec.Fobar[iw] = Fobar[iw];
                l1rec.Fo[iw] = l1rec.Fobar[iw] * l1rec.fsol;
            }
        }

        for (ip = 0; ip < npix; ip++) {
            if (l1rec.Lt[ip * l1rec.nbands + r] > 0.01) {
                if (want_atmocor) {
                    if (input.cirrus_opt)
                        rhoc = l1rec.rho_cirrus[ip]/twc;
                    else
                        rhoc = 0.0;

                    sr_r = l1rec.rhos[ip * l1rec.nbands + r] - rhoc;
                    sr_g = l1rec.rhos[ip * l1rec.nbands + g] - rhoc;
                    sr_b = l1rec.rhos[ip * l1rec.nbands + b] - rhoc;
                    /*
                       printf("r %d %f %f %f\n",ip,l1rec.Lt[ip*NBANDS+r],l1rec.Lr[ip*NBANDS+r],l1rec.rhos[ip*NBANDS+r]);
                       printf("g %d %f %f %f\n",ip,l1rec.Lt[ip*NBANDS+g],l1rec.Lr[ip*NBANDS+g],l1rec.rhos[ip*NBANDS+g]);
                       printf("b %d %f %f %f\n",ip,l1rec.Lt[ip*NBANDS+b],l1rec.Lr[ip*NBANDS+b],l1rec.rhos[ip*NBANDS+b]);
                     */
                } else {
                    sr_r = toa_reflect (&l1rec, ip, r);
                    sr_g = toa_reflect (&l1rec, ip, g);
                    sr_b = toa_reflect (&l1rec, ip, b);
                }
            if (want_linscale == 1) {
                rgb[0] = linscale (sr_r, min, max);
                rgb[1] = linscale (sr_g, min, max);
                rgb[2] = linscale (sr_b, min, max);
            } else {
                rgb[0] = logscale (sr_r, min, max);
                rgb[1] = logscale (sr_g, min, max);
                rgb[2] = logscale (sr_b, min, max);
            }
        } else {
            rgb[0] = rgb[1] = rgb[2] = 255;
        }
        if (l1file.sensorID == OCTS && rgb[2] > 234) {
            if (rgb[0] < 230)
                rgb[0] = 250;
            if (rgb[1] < 230)
                rgb[1] = 251;
            rgb[2] = 252;
        }


            if (want_bin || want_ppm) {
                fwrite (rgb, 1, 3, outfp);
            } else if(want_hdf) {

                if (want_24bit) {
                    i = ip * 3;
                    lineBuffer[i++] = rgb[0];
                    lineBuffer[i++] = rgb[1];
                    lineBuffer[i]   = rgb[2];
                } else if (want_8bit) {
                    pix->red = rgb[0];
                    pix->green = rgb[1];
                    pix->blue = rgb[2];
                    pix->index = pix->length = 0;
                    pix++;
                }

            } else if (want_png) {
                i = ip * 3;
                lineBuffer[i++] = rgb[0];
                lineBuffer[i++] = rgb[1];
                lineBuffer[i]   = rgb[2];
            }

        } // for ip

        if (want_24bit) {
            SDwritedata (sds_id_rgb, start_w, NULL, count_w,
                (VOIDP) lineBuffer);

        } else if (want_png) {
            png_write_row(png_ptr, lineBuffer);
            
        }

    } // for iscan


    /* Write Global Attributes */

    if (want_hdf) {
      idDS ds_id;
      ds_id.deflate = 0;
      ds_id.fid = sd_id_w;
      ds_id.fftype = DS_HDF;
        SetChrGA (ds_id, "Product Name", outfile);

        sprintf (strbuf, "%s Level-1 Browse Data",
                 sensorName[l1file.sensorID]);
        SetChrGA (ds_id, "Title", strbuf);

/*
        sprintf (strbuf,
                 "NASA/GSFC %s Level-1 %s quasi-true-color browse data, day %d, %d",
                 sensorName[l1file.sensorID], strbuf2, *l1rec.day,
                 *l1rec.year);
*/
        sprintf (strbuf,
                 "NASA/GSFC %s Level-1 quasi-true-color browse data, day %d, %d",
                 sensorName[l1file.sensorID], *l1rec.day, *l1rec.year);
        SetChrGA (ds_id, "Legend", strbuf);

        SetChrGA (ds_id, "Processing Time", ydhmsf (now (), 'G'));
        SDsetattr (sd_id_w, "Number of Scan Lines", DFNT_INT32, 1,
                   (VOIDP) & nscan);
        SDsetattr (sd_id_w, "Pixels per Scan Line", DFNT_INT32, 1,
                   (VOIDP) & npix);

        scene_meta_write (ds_id);
    }

    closel1 (&l1file);

    if (want_bin || want_ppm) {
        fclose (outfp);
    } else if(want_hdf) {

        if (want_24bit)
            SDendaccess (sds_id_rgb);

        if (want_8bit) {
            QuantizeImage (img, 256, 0, 0, YIQColorspace, True);

            /* Extract the palette from the image structure. */
            lut = palette;
            clr_ptr = img->colormap;
            for (i = 0; i < img->colors; i++) {
                *lut++ = clr_ptr->red;
                *lut++ = clr_ptr->green;
                *lut++ = clr_ptr->blue;
                clr_ptr++;
            }
            /* Pad the palette to 256 colors if QuantizeImage returned fewer than 256. */
            for (; i < 256; i++) {
                *lut++ = 0;
                *lut++ = 0;
                *lut++ = 0;
            }

            /* Extract the 8-bit raster from the image structure. */
            numoutpix = nscan * npix;
            MALLOC (raster, unsigned char, numoutpix);

            ras_ptr = raster;
            pix = img->pixels;
            k = 0;
            for (i = 0; i < img->packets; i++) {
                for (j = 0; j <= (int) pix->length; j++) {
                    if (++k > numoutpix) {
                        fprintf (stderr, "-E- %s line %d: Too many pixels\n",
                                 __FILE__, __LINE__);
                        exit (EXIT_FAILURE);
                    }
                    *ras_ptr++ = (unsigned char) pix->index;
                }
                pix++;
            }

            /* I don't need the Image structures anymore. */
            free (img->pixels);
            free (img);

            /* Write the image and its palette to the output file. */
            dfr8_addimage (outfile, raster, npix, nscan, palette, "brs_data");
            free (raster);

        }

        SDendaccess (sds_id_lon);
        SDendaccess (sds_id_lat);
        SDend (sd_id_w);
//        Hclose (HDFfid_w);
//        Vend (HDFfid_w);
    } else if(want_png) {
        png_write_end(png_ptr, info_ptr);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        
        
        fclose (outfp);
    }

    exit (0);
}
