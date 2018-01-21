/* ========================================================================
 *
 *   MSl1b2info input-l1b-filename
 *
 * Description:
 * 
 * Modification history:
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Bryan A. Franz   GSC             28 March 1999 Original development
 *   Dan Knowles      SAIC            2004          Complete rewrite to create geoboxes
 *   W. Robinson      SAIC            10 dec 2004   add CZCS sensor
 *   Dan Knowles      SAIC            18 Aug 2005   Changed daynight logic to be based on solar zenith angle
 *   Sean Bailey      NASA            23 Oct 2015   Added logic to determine if sun in SD for VIIRS
 *
 * ======================================================================== */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>

#include "l12_proto.h"
#include <timeutils.h>
#include "version.h"

#define NUM_SPACECRAFT_DIRECTION_MODES 3
#define DEFAULT_SPACECRAFT_DIRECTION 0   /* set to an invalid value */
#define ASCENDING_TRACK  1
#define DESCENDING_TRACK 2

#define NORTH 1
#define SOUTH 0

#ifdef TRUE
#undef TRUE
#endif

#define TRUE 1
#define FALSE 0
#define SECTOR_ROTATION 4

#define TMP_FILENAME_MAX 255

#define DEFAULT_COORD_VALUE -999.   /* set to an invalid value */

#define DEFAULT_DAYNIGHT_MODE 0
#define DAY_MODE (1<<0)
#define NIGHT_MODE (1<<1)

/* equatorial radius, didn't bother to use polar value also since precision not that important for current use */
#define EARTH_RADIUS_EQUATORIAL 6378   

typedef struct {
    float32 north_lat;
    float32 south_lat;
    float32 west_lon;
    float32 east_lon;
    unsigned char daynightflag;
    uint32_t pixel_count;
} box_t;

typedef struct {
    int sensor_id;
    int day_node;
} day_node_t;


void set_north_south_boundaries(float32, float32 *, float32 *);
void set_west_east_boundaries(float32, float32 *, float32 *);
int check_if_in_west_east_boundaries(float32, float32, float32);
void get_box_square_km(double, double, double, double, double *);
void get_box_square_km_sublevel(double, double, double, double, double *);
int check_if_in_box(float32, float32, float32, float32, float32, float32);
double get_lon_distance(double, double);

#define USAGESTR \
"%s %d.%d.%d-r%d (%s %s)\n\n\
Usage: %s [-n number-of-boxes] [-s] L1-filename [met-file]\n\n\
Where:\n\
\t-d degree limits to set on boxes (use instead of -n)\n\
\t-i set subsampling increment (default=1 gets all scan lines and pixels)\n\
\t-n generates up to this number of lat/lon boxes per orbit-side\n\
\t-s prints only data needed for database fields\n\
\t-v verify number of pixels that match in boxed region\n\
\t-o [output file]\n\n\
Return status codes:\n\
\t0 - good\n\
\t1 - general error\n\
\t100 - all pixels have NAVFAIL or NAVWARN flag set\n\
\t101 - no 'good' navigation was found\n\
\t102 - direction (ascending/descending) not set\n\
\t103 - master day night flag not set\n\
\t104 - total box count = 0\n\
(For multiple conditions, lowest value status code takes priority)\n\
\n\
"


#define PRINT_USAGE(x)  printf(USAGESTR,(x),L2GEN_VERSION_MAJOR,L2GEN_VERSION_MINOR,L2GEN_VERSION_PATCH_LEVEL,SVN_REVISION,__DATE__,__TIME__,(x))

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    l1str l1rec; /* generic level-1b scan structure      */
    filehandle l1file; /* input file handle                    */
    instr input;
    float64 utime;
    char ofile[FILENAME_MAX];
    FILE *fp;
    char *isodatetime;
#if 0
    char geofile[FILENAME_MAX] = "";
    char geofile[TMP_FILENAME_MAX] = "";
#endif
    int16 year, month, day, hour, minute;
    double second;
    float32 mem_st_lat, mem_st_lon, mem_en_lat, mem_en_lon;

    int32_t spix = 0; /* start pixel for subscene process   */
    int32_t cpix = -1; /* center pixel of subscene scan      */
    int32_t epix = -1; /* end pixel for subscene process     */
    int32_t sscan = 0; /* start scan for subscene process    */
    int32_t cscan = -1; /* center scan for subscene process   */
    int32_t escan = -1; /* end scan for subscene process      */
    int32_t iscan;
    int32_t ip;
    int32_t num_match = 0;
    int32_t num_no_match = 0;
    int32_t num_nav_fail = 0;
    int32_t num_bad_lat_lon_pixels = 0;
    int32_t num_pixels = 0;
    int32_t navwarn_all_set = 1;
    int32_t good_spix, good_epix, good_cpix;

    int i;
    int box_index;
    int spacecraft_direction_index;
    int curr_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int prev_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int initial_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int final_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
    int match;
    int in_box;
    int num_boxes = 0;
    int num_degrees = 0;
    int show_diagnostic_info;
    int show_standard_info;
    int show_sdps_info;
    int case_n = FALSE;
    int case_s = FALSE;
    int case_d = FALSE;
    int case_v = FALSE;
    int case_o = FALSE;
    int total_box_count = 0;
    int user_defined_increment = 1;
    int exitflag = 0;
    int line_nav_ok;
    int start_printed;
    int np;
    static int first_good_scan = -1;

    double total_box_square_km = 0;
    double box_square_km;
    double percent_match;

    unsigned char curr_daynight_mode = DEFAULT_DAYNIGHT_MODE;
    unsigned char master_daynightflag = DEFAULT_DAYNIGHT_MODE;

    static char *node[3] = {"Unidentified", "Ascending", "Descending"};
    static char *daynight_string[4] = {"Unidentified", "Day", "Night", "Both"};

    float32 epix_northern_lat = DEFAULT_COORD_VALUE;
    float32 epix_southern_lat = DEFAULT_COORD_VALUE;
    float32 spix_northern_lat = DEFAULT_COORD_VALUE;
    float32 spix_southern_lat = DEFAULT_COORD_VALUE;
    float32 northern_lat = DEFAULT_COORD_VALUE;
    float32 southern_lat = DEFAULT_COORD_VALUE;
    float32 eastern_lon = DEFAULT_COORD_VALUE;
    float32 western_lon = DEFAULT_COORD_VALUE;
    float32 center_lon = DEFAULT_COORD_VALUE;
    float32 center_lat = DEFAULT_COORD_VALUE;
    float32 prev_lat_cpix = DEFAULT_COORD_VALUE;
    float32 tmp_north_lat;
    float32 tmp_south_lat;
    float32 tmp_west_lon;
    float32 tmp_east_lon;
    float32 lat_breakpoint;
    float32 increment;


    box_t **box = NULL; /* will become 2 dimensional box array in [boxes][direction]  */
    int counter;
    int num_args;

    char *option_string = "d:hi:n:sv:o:";
    int options = 0;
    while ((options = getopt(argc, argv, option_string)) != -1) {
        switch (options) {
            case 'd':
                case_d = TRUE;
                num_degrees = atoi(optarg);
                break;
            case 'i':
                user_defined_increment = atoi(optarg);
                break;
            case 'h':
                PRINT_USAGE(argv[0]);
                exit(FATAL_ERROR);
                break;
            case 'n':
                case_n = TRUE;
                num_boxes = atoi(optarg);
                break;
            case 's':
                case_s = TRUE;
                break;
            case 'v':
                case_v = TRUE;
                break;
            case 'o':
                case_o = TRUE;
                snprintf(ofile, FILENAME_MAX, "%s", optarg);
                break;
            default:
                break;
        }
    }


    /*******************************************************************************************
     *    Check options and set num_boxes
     *******************************************************************************************/

    if (case_d == TRUE) {
        if (num_degrees >= 1 && num_degrees <= 180) {
            if (case_n != TRUE) {
                num_boxes = (int) 180 / num_degrees;
                printf("num_boxes=%d\n", num_boxes);
            } else {
                printf("INPUT ERROR: Cannot choose both -d and -n options\n");
                PRINT_USAGE(argv[0]);
                exit(FATAL_ERROR);
            }
        } else {
            printf("INPUT ERROR: -d (degrees lat) option must be between 1 and 180\n");
            PRINT_USAGE(argv[0]);
            exit(FATAL_ERROR);
        }
    }
    else if (case_n == TRUE) {
        if (num_boxes < 0) {
            printf("INPUT ERROR: -n (number of lat divisions) option must be greater than or equal to 0\n");
            PRINT_USAGE(argv[0]);
            exit(FATAL_ERROR);
        }
    } else {
        num_boxes = 0;
    }
    /*    else if (case_s == TRUE)
        {
            if (case_n != TRUE && case_d != TRUE)
            {
                printf("INPUT ERROR: -s option requires -d or -n option also\n");
                PRINT_USAGE(argv[0]);
                exit(FATAL_ERROR);
            }
        }
     */

    /*******************************************************************************************
     *    Setup print options
     *******************************************************************************************/

    if (case_s == TRUE) {
        show_sdps_info = TRUE;
        show_standard_info = FALSE;
        show_diagnostic_info = FALSE;
    }
    else if (case_n == TRUE || case_d == TRUE) {
        show_sdps_info = FALSE;
        show_standard_info = TRUE;
        show_diagnostic_info = TRUE;
    }
    else {
        show_sdps_info = FALSE;
        show_standard_info = TRUE;
        show_diagnostic_info = FALSE;
    }
    if (case_o == TRUE) {
        fp = fopen(ofile, "w");
    } else {
        fp = stdout;
    }


    num_args = argc - optind;

    // make the l1 readers quiet
    want_verbose = 0;

    cdata_();
    msl12_input_init(&input);
    filehandle_init(&l1file);
    l1file.input = &input;

    /*									*/
    /* Process required command-line arguments                          */
    /*									*/
    switch (num_args) {
        case 2:
            strcpy(l1file.name, argv[optind + 0]);
            /* strcpy(geofile,argv[optind+1]);*/
            l1file.geofile = argv[optind + 1];
            break;
        case 1:
            strcpy(l1file.name, argv[optind + 0]);
            break;
        default:
            PRINT_USAGE(argv[0]);
            exit(FATAL_ERROR);
            break;
    }

    /* Set default input parameters */
    if (msl12_input_defaults(&l1file, &input) != 0) {
        printf("-E- %s: Error parsing input parameters.\n", argv[0]);
        exit(FATAL_ERROR);
    }
    strcpy(input.program_name, "l1info");
    /* for VIIRS, change default to no cal (avoids trouble on previously
       calibrated files */
    if (input.sensorID == VIIRS) strcpy(input.calfile, "");
    input.proc_sst = 0;
    input.proc_ocean = 0;
    input.atmocor = 0;
    input.proc_land = 0;
    
    // The bowtie effect...if an extract is small enough and on the scene
    // edge, the bowtie effect can confuse the geobox logic, so bump the
    // subsampling increment to make things happy
    if ((input.sensorID == HMODISA || input.sensorID == HMODIST) &&
        (user_defined_increment < 10)) {
        user_defined_increment = 10;
    }

    /*
     * The following tidbit calculates if the Sun is visible to the solar
     * diffuser for VIIRS so that the OBC files can be flagged to keep
     * Since this requires reading from the geolocation file, this bit is
     * done prior to opening the L1 file (which would also open the geo file)
     */
    int sun_in_sd = FALSE;
    int sector_rotation = FALSE;
    int filledScans = 0;
    if (input.sensorID == VIIRS) {
        if (user_defined_increment < 16)
            user_defined_increment = 16;
        // Get the ACTUAL scans for printing the real number not the 
        // arbitrary dimension size
        int32_t ncid, status;
        if ((nc_open(input.ifile[0], NC_NOWRITE, &ncid)) == NC_NOERR) {
            nc_get_att_int(ncid, NC_GLOBAL, "number_of_filled_scans", &filledScans);
            status = nc_close(ncid);
            if (status != NC_NOERR) {
                printf("-E- %s: Error closing %s.\n", argv[0], input.ifile[0]);
                exit(FATAL_ERROR);
            };
        }
        int32_t geoid;
        //Only bother to do this if provided a geolocation file
        if ((nc_open(input.geofile, NC_NOWRITE, &geoid)) == NC_NOERR) {

            int32_t grpid, varid, dimid;
            size_t nscans, nvectors;
                
            // don't try this for the SDR files...only the NASA L1 format files
            if ((nc_inq_ncid(geoid, "All_Data", &grpid)) != NC_NOERR) {
                status = nc_inq_dimid(geoid, "number_of_scans", &dimid);
                if (status != NC_NOERR) {
                    printf("-E- Error reading number_of_scans attribute.\n");
                    exit(FATAL_ERROR);
                };
                nc_inq_dimlen(geoid, dimid, &nscans);
                status = nc_inq_dimid(geoid, "vector_elements", &dimid);
                if (status != NC_NOERR) {
                    printf("-E- Error reading vector_elements attribute.\n");
                    exit(FATAL_ERROR);
                };
                nc_inq_dimlen(geoid, dimid, &nvectors);

                if ((nc_inq_ncid(geoid, "navigation_data", &grpid)) == NC_NOERR) {
                    status = nc_inq_varid(grpid, "solar_inst", &varid);
                    if (status != NC_NOERR) {
                        printf("-E- Error finding solar_inst variable.\n");
                        exit(FATAL_ERROR);
                    };

                    float **solar_inst; //(number_of_scans, vector_elements);
                    if ((solar_inst = malloc(nscans * sizeof (float *))) == NULL) {
                        printf("%s -Error: Cannot allocate memory to solar_inst data\n",
                                __FILE__);
                        exit(FATAL_ERROR);
                    }
                    float *tmp_arr;
                    if ((tmp_arr = malloc(nscans * nvectors * sizeof (float))) == NULL) {
                        printf("%s -Error: Cannot allocate memory to solar_inst data\n",
                                __FILE__);
                        exit(FATAL_ERROR);
                    }
                    int i;
                    for (i = 0; i < nscans; i++) {
                        solar_inst[i] = tmp_arr + i *nvectors;

                    }
                    size_t start[] = {0, 0}; /* start at first value */
                    size_t cnt[] = {nscans, nvectors};
                    status = nc_get_vara_float(grpid, varid, start, cnt, &solar_inst[0][0]);
                    if (status != NC_NOERR) {
                        printf("-E- Error reading solar_inst variable.\n");
                        exit(FATAL_ERROR);
                    };

                    for (i = 0; i < nscans; i++) {
                        if (solar_inst[i][0] >= -1 && solar_inst[i][0] <= 1) {
                            float solar_inst_az = RADEG * atan2f(solar_inst[i][0], -solar_inst[i][2]);
                            if (solar_inst_az >= 100.0 && solar_inst_az <= 110.0) {
                                sun_in_sd = TRUE;
                                break;
                            }
                        }
                    }

                    free(solar_inst);
                    free(tmp_arr);
                }
                if ((nc_inq_ncid(geoid, "scan_line_attributes", &grpid)) == NC_NOERR) {
                    size_t start[] = {0}; /* start at first value */
                    size_t cnt[] = {nscans};
                    status = nc_inq_varid(grpid, "scan_quality", &varid);
                    if (status != NC_NOERR) {
                        printf("-E- Error finding scan_quality variable.\n");
                        exit(FATAL_ERROR);
                    };
                    int16_t *scan_quality; //(number_of_scans);
                    if ((scan_quality = malloc(nscans * sizeof (int16_t))) == NULL) {
                        printf("%s -Error: Cannot allocate memory to scan_quality data\n",
                                __FILE__);
                        exit(FATAL_ERROR);
                    }
                    status = nc_get_vara_short(grpid, varid, start, cnt, scan_quality);
                    if (status != NC_NOERR) {
                        printf("-E- Error reading scan_quality variable.\n");
                        exit(FATAL_ERROR);
                    };
                    for (i = 0; i < nscans; i++) {
                        if (scan_quality[i] & SECTOR_ROTATION) {
                            sector_rotation = TRUE;
                            break;  
                        }
                    }

                    free(scan_quality);
                }
            }

            status = nc_close(geoid);
            if (status != NC_NOERR) {
                printf("-E- %s: Error closing %s.\n", argv[0], input.geofile);
                exit(FATAL_ERROR);
            };
        }
    }

    /*									*/
    /* Open input file and get sensor and scan information from handle  */
    /*                                                                  */
    /* done in msl12_input_defaults()
    if (getFormat(&l1file) == -1) {
        printf("-E- %s: Unrecognized input file format.\n",argv[0]);
        exit(FATAL_ERROR);
    }
     */
    if (openl1(&l1file) != 0) {
        printf("-E- %s: Error opening %s for reading.\n", argv[0], l1file.name);
        exit(FATAL_ERROR);
    }

    /*									*/
    /* Allocate memory for L1 scan data			  	        */
    /*									*/
    if (alloc_l1(l1file.npix, input.nbands, NBANDSIR, l1file.n_inprods, &l1rec) == 0) {
        printf("-E- %s: Unable to allocate L1 record.\n", argv[0]);
        exit(FATAL_ERROR);
    }

    spix = 0;
    epix = l1file.npix - 1;
    cpix = spix + (epix - spix + 1) / 2;

    sscan = 0;
    escan = l1file.nscan - 1;
    cscan = sscan + (escan - sscan + 1) / 2;

    l1file.spix = spix;
    l1file.epix = epix;
    l1file.dpix = 1;
    
    if (show_standard_info == TRUE) {
        fprintf(fp, "Sensor=%s\n", sensorName[l1file.sensorID]);
    }

    if (show_sdps_info == TRUE || show_standard_info == TRUE) {
        fprintf(fp, "Orbit_Number=%d\n", l1file.orbit_number);
        if (l1file.format == FMT_VIIRSL1A || l1file.format == FMT_VIIRSL1BNC) {
            fprintf(fp, "Number_of_Scans=%d\n", filledScans);
        } else {
            fprintf(fp, "Number_of_Scans=%d\n", l1file.nscan);
        }            
    }

    if (show_standard_info == TRUE) {
        fprintf(fp, "Number_of_Pixels_per_Scan=%d\n", l1file.npix);
        /*
        printf("Percent_Clouds=%d\n",l1file.percent_cloud);
        printf("Percent_Land=%d\n",l1file.percent_land);
        printf("Percent_Water=%d\n",l1file.percent_water);
         */
    }


    if (num_boxes > 0) {
        box = (box_t **) calloc(num_boxes, sizeof (box_t*));

        if (box == NULL) {
            printf("out of memory\n");
            exit(FATAL_ERROR);
        }
    }

    for (counter = 0; counter < num_boxes; counter++) {
        box[counter] = calloc(NUM_SPACECRAFT_DIRECTION_MODES, sizeof (box_t));
        if (box[counter] == NULL) {
            printf("out of memory\n");
            exit(FATAL_ERROR);
        }
    }



    /***********************************************************************************************
     *    Initialize boxes
     ***********************************************************************************************/

    for (box_index = 0; box_index < num_boxes; box_index++) {
        for (spacecraft_direction_index = 0; spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; spacecraft_direction_index++) {
            box[box_index][spacecraft_direction_index].north_lat = DEFAULT_COORD_VALUE;
            box[box_index][spacecraft_direction_index].south_lat = DEFAULT_COORD_VALUE;
            box[box_index][spacecraft_direction_index].west_lon = DEFAULT_COORD_VALUE;
            box[box_index][spacecraft_direction_index].east_lon = DEFAULT_COORD_VALUE;
            box[box_index][spacecraft_direction_index].daynightflag = DEFAULT_DAYNIGHT_MODE;
        }
    }


    /***************************************************************************
     *    Loop through all scan lines
     ***************************************************************************/

    /*
        for (iscan=sscan; iscan<=escan; iscan++) 
        {  */

    iscan = sscan;
    line_nav_ok = 0;
    start_printed = FALSE;
    while (iscan <= escan) {
        good_spix = -1;
        good_epix = -1;
        good_cpix = -1;
        readl1(&l1file, iscan, &l1rec);
        /*
         *  get info for each line and remember the info if end 
         *  and center geolocation are OK
         */
        if (*l1rec.year > 0) {
            utime = yds2unix(*l1rec.year, *l1rec.day, *l1rec.msec / ((double) 1000.0));
            isodatetime = unix2isodate(utime, 'G');
            unix2ymdhms(utime, &year, &month, &day, &hour, &minute, &second);
            if (start_printed == FALSE) {
                if (show_standard_info == TRUE) {
                    fprintf(fp, "\n");
                    fprintf(fp, "INFO: FIRST SCAN LINE\n");
                    fprintf(fp, "Start_Date=%s\n", isodatetime);
                }
                if (show_sdps_info == TRUE) {
                    fprintf(fp, "Start_Date=%4d-%02d-%02d\n", year, month, day);
                    fprintf(fp, "Start_Time=%02d:%02d:%06.3f\n", hour, minute, second);
                }
                start_printed = TRUE;
            }
        }
        /*
         * establish good start, end, center pixels if any
         */
        for (ip = spix; ip <= epix; ip++) {
            if ((l1rec.flags[ip] & NAVFAIL) == 0) {
                if (good_spix == -1)
                    good_spix = ip;
                good_epix = ip;
                if ((ip <= cpix) || (good_cpix == -1))
                    good_cpix = ip;
                if ((l1rec.flags[ip] & NAVWARN) == 0)
                    navwarn_all_set = 0;
            }
        }
        if (good_spix != -1) {
            line_nav_ok = 1;
            if (first_good_scan == -1) {
                first_good_scan = iscan;
                if (show_standard_info == TRUE) {
                    fprintf(fp, "Upper_Left_Lon=%f\n", l1rec.lon[good_spix]);
                    fprintf(fp, "Upper_Left_Lat=%f\n", l1rec.lat[good_spix]);
                    fprintf(fp, "Upper_Right_Lon=%f\n", l1rec.lon[good_epix]);
                    fprintf(fp, "Upper_Right_Lat=%f\n", l1rec.lat[good_epix]);
                }

            }

            mem_st_lat = l1rec.lat[good_spix];
            mem_st_lon = l1rec.lon[good_spix];
            mem_en_lat = l1rec.lat[good_epix];
            mem_en_lon = l1rec.lon[good_epix];

        }

        if (iscan == escan) {
            if (show_standard_info == TRUE) {
                fprintf(fp, "\n");
                fprintf(fp, "INFO: LAST SCAN LINE\n");
                fprintf(fp, "End_Date=%s\n", isodatetime);
                fprintf(fp, "Lower_Left_Lon=%f\n", mem_st_lon);
                fprintf(fp, "Lower_Left_Lat=%f\n", mem_st_lat);
                fprintf(fp, "Lower_Right_Lon=%f\n", mem_en_lon);
                fprintf(fp, "Lower_Right_Lat=%f\n", mem_en_lat);
            }

            if (show_sdps_info == TRUE) {
                fprintf(fp, "End_Date=%4d-%02d-%02d\n", year, month, day);
                fprintf(fp, "End_Time=%02d:%02d:%06.3f\n", hour, minute, second);
            }
        }


        /***************************************************************************************************
         *    Determine MAX and MIN values for epix_lat and spix_lat based on current and all preceeding scan lines
         ***************************************************************************************************/


        if (line_nav_ok) {
            set_north_south_boundaries(l1rec.lat[good_spix], &spix_northern_lat, &spix_southern_lat);
            set_north_south_boundaries(l1rec.lat[good_epix], &epix_northern_lat, &epix_southern_lat);
        }

        /**************************************************************
         *    Find lat and lon of the center pixel of the middle scan
         **************************************************************/

        if (center_lat == DEFAULT_COORD_VALUE || center_lon == DEFAULT_COORD_VALUE) {
            if (iscan >= cscan && (l1rec.flags[good_cpix] & NAVFAIL) == 0) {
                center_lat = l1rec.lat[good_cpix];
                center_lon = l1rec.lon[good_cpix];
            }
        }


        /***************************************************************************
         *    Determine spacecraft direction based on center pixel of the scan
         ***************************************************************************/

        if ((l1rec.flags[good_cpix] & NAVFAIL) == 0) {
            if (prev_lat_cpix != DEFAULT_COORD_VALUE) {
                if (l1rec.lat[good_cpix] > prev_lat_cpix) {
                    curr_spacecraft_direction = ASCENDING_TRACK;
                } else if (l1rec.lat[good_cpix] < prev_lat_cpix) {

                    curr_spacecraft_direction = DESCENDING_TRACK;
                } else {
                    curr_spacecraft_direction = prev_spacecraft_direction;
                }

                if (initial_spacecraft_direction == DEFAULT_SPACECRAFT_DIRECTION &&
                        curr_spacecraft_direction != DEFAULT_SPACECRAFT_DIRECTION) {
                    initial_spacecraft_direction = curr_spacecraft_direction;
                }

            }
            prev_lat_cpix = l1rec.lat[good_cpix];
            prev_spacecraft_direction = curr_spacecraft_direction;

        }


        /***************************************************************************
         *    Loop through all pixels 
         ***************************************************************************/

        /* for (ip=spix; ip<=epix; ip++) */

        ip = good_spix;
        
        while (ip <= good_epix) {

            if ((l1rec.flags[ip] & NAVFAIL) == 0) {
                /***************************************************************************
                 *    Determine day/night mode of pixel 
                 ***************************************************************************/
                // needed to test for valid solz - invalid ones were incorrectly
                // changing the daynight mode
                if (l1rec.solz[ip] > 0. && l1rec.solz[ip] < 180.){
                    if (l1rec.solz[ip] < 90) {
                        curr_daynight_mode = DAY_MODE;
                    } else {
                        curr_daynight_mode = NIGHT_MODE;
                    }
                }

                master_daynightflag |= curr_daynight_mode;


                /***************************************************************************************************
                 *    Determine MAX and MIN values for both lat and lon based on current and all preceeding scan lines 
                 ***************************************************************************************************/
                if (line_nav_ok) {
                    set_north_south_boundaries(l1rec.lat[ip], &northern_lat, &southern_lat);
                    set_west_east_boundaries(l1rec.lon[ip], &western_lon, &eastern_lon);
                }


                if (num_boxes > 0) {
                    /***************************************************************************
                     *    Determine and set boundaries for regular boxes
                     ***************************************************************************/

                    if (l1rec.lat[ip] >= -90. && l1rec.lat[ip] <= 90. && l1rec.lon[ip] >= -180. && l1rec.lon[ip] <= 180.) {
                        i = 0;
                        box_index = num_boxes;
                        increment = 180. / num_boxes;

                        for (lat_breakpoint = 90. - increment; lat_breakpoint >= -90.; lat_breakpoint -= increment) {
                            if (l1rec.lat[ip] >= lat_breakpoint && l1rec.lat[ip] < (lat_breakpoint + increment)) {
                                box_index = i;
                            }
                            i++;

                        }

                        if (box_index < num_boxes) {
                            set_north_south_boundaries(l1rec.lat[ip],
                                    &box[box_index][curr_spacecraft_direction].north_lat,
                                    &box[box_index][curr_spacecraft_direction].south_lat);
                            set_west_east_boundaries(l1rec.lon[ip],
                                    &box[box_index][curr_spacecraft_direction].west_lon,
                                    &box[box_index][curr_spacecraft_direction].east_lon);
                            box[box_index][curr_spacecraft_direction].daynightflag |= curr_daynight_mode;
                            box[box_index][curr_spacecraft_direction].pixel_count++;
                        }
                    } else {
                        num_bad_lat_lon_pixels++;

                        if (l1rec.lat[ip] < -90. || l1rec.lat[ip] > 90.) {
                            fprintf(stderr, "WARNING_MSG: lat pixel out of range: lat=%f, pixel_index=%d", l1rec.lat[ip], ip);
                        }

                        if (l1rec.lon[ip] < -180. || l1rec.lon[ip] > 180.) {
                            fprintf(stderr, "WARNING_MSG: lon pixel out of range: lon=%f, pixel_index=%d", l1rec.lon[ip], ip);
                        }
                    }
                }
            } else {
                num_nav_fail++; /* count failed pixels  */
            }

            num_pixels++; /* count total pixels */

            if (ip < epix - user_defined_increment) {
                ip += user_defined_increment;
            } else {
                ip++;
            }
        }


        if (iscan < escan - user_defined_increment) {
            iscan += user_defined_increment;
        } else {
            iscan++;
        }
    }

    final_spacecraft_direction = curr_spacecraft_direction;


    if (num_bad_lat_lon_pixels > 0) {
        fprintf(stderr, "WARNING_MSG: found %d lat/lon pixel(s) which were out of range", num_bad_lat_lon_pixels);
    }



    /***************************************************************************************************
     *    Merge initial scans box with appropriate box now that we know the initial_spacecraft_direction 
     ***************************************************************************************************/

    for (box_index = 0; box_index < num_boxes; box_index++) {
        if (box[box_index][DEFAULT_SPACECRAFT_DIRECTION].north_lat != DEFAULT_COORD_VALUE) {
            set_north_south_boundaries(box[box_index][DEFAULT_SPACECRAFT_DIRECTION].north_lat,
                    &box[box_index][initial_spacecraft_direction].north_lat,
                    &box[box_index][initial_spacecraft_direction].south_lat);
            set_north_south_boundaries(box[box_index][DEFAULT_SPACECRAFT_DIRECTION].south_lat,
                    &box[box_index][initial_spacecraft_direction].north_lat,
                    &box[box_index][initial_spacecraft_direction].south_lat);
            set_west_east_boundaries(box[box_index][DEFAULT_SPACECRAFT_DIRECTION].west_lon,
                    &box[box_index][initial_spacecraft_direction].west_lon,
                    &box[box_index][initial_spacecraft_direction].east_lon);
            set_west_east_boundaries(box[box_index][DEFAULT_SPACECRAFT_DIRECTION].east_lon,
                    &box[box_index][initial_spacecraft_direction].west_lon,
                    &box[box_index][initial_spacecraft_direction].east_lon);
            box[box_index][initial_spacecraft_direction].pixel_count += box[box_index][DEFAULT_SPACECRAFT_DIRECTION].pixel_count;
            box[box_index][initial_spacecraft_direction].daynightflag |= box[box_index][DEFAULT_SPACECRAFT_DIRECTION].daynightflag;

            /* re initialize default spacecraft direction box */
            box[box_index][DEFAULT_SPACECRAFT_DIRECTION].north_lat = DEFAULT_COORD_VALUE;
            box[box_index][DEFAULT_SPACECRAFT_DIRECTION].south_lat = DEFAULT_COORD_VALUE;
            box[box_index][DEFAULT_SPACECRAFT_DIRECTION].west_lon = DEFAULT_COORD_VALUE;
            box[box_index][DEFAULT_SPACECRAFT_DIRECTION].east_lon = DEFAULT_COORD_VALUE;
            box[box_index][DEFAULT_SPACECRAFT_DIRECTION].daynightflag = DEFAULT_DAYNIGHT_MODE;
        }
    }

    if (num_boxes != 0) {
        if (num_pixels == num_nav_fail) {
            fprintf(stderr, "ERROR_MSG: all %d pixels have navigation failure\n", num_pixels);
            exitflag = 100;
        }
        else if (northern_lat == DEFAULT_COORD_VALUE ||
                southern_lat == DEFAULT_COORD_VALUE ||
                eastern_lon == DEFAULT_COORD_VALUE ||
                western_lon == DEFAULT_COORD_VALUE) {
            fprintf(stderr, "ERROR_MSG: lat/lon coordinates never set in this program\n");
            exitflag = 101;
        }
        else if (case_n == TRUE && (initial_spacecraft_direction == DEFAULT_SPACECRAFT_DIRECTION ||
                final_spacecraft_direction == DEFAULT_SPACECRAFT_DIRECTION)
                ) {
            fprintf(stderr, "ERROR_MSG: initial/final direction never set in this program\n");
            exitflag = 102;
        }

    }
    if (navwarn_all_set == 1) {
        fprintf(stderr,
                "ERROR_MSG: All pixels in the granule have NAVWARN or NAVFAIL set\n");
        exitflag = 100;
    }

    /***************************************************************************************************
     *    Print summary info 
     ***************************************************************************************************/

    if (show_standard_info == TRUE) {
        fprintf(fp, "\n");
        fprintf(fp, "SUMMARY STATS\n");
        fprintf(fp, "Northernmost_Lat=%f\n", northern_lat);
        fprintf(fp, "Southernmost_Lat=%f\n", southern_lat);
        fprintf(fp, "Easternmost_Lon=%f\n", eastern_lon);
        fprintf(fp, "Westernmost_Lon=%f\n", western_lon);
        fprintf(fp, "Center_Lat=%f\n", center_lat);
        fprintf(fp, "Center_Lon=%f\n", center_lon);
        fprintf(fp, "Start_Node=%s\n", node[initial_spacecraft_direction]);
        fprintf(fp, "End_Node=%s\n", node[final_spacecraft_direction]);
        fprintf(fp, "Daynight=%s\n", daynight_string[master_daynightflag]);
        fprintf(fp, "Moon_in_SV=%1d\n", l1file.sv_with_moon);
        if (input.sensorID == VIIRS)
            fprintf(fp, "Sun_in_SD=%1d\n", sun_in_sd);
            fprintf(fp, "Sector_Rotation=%1d\n", sector_rotation);
    }


    if (show_diagnostic_info == TRUE) {
        fprintf(fp, "Total_Num_Pixels=%d\n", num_pixels);
        fprintf(fp, "Num_Nav_Fail_Pixels=%d\n", num_nav_fail);
        fprintf(fp, "\n");
        fprintf(fp, "End_of_Scan_Northernmost_Lat=%f\n", epix_northern_lat);
        fprintf(fp, "End_of_Scan_Southernmost_Lat=%f\n", epix_southern_lat);
        fprintf(fp, "Start_of_Scan_Northernmost_Lat=%f\n", spix_northern_lat);
        fprintf(fp, "Start_of_Scan_Southernmost_Lat=%f\n", spix_southern_lat);
        fprintf(fp, "\n");
    }



    if (show_sdps_info == TRUE) {
        fprintf(fp, "DayNightFlag=%d\n", master_daynightflag);
        fprintf(fp, "Moon_in_SV=%1d\n", l1file.sv_with_moon);
        if (input.sensorID == VIIRS)
            fprintf(fp, "Sun_in_SD=%1d\n", sun_in_sd);
            fprintf(fp, "Sector_Rotation=%1d\n", sector_rotation);
    }

    if (exitflag == 0 && master_daynightflag == 0) {
        exitflag = 103;
        fprintf(stderr, "ERROR_MSG: master_daynightflag never set in this program\n");
    }


    /***************************************************************************************************
     *    Print box info 
     ***************************************************************************************************/

    for (spacecraft_direction_index = 0; spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; spacecraft_direction_index++) {
        for (box_index = 0; box_index < num_boxes; box_index++) {
            if (box[box_index][spacecraft_direction_index].north_lat != DEFAULT_COORD_VALUE) {
                total_box_count++;

                tmp_north_lat = box[box_index][spacecraft_direction_index].north_lat;
                tmp_south_lat = box[box_index][spacecraft_direction_index].south_lat;
                tmp_west_lon = box[box_index][spacecraft_direction_index].west_lon;
                tmp_east_lon = box[box_index][spacecraft_direction_index].east_lon;

                get_box_square_km((double) tmp_north_lat,
                        (double) tmp_south_lat,
                        (double) tmp_west_lon,
                        (double) tmp_east_lon,
                        &box_square_km);

                total_box_square_km += box_square_km;

                if (show_diagnostic_info == TRUE) {
                    fprintf(fp, "Box_Summary\n");
                    fprintf(fp, "Box_Northern_boundary=%f\n", box[box_index][spacecraft_direction_index].north_lat);
                    fprintf(fp, "Box_Southern_boundary=%f\n", box[box_index][spacecraft_direction_index].south_lat);
                    fprintf(fp, "Box_Western_boundary=%f\n", box[box_index][spacecraft_direction_index].west_lon);
                    fprintf(fp, "Box_Eastern_boundary=%f\n", box[box_index][spacecraft_direction_index].east_lon);
                    fprintf(fp, "Box_Daynightflag=%s\n", daynight_string[box[box_index][spacecraft_direction_index].daynightflag]);
                    fprintf(fp, "Box_Pixel_count=%d\n", box[box_index][spacecraft_direction_index].pixel_count);
                    fprintf(fp, "Box_size_in_square_million_km=%.2f\n", box_square_km / 1000000.);
                    fprintf(fp, "\n");
                }
                if (show_sdps_info == TRUE) {
                    fprintf(fp, "GeoBox_%d=%f,%f,%f,%f,%d\n",
                            total_box_count,
                            box[box_index][spacecraft_direction_index].north_lat,
                            box[box_index][spacecraft_direction_index].south_lat,
                            box[box_index][spacecraft_direction_index].west_lon,
                            box[box_index][spacecraft_direction_index].east_lon,
                            box[box_index][spacecraft_direction_index].daynightflag);
                }


            }
        }
    }






    if (show_diagnostic_info == TRUE) {
        fprintf(fp, "Total_box_count=%d\n", total_box_count);
        fprintf(fp, "Total_box_size_square_million_km=%.2f\n", total_box_square_km / 1000000.);
    }

    if (exitflag == 0 && case_n == TRUE && total_box_count == 0 && num_boxes != 0) {
        exitflag = 104;
        fprintf(stderr, "ERROR_MSG: total_box_count = 0\n");
    }


    if (case_v == TRUE) {

        for (iscan = sscan; iscan <= escan; iscan++) {
            readl1(&l1file, iscan, &l1rec);

            for (ip = spix; ip <= epix; ip++) {
                if ((l1rec.flags[ip] & NAVFAIL) == 0) {
                    match = FALSE;

                    for (box_index = 0; box_index < num_boxes; box_index++) {
                        for (spacecraft_direction_index = 0; spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; spacecraft_direction_index++) {
                            in_box = check_if_in_box(l1rec.lat[ip],
                                    l1rec.lon[ip],
                                    box[box_index][spacecraft_direction_index].north_lat,
                                    box[box_index][spacecraft_direction_index].south_lat,
                                    box[box_index][spacecraft_direction_index].west_lon,
                                    box[box_index][spacecraft_direction_index].east_lon);
                            if (in_box == TRUE) {
                                match = TRUE;
                            }
                        }
                    }

                    if (match == TRUE) {
                        num_match++;
                    } else {
                        num_no_match++;
                    }

                }
            }
        }

        if ((num_match + num_no_match) > 0) {
            percent_match = 100. * num_match / (num_match + num_no_match);
        } else {
            percent_match = 0;
        }

        fprintf(fp, "percent_match=%f, num_match=%d, num_no_match=%d\n", percent_match, num_match, num_no_match);
    }


    /* free the dynamically allocated memory */

    if (num_boxes > 0) {
        for (counter = 0; counter < num_boxes; counter++) {
            free(box[counter]);
        }
        /* and finally free ptr to ptr....the space to store ptr to box_t*/
        free(box);
    }
    if (case_o == TRUE)
        fclose(fp);

    exit(exitflag);
}

/* take first arg (lat) and adjust northern_boundary or southern_boundary if needed  */

void
set_north_south_boundaries(float32 lat, float32 *northern_boundary, float32 *southern_boundary) {
    if (lat > 90.) {
        lat = 90.;
    }

    if (lat < -90.) {
        lat = -90.;
    }

    if (*northern_boundary != DEFAULT_COORD_VALUE) {
        *northern_boundary = MAX(lat, *northern_boundary);
    } else {
        *northern_boundary = lat;
    }

    if (*southern_boundary != DEFAULT_COORD_VALUE) {
        *southern_boundary = MIN(lat, *southern_boundary);
    } else {
        *southern_boundary = lat;
    }

}

double
get_lon_distance(double lon1, double lon2) {
    double results;
    double distance_1;
    double distance_2;

    distance_1 = lon1 - lon2;
    distance_2 = lon2 - lon1;
    distance_1 = MAX(distance_1, distance_2);

    distance_2 = 360. - distance_1;

    results = MIN(distance_1, distance_2);

    return results;
}

int
check_if_in_box(float32 lat, float32 lon, float32 northern_boundary, float32 southern_boundary,
        float32 western_boundary, float32 eastern_boundary) {
    int results;
    int east_west_results;
    int north_south_results;

    east_west_results = check_if_in_west_east_boundaries(lon, western_boundary, eastern_boundary);

    if (lat >= southern_boundary && lat <= northern_boundary) {
        north_south_results = TRUE;
    } else {
        north_south_results = FALSE;
    }

    if (north_south_results == TRUE && east_west_results == TRUE) {
        results = TRUE;
    } else {
        results = FALSE;
    }

    return results;
}

int
check_if_in_west_east_boundaries(float32 lon, float32 western_boundary, float32 eastern_boundary) {

    int results = FALSE;

    if (eastern_boundary >= western_boundary) {
        /* no date line crossing */

        if (lon >= western_boundary &&
                lon <= eastern_boundary) {
            results = TRUE;
        }
    } else {
        /* date line crossing */

        if (lon >= western_boundary ||
                lon <= eastern_boundary) {
            results = TRUE;
        }
    }

    return results;
}

void
get_box_square_km(double northern_boundary,
        double southern_boundary,
        double western_boundary,
        double eastern_boundary,
        double * box_square_km) {

    double northern_box_square_km;
    double southern_box_square_km;


    if (northern_boundary > 0. && southern_boundary < .0) {
        /*************************************************************************************
         *    Divide into two boxes since equator is crossed
         *************************************************************************************/

        get_box_square_km_sublevel(northern_boundary,
                (double) 0.,
                western_boundary,
                eastern_boundary,
                &northern_box_square_km);

        get_box_square_km_sublevel((double) 0.,
                southern_boundary,
                western_boundary,
                eastern_boundary,
                &southern_box_square_km);

        *box_square_km = northern_box_square_km + southern_box_square_km;
    } else {
        get_box_square_km_sublevel(northern_boundary,
                southern_boundary,
                western_boundary,
                eastern_boundary,
                box_square_km);
    }

}

void
get_box_square_km_sublevel(double northern_boundary,
        double southern_boundary,
        double western_boundary,
        double eastern_boundary,
        double * box_square_km) {
    double lat_width_degrees;
    double lon_width_degrees;
    double north_lon_radius;
    double south_lon_radius;
    double north_lon_width_km;
    double south_lon_width_km;
    double avg_lon_width_km;
    double lat_width_km;

    if (eastern_boundary >= western_boundary) {
        /* no date line crossing */
        lon_width_degrees = eastern_boundary - western_boundary;
    } else {
        /* date line crossing */
        lon_width_degrees = 360. + eastern_boundary - western_boundary;
    }

    north_lon_radius = EARTH_RADIUS_EQUATORIAL * fabs(cos(northern_boundary * (PI / 180.)));
    north_lon_width_km = (PI / 180.) * lon_width_degrees * north_lon_radius;

    south_lon_radius = EARTH_RADIUS_EQUATORIAL * fabs(cos(southern_boundary * (PI / 180.)));
    south_lon_width_km = (PI / 180.) * lon_width_degrees * south_lon_radius;

    avg_lon_width_km = (north_lon_width_km + south_lon_width_km) / 2.;

    lat_width_degrees = northern_boundary - southern_boundary;
    lat_width_km = (PI / 180.) * lat_width_degrees * EARTH_RADIUS_EQUATORIAL;

    *box_square_km = lat_width_km * avg_lon_width_km;

}

void
set_west_east_boundaries(float32 lon, float32 * western_boundary, float32 * eastern_boundary) {
    float32 boundary_width;
    float32 width_to_west;
    float32 width_to_east;

    if (*western_boundary != -180. || *eastern_boundary != 180.) {
        if (*eastern_boundary != DEFAULT_COORD_VALUE && *western_boundary != DEFAULT_COORD_VALUE) {
            if (check_if_in_west_east_boundaries(lon, *western_boundary, *eastern_boundary) != TRUE) {
                /************************************************************************
                 *    Determine longitude width to east of current boundary
                 ************************************************************************/

                if (lon > *eastern_boundary) {
                    /* dateline not crossed */
                    width_to_east = lon - *eastern_boundary;
                } else {
                    /* dateline crossed */
                    width_to_east = 360. + lon - *eastern_boundary;
                }

                /************************************************************************
                 *    Determine longitude width to west of current boundary
                 ************************************************************************/

                if (*western_boundary > lon) {
                    /* dateline not crossed */
                    width_to_west = *western_boundary - lon;
                } else {
                    /* dateline crossed */
                    width_to_west = *western_boundary + 360. - lon;
                }

                /************************************************************************
                 *    Set closest west-east boundary
                 ************************************************************************/

                if (fabs(width_to_west) <= fabs(width_to_east)) {
                    *western_boundary = lon;
                } else {
                    *eastern_boundary = lon;
                }

                /************************************************************************
                 *    Determine longitude width between western_boundary and eastern_boundary
                 ************************************************************************/

                if (*eastern_boundary >= *western_boundary) {
                    /* no date line crossing */
                    boundary_width = *eastern_boundary - *western_boundary;
                } else {
                    /* date line crossing */
                    boundary_width = 360. + *eastern_boundary - *western_boundary;
                }

                /************************************************************************
                 *    if west-to-east span > 355 then just set span to -180 to 180
                 ************************************************************************/

                if (boundary_width > 355.) {
                    *western_boundary = -180.;
                    *eastern_boundary = 180.;
                }
            }
        } else {
            *eastern_boundary = lon;
            *western_boundary = lon;
        }
    }
}
