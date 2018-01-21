#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <list>
#include "hdf.h"
#include "mfhdf.h"
#include "l1info_aquarius.h"


#define NUM_SPACECRAFT_DIRECTION_MODES 3
#define DEFAULT_SPACECRAFT_DIRECTION 0   /* set to an invalid value */
#define ASCENDING 1
#define DESCENDING 2

#define NORTH 1
#define SOUTH 0
//#define TRUE 1
#define FALSE 0

#define TMP_FILENAME_MAX 255

#define DEFAULT_COORD_VALUE -999.   /* set to an invalid value */

#define DEFAULT_DAYNIGHT_MODE 0
#define DAY_MODE (1<<0)
#define NIGHT_MODE (1<<1)
#define FATAL_ERROR 110

/* equatorial radius, didn't bother to use polar value also since precision not that important for current use */
#define EARTH_RADIUS_EQUATORIAL 6378   

#define VERSION "1.03"

//    Modification history:
//
//     Programmer     Organization      Date      Description of change
//   --------------   ------------    --------    ---------------------
//   Joel Gales       FutureTech      11/12/09    Add support for InT L1A
//   Joel Gales       FutureTech      04/12/10    Add support for V2 geo
//   Joel Gales       FutureTech      12/02/10    Add support for V3 geo
//   Joel Gales       FutureTech      08/10/11    Update to V 1.00
//   Joel Gales       FutureTech      07/05/12    Add zeroed rpy_adj
//                                                parameter to geolocation
//   Joel Gales       FutureTech 1.03 12/06/12    Add gsl_interp_init
//                                                for compatibility with
//                                                Turn off gsl error handler
//                                                for attitude interpolation

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
int check_if_in_west_east_boundaries(float32, float32, float32);
void set_west_east_boundaries(float32, float32 *, float32 *);
double get_lon_distance(double,double);

double gpstai2unix( double gpstai);
double gpstai2utc2000( double gpstai);

using namespace std;

int main (int argc, char* argv[])
{
  int num_boxes = 20;
  const char *option_string = "n:s";
  int options = 0;

  cout << "l1info_aquarius " << VERSION << " (" 
       <<  __DATE__ << " " << __TIME__ << ")" << endl;

  if ( argc == 1) {
    cout << endl << "l1info_aquarius [-n nboxes] [-s] l1a_filename"
	 << endl;
    cout << "     -n:  Number of geoboxes [Default: 20]" << endl;
    cout << "     -s:  Display info for database" << endl;
    return 0;
  }

  bool case_s = false;
  while ((options = getopt(argc, argv, option_string)) != -1)
    {
      switch (options)
        {
	case 'n':
	  num_boxes = atoi(optarg);
	  break;
	case 's':
	  case_s = true;
	  break;
	default:
	  break;
        }
    }                          

  static Hdf::hdf5_Aquarius l1afile;
  l1afile.openl1(argv[optind+0], H5F_ACC_RDONLY);

  double granuleStart, granuleStop;
  l1afile.getGranuleTimes( &granuleStart, &granuleStop);

  time_t rawtime;
  tm * ptm;
  char buffer [80];

  cout << endl;
  cout << "Number_of_Blocks=" << l1afile.nBlks() << endl;

  rawtime = (unsigned long) gpstai2unix( granuleStart);
  ptm = gmtime ( &rawtime );
  strftime (buffer, 80, "%Y-%m-%d", ptm);
  cout << "Start_Date=" << buffer << endl;
  strftime (buffer, 80, "%H:%M:%S", ptm);
  cout << "Start_Time=" << buffer << endl;

  rawtime = (unsigned long) gpstai2unix( granuleStop);
  ptm = gmtime ( &rawtime );
  strftime (buffer, 80, "%Y-%m-%d", ptm);
  cout << "End_Date=" << buffer << endl;
  strftime (buffer, 80, "%H:%M:%S", ptm);
  cout << "End_Time=" << buffer << endl;

  cout << endl;
  cout << endl;

  // Check if InT L1A (no eph/att)
  uint32_t numOrbVec;
  if ( l1afile.readl1_eph( &numOrbVec) == -1) return 0;

  if ( case_s) return 0;


  double nodeCrossingTime;
  float nodeLongitude;
  l1afile.getNodeInfo( &nodeCrossingTime, &nodeLongitude);


  float *cellon, *cellat, *celtht, *celphi, *suntht, *sunphi, 
    *glxlon, *glxlat;
  double *zang;
  cellon = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  cellat = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  celtht = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  celphi = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  suntht = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  sunphi = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  glxlon = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  glxlat = new float[l1afile.nBlks()*NUMBER_OF_BEAMS];
  zang = new double[l1afile.nBlks()];

  getGeoNavSun( &l1afile,
		cellon, cellat, celtht, celphi,
		suntht, sunphi, glxlon, glxlat,
		zang);

  int counter, box_index;
  int   spacecraft_direction_index;
  int   curr_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
  int   prev_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
  int   initial_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;
  int   final_spacecraft_direction = DEFAULT_SPACECRAFT_DIRECTION;

  float32 prev_lat_cpix = DEFAULT_COORD_VALUE;
  float32 lat_breakpoint;

  box_t **box = NULL;
  box = (box_t **) calloc(num_boxes,sizeof(box_t*));
  if(box == NULL)
    {
      printf("out of memory\n");
      exit(FATAL_ERROR);
    }

  for(counter=0;counter<num_boxes;counter++) {
    box[counter] = (box_t*) calloc(NUM_SPACECRAFT_DIRECTION_MODES,sizeof(box_t));
        if(box[counter] == NULL) {
	  printf("out of memory\n");
	  exit(FATAL_ERROR);
        }
  }

  for (box_index = 0; box_index < num_boxes; box_index++) {
    for (spacecraft_direction_index = 0; 
	 spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; 
	 spacecraft_direction_index++) {
      box[box_index][spacecraft_direction_index].north_lat = DEFAULT_COORD_VALUE;
      box[box_index][spacecraft_direction_index].south_lat = DEFAULT_COORD_VALUE;
      box[box_index][spacecraft_direction_index].west_lon = DEFAULT_COORD_VALUE;
      box[box_index][spacecraft_direction_index].east_lon = DEFAULT_COORD_VALUE;
    }
  }


  // Main Loop
  for (uint32_t iblk=0; iblk<l1afile.nBlks(); iblk++) {

    double blkSec;

    if ((iblk % 500) == 0) cout << "iblk: " << iblk << endl;

    l1afile.readl1_radiometer(iblk, 0, 0, &blkSec, NULL, NULL, NULL, NULL);

    if ( blkSec == -1) continue;

    /***************************************************************************
     *    Determine spacecraft direction based on center pixel of the scan
     ***************************************************************************/
    float lat_cpix = cellat[1+NUMBER_OF_BEAMS*iblk];
    if (prev_lat_cpix != DEFAULT_COORD_VALUE)
      {
	if (lat_cpix > prev_lat_cpix)
	  {
	    curr_spacecraft_direction = ASCENDING;
	  }
	else if (lat_cpix < prev_lat_cpix)
	  {
	    curr_spacecraft_direction = DESCENDING;
	  }
	else
	  {
	    curr_spacecraft_direction = prev_spacecraft_direction;
	  }

	if (initial_spacecraft_direction == DEFAULT_SPACECRAFT_DIRECTION &&
	    curr_spacecraft_direction != DEFAULT_SPACECRAFT_DIRECTION)
	  {
	    initial_spacecraft_direction = curr_spacecraft_direction;
	  }
      }
    prev_lat_cpix  = lat_cpix;
    prev_spacecraft_direction = curr_spacecraft_direction;
    
    float clat, clon;
    // Beam (Horn) Loop
    for (int32_t ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {

      clat = cellat[ibeam+NUMBER_OF_BEAMS*iblk];
      clon = cellon[ibeam+NUMBER_OF_BEAMS*iblk];

      if (clat >= -90. && clat <= 90. && clon >= -180. && 
	  clon <= 180.) {
	int i = 0; 
	box_index = num_boxes; 
	float increment = 180. / num_boxes;
    
	for (lat_breakpoint = 90. - increment; lat_breakpoint >= -90.; lat_breakpoint -= increment)
	  { 
	    if (clat >= lat_breakpoint && clat < (lat_breakpoint + increment) )
	      {
		box_index = i;
	      }
	    i++;
	  }
 
	if (box_index < num_boxes) {
	  set_north_south_boundaries(clat,
				     &box[box_index][curr_spacecraft_direction].north_lat,
				     &box[box_index][curr_spacecraft_direction].south_lat);
	  set_west_east_boundaries(clon,
				   &box[box_index][curr_spacecraft_direction].west_lon,
				   &box[box_index][curr_spacecraft_direction].east_lon);
	  box[box_index][curr_spacecraft_direction].pixel_count++;
	}
      } // ibeam loop
      
    } // iblk (End of Main Loop)

    final_spacecraft_direction = curr_spacecraft_direction;

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
	box[box_index][initial_spacecraft_direction].pixel_count += 
	  box[box_index][DEFAULT_SPACECRAFT_DIRECTION].pixel_count;

	/* re initialize default spacecraft direction box */
	box[box_index][DEFAULT_SPACECRAFT_DIRECTION].north_lat = DEFAULT_COORD_VALUE;
	box[box_index][DEFAULT_SPACECRAFT_DIRECTION].south_lat = DEFAULT_COORD_VALUE;
	box[box_index][DEFAULT_SPACECRAFT_DIRECTION].west_lon = DEFAULT_COORD_VALUE;
	box[box_index][DEFAULT_SPACECRAFT_DIRECTION].east_lon = DEFAULT_COORD_VALUE;
      }
    }
  }

  /***************************************************************************************************
   *    Print box info 
   ***************************************************************************************************/
  for (spacecraft_direction_index = 0; spacecraft_direction_index < NUM_SPACECRAFT_DIRECTION_MODES; spacecraft_direction_index++)
    {
      for (box_index = 0; box_index < num_boxes; box_index++)
        {
	  if (box[box_index][spacecraft_direction_index].north_lat != DEFAULT_COORD_VALUE)
            {
	      printf("GeoBox_%d=%f,%f,%f,%f,%d\n",
		     box_index,
		     box[box_index][spacecraft_direction_index].north_lat,
		     box[box_index][spacecraft_direction_index].south_lat,
		     box[box_index][spacecraft_direction_index].west_lon,
		     box[box_index][spacecraft_direction_index].east_lon,
		     curr_spacecraft_direction);
	    }
	}
    }


  l1afile.closel1();

  delete[] cellon;
  delete[] cellat;
  delete[] celtht;
  delete[] celphi;
  delete[] suntht;
  delete[] sunphi;
  delete[] glxlon;
  delete[] glxlat;
  delete[] zang;

  for(counter = 0; counter < num_boxes; counter++) {
    free(box[counter]);
  }
  /* and finally free ptr to ptr....the space to store ptr to box_t*/
  free(box);

  return 0;
}


int getGeoNavSun( Hdf::hdf5_Aquarius *l1afile, 
		  float *cellon, float *cellat, float *celtht, float *celphi,
		  float *suntht, float *sunphi, float *glxlon, float *glxlat,
		  double *zang)
{
  // Read Orbital Data
  uint32_t numOrbVec;
  l1afile->readl1_eph( &numOrbVec);
  double *orbTime, *orbPos, *orbVel, *dum;
  orbTime = new double[numOrbVec];
  orbPos = new double[numOrbVec*3];
  orbVel = new double[numOrbVec*3];
  dum = new double[numOrbVec];
  l1afile->readl1_eph( orbTime, orbPos, orbVel);

  double *Time, *Pos, *Vel;
  Time = new double[l1afile->nBlks()];
  Pos = new double[l1afile->nBlks()*3];
  Vel = new double[l1afile->nBlks()*3];

  double blkSec;

  // Cubic Interpolation Using orbPos and orbVel
  gsl_matrix *time_mat = gsl_matrix_alloc(4,4);
  gsl_vector *b = gsl_vector_alloc(4);

  gsl_permutation *p = gsl_permutation_alloc(4);
  int s;
  double C, D;

  for (size_t i=0; i<l1afile->nBlks(); i++) {
    l1afile->readl1_radiometer(i, 0, 0, &blkSec, NULL, NULL, NULL, NULL);
    Time[i] = blkSec;
    if ( blkSec == -1) continue;
    for (size_t j=0; j<numOrbVec-1; j++) {
      if ( Time[i] < orbTime[0]) {
	Time[i] = -1;
	break;
      }
      if ( orbTime[j] <= Time[i] && orbTime[j+1] > Time[i]) {
	double time1_1 = orbTime[j] - blkSec;
	double time2_1 = orbTime[j+1] - blkSec;
	double time1_2 = time1_1 * time1_1;
	double time2_2 = time2_1 * time2_1;
	double time1_3 = time1_2 * time1_1;
	double time2_3 = time2_2 * time2_1;
	
	gsl_matrix_set( time_mat, 0, 0, time1_3);
	gsl_matrix_set( time_mat, 0, 1, time1_2);
	gsl_matrix_set( time_mat, 0, 2, time1_1);
	gsl_matrix_set( time_mat, 0, 3, 1.0);

	gsl_matrix_set( time_mat, 1, 0, time2_3);
	gsl_matrix_set( time_mat, 1, 1, time2_2);
	gsl_matrix_set( time_mat, 1, 2, time2_1);
	gsl_matrix_set( time_mat, 1, 3, 1.0);

	gsl_matrix_set( time_mat, 2, 0, 3*time1_2);
	gsl_matrix_set( time_mat, 2, 1, 2*time1_1);
	gsl_matrix_set( time_mat, 2, 2, 1.0);
	gsl_matrix_set( time_mat, 2, 3, 0.0);

	gsl_matrix_set( time_mat, 3, 0, 3*time2_2);
	gsl_matrix_set( time_mat, 3, 1, 2*time2_1);
	gsl_matrix_set( time_mat, 3, 2, 1.0);
	gsl_matrix_set( time_mat, 3, 3, 0.0);

	gsl_linalg_LU_decomp( time_mat, p, &s);


	// X, Vx
	gsl_vector_set( b, 0, orbPos[3*j]);
	gsl_vector_set( b, 1, orbPos[3*(j+1)]);
	gsl_vector_set( b, 2, orbVel[3*j]);
	gsl_vector_set( b, 3, orbVel[3*(j+1)]);
	gsl_linalg_LU_svx( time_mat, p, b);

	C = gsl_vector_get( b, 2);
	D = gsl_vector_get( b, 3);

	Pos[3*i] = D;
	Vel[3*i] = C;

	// Y, Vy
	gsl_vector_set( b, 0, orbPos[3*j+1]);
	gsl_vector_set( b, 1, orbPos[3*(j+1)+1]);
	gsl_vector_set( b, 2, orbVel[3*j+1]);
	gsl_vector_set( b, 3, orbVel[3*(j+1)+1]);
	gsl_linalg_LU_svx( time_mat, p, b);

	C = gsl_vector_get( b, 2);
	D = gsl_vector_get( b, 3);

	Pos[3*i+1] = D;
	Vel[3*i+1] = C;

	// Z, Vz
	gsl_vector_set( b, 0, orbPos[3*j+2]);
	gsl_vector_set( b, 1, orbPos[3*(j+1)+2]);
	gsl_vector_set( b, 2, orbVel[3*j+2]);
	gsl_vector_set( b, 3, orbVel[3*(j+1)+2]);
	gsl_linalg_LU_svx( time_mat, p, b);

	C = gsl_vector_get( b, 2);
	D = gsl_vector_get( b, 3);

	Pos[3*i+2] = D;
	Vel[3*i+2] = C;
      }
    }
  }




  // Read Attitude Data
  uint32_t numAttSamp;
  l1afile->readl1_att( &numAttSamp);
  double *attTime, *attAng;
  attTime = new double[numAttSamp];
  attAng = new double[numAttSamp*3];
  dum = new double[numAttSamp];
  l1afile->readl1_att( attTime, attAng, NULL);

  double *rpy;
  rpy = new double[l1afile->nBlks()*3];


  // Set up interpolation
  gsl_interp_accel *acc;
  gsl_interp *interp;
  acc = gsl_interp_accel_alloc ();
  interp = gsl_interp_alloc (gsl_interp_linear, numAttSamp);

  // Turn off error handler to avoid potential aborts at last attTime
  gsl_error_handler_t *old_gsl_handler = gsl_set_error_handler_off();

  // Roll
  for (size_t i=0; i<numAttSamp; i++) dum[i] = attAng[3*i];
  gsl_interp_init (interp, attTime, dum, numAttSamp);
  for (size_t i=0; i<l1afile->nBlks(); i++)
    rpy[3*i] = gsl_interp_eval (interp, attTime, dum, Time[i], acc);

  // Pitch
  for (size_t i=0; i<numAttSamp; i++) dum[i] = attAng[3*i+1];
  gsl_interp_init (interp, attTime, dum, numAttSamp);
  for (size_t i=0; i<l1afile->nBlks(); i++)
    rpy[3*i+1] = gsl_interp_eval (interp, attTime, dum, Time[i], acc);

  // Yaw
  for (size_t i=0; i<numAttSamp; i++) dum[i] = attAng[3*i+2];
  gsl_interp_init (interp, attTime, dum, numAttSamp);
  for (size_t i=0; i<l1afile->nBlks(); i++)
    rpy[3*i+2] = gsl_interp_eval (interp, attTime, dum, Time[i], acc);

  gsl_set_error_handler( old_gsl_handler);

  gsl_interp_accel_free (acc);
  gsl_interp_free (interp);


  cout << "Compute lon/lat" << endl;
  double time_utc_2000;
  float dumfoot[4*3], dumarr[3];
  double dbldumarr[3], dbldumarr2[3], dumsight[3*3];

  double rpy_adj[3]={0.0,0.0,0.0};

  for (size_t i=0; i<l1afile->nBlks(); i++) {
    if ( Time[i] == -1) continue;
    //    if ((i % 500) == 0) cout << i << endl;
    time_utc_2000 = gpstai2utc2000( Time[i]);

    geolocation_( rpy_adj, &time_utc_2000, &Pos[3*i], &Vel[3*i], &rpy[3*i], 
		  &cellat[NUMBER_OF_BEAMS*i], &cellon[NUMBER_OF_BEAMS*i], 
		  &celtht[NUMBER_OF_BEAMS*i], &celphi[NUMBER_OF_BEAMS*i], 
		  &suntht[NUMBER_OF_BEAMS*i], &sunphi[NUMBER_OF_BEAMS*i], 
		  &dumarr[0], &dumarr[0], 
		  &glxlat[NUMBER_OF_BEAMS*i], &glxlon[NUMBER_OF_BEAMS*i],
		  &zang[i], 
		  &dbldumarr[0], &dbldumarr[0], &dbldumarr[0], &dbldumarr[0], 
		  &dumfoot[0], 
		  &dumfoot[0],
		  &dbldumarr2[0], &dbldumarr2[0], 
		  &dbldumarr2[0], &dbldumarr2[0],
		  &dumsight[0]);
  }

  for (size_t i=0; i<NUMBER_OF_BEAMS*l1afile->nBlks(); i++) {
    cellon[i] = 
      cellon[i]*(cellon[i] <= 180.0) + (cellon[i]-360.)*(cellon[i] > 180.0);
  }


  delete[] orbTime;
  delete[] orbPos;
  delete[] orbVel;
  delete[] attTime;
  delete[] attAng;
  delete[] dum;
  delete[] Pos;
  delete[] Vel;
  delete[] Time;
  delete[] rpy;

  return 0;
}


/* take first arg (lat) and adjust northern_boundary or southern_boundary if needed  */

void
set_north_south_boundaries(float32 lat, float32 *northern_boundary, float32 *southern_boundary)
{
    if (lat > 90.)
    {
        lat = 90.;
    }
    
    if (lat < -90.)
    {
        lat = -90.;
    }

    if (*northern_boundary != DEFAULT_COORD_VALUE)
    {
        *northern_boundary = MAX(lat,*northern_boundary);
    }
    else
    {
        *northern_boundary = lat;
    }

    if (*southern_boundary != DEFAULT_COORD_VALUE)
    {
        *southern_boundary = MIN(lat,*southern_boundary);
    }
    else
    {
        *southern_boundary = lat;
    }

}

double 
get_lon_distance(double lon1,double lon2)
{
    double results;
    double distance_1;
    double distance_2;

    distance_1 = lon1 - lon2;
    distance_2 = lon2 - lon1;
    distance_1 = MAX(distance_1,distance_2);

    distance_2 = 360. - distance_1;

    results = MIN(distance_1,distance_2);

    return results;
}


int
check_if_in_west_east_boundaries(float32 lon, float32 western_boundary, float32 eastern_boundary)
{

    int results = FALSE;

    if (eastern_boundary >= western_boundary)
    {
        /* no date line crossing */

        if (lon >= western_boundary && 
            lon <= eastern_boundary) 
        {
            results = TRUE;
        }
    }
    else
    {
        /* date line crossing */

        if (lon >= western_boundary ||
            lon <= eastern_boundary)
        {
            results = TRUE;
        }
    }

    return results;
}


void
set_west_east_boundaries(float32 lon, float32 * western_boundary, float32 * eastern_boundary)
{
    float32 boundary_width;
    float32 width_to_west;
    float32 width_to_east;

    if (*western_boundary != -180. || *eastern_boundary != 180.) 
    {
        if (*eastern_boundary != DEFAULT_COORD_VALUE && *western_boundary != DEFAULT_COORD_VALUE)
        {
            if (check_if_in_west_east_boundaries(lon,*western_boundary,*eastern_boundary) != TRUE)
            {
                /************************************************************************
                *    Determine longitude width to east of current boundary
                ************************************************************************/
     
                if (lon > *eastern_boundary)
                {
                    /* dateline not crossed */
                    width_to_east = lon - *eastern_boundary;
                }
                else
                {
                    /* dateline crossed */
                    width_to_east = 360. + lon - *eastern_boundary;
                }
    
                /************************************************************************
                *    Determine longitude width to west of current boundary
                ************************************************************************/
     
                if (*western_boundary > lon)
                {
                    /* dateline not crossed */
                    width_to_west = *western_boundary - lon;
                }
                else
                {
                    /* dateline crossed */
                    width_to_west = *western_boundary + 360. - lon;
                }
    
                /************************************************************************
                *    Set closest west-east boundary
                ************************************************************************/
     
                if (fabs(width_to_west) <= fabs(width_to_east))
                {
                    *western_boundary = lon;
                }
                else
                {
                    *eastern_boundary = lon;
                }

                /************************************************************************
                *    Determine longitude width between western_boundary and eastern_boundary
                ************************************************************************/
     
                if (*eastern_boundary >= *western_boundary)
                {
                    /* no date line crossing */
                    boundary_width = *eastern_boundary - *western_boundary;
                }
                else
                {
                    /* date line crossing */
                    boundary_width = 360. + *eastern_boundary - *western_boundary;
                }

                /************************************************************************
                *    if west-to-east span > 355 then just set span to -180 to 180
                ************************************************************************/
     
                if (boundary_width > 355.)
                {
                    *western_boundary = -180.;
                    *eastern_boundary = 180.;
                }
            }
        }
        else
        {
            *eastern_boundary = lon;
            *western_boundary = lon;
        }
    }
}
