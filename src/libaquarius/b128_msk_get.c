#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXPIX      5500
#define NITEMPTR 1024 * 64
#define BOX_MAX 100
#define MAX_PIX_PER_SCENE MAXPIX*10000

static FILE *stream[2] = { NULL, NULL };
static short ptr_arr[2][ NITEMPTR ];
static int access_cnt[2]   = {  0,  0 };
static int last_box_num[2] = { -1, -1 };
static int boxes_active[2] = {  0,  0 };

/* this is an array of structures to hold info on each 1 degree box
   recently needed */
static struct box_info_str
  {
  int box_num;  /* box number stored in this slot */
  int last_access_cnt;  /* last access count for this boxes data */
  unsigned short *box_ptr; /* pointer to the box's data */
  } box_info[2][BOX_MAX];


int endianess( void )
{
    int x = 1;
    return ( *(char *)&x == 1 );
}


int swapc_bytes(in, nbyte, ntime)
char *in;
int  nbyte, ntime;
{
  char *tmpbuf, *ptr;
  int  i, j, k;

  tmpbuf = (char *) malloc(nbyte+1);

  for (j=0; j<ntime; j++)
  {
    ptr = in + j * nbyte;
    memcpy(tmpbuf, ptr, nbyte);

    for (i=0, k=nbyte-1; i<nbyte; i++, k--)
      tmpbuf[i] = ptr[k];

    memcpy(ptr, tmpbuf, nbyte);
  }
  free(tmpbuf);

  return 0;
}


int b128_msk_init( char *landfile, int msknum )
/*******************************************************************

   b128_msk_init

   purpose: initialize the 218 x 128 per degree resolution processing
        to get the a mask

   Returns type: int - 0 is good, else an open or read error

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            landfile         I      name of mask file to use
      int               msknum           I      number of the mask (0 or 1)

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       8-Oct-1997      Original development
      W. Robinson       19-Dec-1997     adapt to use 2 mask files

*******************************************************************/
  {
  short rec_arr[1024];
 /*
  *  make sure valid mask # is input
  */
  if( msknum < 0 || msknum > 1 )
    {
    fprintf(stderr, "b128_msk_init: Fatal error. only mask 0 or 1 can be used\n" );
    fprintf(stderr, "               current mask: %d\n", msknum );
    fprintf(stderr, "Exiting\n" );
    return -1;
    }
 /*
  *  open the mask file
  */
  if( ( stream[msknum] = fopen( landfile, "rb" ) ) == NULL )
    {
    fprintf(stderr, "b128_msk_init: Fatal error. failed to open mask file # %d:\n",
         msknum );
    fprintf(stderr, "    '%s'\n", landfile );
    fprintf(stderr, "Exiting\n" );
    return -1;
    }
 /*
  *  read in first record and check it out
  */
  if( fread( rec_arr, sizeof(short), 1024, stream[msknum] ) != 1024 )
    {
    fprintf(stderr, "b128_msk_init: read error for first record of file # %d:\n",
         msknum );
    fprintf(stderr, "    '%s'\n", landfile );
    fprintf(stderr, "Exiting\n" );
    return -1;
    }

  if ( endianess() == 1 )
      swapc_bytes((char*)rec_arr, 2, 1024);

  if( rec_arr[0] != 128 || rec_arr[2] != 2048 ||
      rec_arr[3] != -180 || rec_arr[4] != 180 ||
      rec_arr[5] != -90 ||rec_arr[6] != 90 )
    {
    fprintf(stderr, "b128_msk_init: header record does not match expected values\n" );
    fprintf(stderr, "               for mask # %d\n", msknum );
    fprintf(stderr, "\nvalue   expected   read\n" );
    fprintf(stderr,  "0          128      %d\n", rec_arr[0] );
    fprintf(stderr,  "2         2048      %d\n", rec_arr[2] );
    fprintf(stderr,  "3         -180      %d\n", rec_arr[3] );
    fprintf(stderr,  "4          180      %d\n", rec_arr[4] );
    fprintf(stderr,  "5          -90      %d\n", rec_arr[5] );
    fprintf(stderr,  "6           90      %d\n\n", rec_arr[6] );
    fprintf(stderr, "Exiting\n" );
    return -1;
    }
 /*
  *  read the next 64 1024 short records as pointers to the state
  *  of the 1 degree box
  */
  if( fread( ptr_arr[msknum], sizeof(short), NITEMPTR, stream[msknum] ) != 
                                                    NITEMPTR )
    {
    fprintf(stderr, "b128_msk_init: read error for pointer array of file:\n" );
    fprintf(stderr, "    '%s'\n", landfile );
    fprintf(stderr, "               mask # %d\n", msknum );
    fprintf(stderr, "Exiting\n" );
    return -1;
    }

  if ( endianess() == 1 )
      swapc_bytes((char*)ptr_arr[msknum], 2, NITEMPTR);

 /*
  *  initialize some items
  */
  access_cnt[msknum] = 0;  /* current count of mask accesses in the program */
  boxes_active[msknum] = 0;  /* # boxes read into the storage structure */

  return 0;
  }


int b128_msk_get( float lat, float lon, int msknum )
/*******************************************************************

   b128_msk_get

   purpose: For an incoming latitude, longitude find if it is over
      mask or not using the 128 by 128 / degree mask file

   Returns type: int - 1 for mask on, 0 for mask off, -1 for problem

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             lat              I      latitude of point from
                                                -90. to 90.
      float             lon              I      longitudeof point from
                                                -180 to 180.
      int               msknum           I      mask to use ( 0 or 1 )

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       8-Oct-1997      Original development
      W. Robinson       19-Dec-1997     adapt to use 2 mask files

*******************************************************************/
  {
  int land_fl, box_num, box_index, low_access_cnt,
      low_access_index, i, box_wd, box_bit;
  float lat_off, lon_off;
  unsigned short mask[16] = { 0x8000, 0x4000, 0x2000, 0x1000,
                              0x800,  0x400,  0x200,  0x100,
                              0x80,   0x40,   0x20,   0x10,
                              0x8,    0x4,    0x2,    0x1 };
  static int last_box_index[2];
 /*  declare subroutines  */
  int b128_box_num( float, float, float *, float * );
  int b128_wd_bit( float, float, int *, int * );

 /*
  *  check mask # for validity
  */
  if( msknum < 0 || msknum > 1 )
    {
    fprintf(stderr, "b128_msk_init: Fatal error. only mask 0 or 1 can be used\n" );
    fprintf(stderr, "               current mask: %d\n", msknum );
    fprintf(stderr, "Exiting\n" );
    return -1;
    }
 /*
  *  check that b128_msk_init was called - the stream should be non-null
  */
  access_cnt[msknum]++;  /* incriment count of times called */

  if (access_cnt[msknum] >= MAX_PIX_PER_SCENE) {
      printf("b128_msk_get: Pixels Per Scene Count Exceeded.\n");
      exit(1);
  }

  if( stream[msknum] == NULL ) 
    {
    fprintf(stderr, "b128_msk_get: routine b128_msk_init not called yet\n" );
    fprintf(stderr, "              for mask # %d\n", msknum );
    fprintf(stderr, "Exiting\n" );
    land_fl = -1;
    return land_fl;
    }
 /*
  *  get the box number for this lat, lon
  */
  box_num = b128_box_num( lat, lon, &lat_off, &lon_off );
 /*
  *  maybe the box is all mask or non-mask and we can just return that info
  */
  if( *( ptr_arr[msknum] + box_num ) < 2 )
    {
    land_fl = *( ptr_arr[msknum] + box_num );
    }
  else
    {
   /*
    *  the box contains a mix of land, water. see if we're looking at last box 
    */
    if( box_num == last_box_num[msknum] )
      {
      box_index = last_box_index[msknum]; /* index to locally stored box info */
      }
    else
      {
     /*
      *  not latest box, it still may be in another array position
      */
      box_index = -1;
      low_access_cnt = MAX_PIX_PER_SCENE;  /* while searching, record the lowest */
      low_access_index = -1;          /* access count and its location */

      for( i = 0; i < boxes_active[msknum]; i++ )
        {
        if( box_info[msknum][i].box_num == box_num )
          {
         /*
          *  we found the box we need in array
          */
          box_index = i;
          box_info[msknum][i].last_access_cnt = access_cnt[msknum];
          break;
          }
        if( box_info[msknum][i].last_access_cnt < low_access_cnt )
          {
          low_access_cnt = box_info[msknum][i].last_access_cnt;
          low_access_index = i;  /* record least recently used updated */
          }
        }

     /*
      *  if we didn't find the box in the struct array, set one up
      */
      if( box_index == -1 )
        {
       /*
        *  first case: add a box to the array in an empty slot
        */
        if( boxes_active[msknum] < BOX_MAX - 1 )
          {
          boxes_active[msknum]++;
          box_index = boxes_active[msknum] - 1;
          box_info[msknum][box_index].box_ptr = 
                (unsigned short *)malloc( 1024 * sizeof( short ) );
         /*  for debug / info
          printf( "Incrimenting boxes_active to %d\n", boxes_active );
          printf( "          lat: %f   lon: %f\n", lat, lon );
          */
          }
       /*
        *  2nd case is to re-use the least used array location
        */
        else
          {
          box_index = low_access_index;
         /*  for debug / info
          printf( "Re-using box index %d\n", box_index );
          printf( "          lat: %f   lon: %f\n", lat, lon );
          */
          }
       /*
        *  in either case, read in the proper record into the slot
        */
        box_info[msknum][box_index].box_num = box_num;
        box_info[msknum][box_index].last_access_cnt = access_cnt[msknum];
        if( fseek( stream[msknum], 
              ptr_arr[msknum][box_num] * 1024 * sizeof( short ), SEEK_SET )
              != 0 )
          {
          fprintf(stderr, "b128_msk_get: Failure on mask file record seek\n" );
          fprintf(stderr, "       record # %d   mask # %d\n", 
               ptr_arr[msknum][ box_num ], msknum );
          fprintf(stderr, "  Exiting\n" );
          land_fl = -1;
          }
        else
          {
          if( fread( box_info[msknum][box_index].box_ptr, sizeof( short ), 
                   1024, stream[msknum] ) != 1024 )
            {
            fprintf(stderr, "b128_msk_get: Failure on mask record read\n" );
            fprintf(stderr, "       record # %d  mask # %d\n", 
                ptr_arr[msknum][ box_num ], msknum );
            fprintf(stderr, "  Exiting\n" );
            land_fl = -1;
            }

          if ( endianess() == 1 ) 
              swapc_bytes((char*)box_info[msknum][box_index].box_ptr, 2, 1024);


          }
        }
      }
   /*
    *  we have the box, get the exact word and bit to extract
    */
    b128_wd_bit( lat_off, lon_off, &box_wd, &box_bit );
    land_fl = ( ( *( box_info[msknum][ box_index ].box_ptr + box_wd ) &
                                          mask[ box_bit ] ) == 0 ) ? 0 : 1;
    last_box_index[msknum] = box_index;
    last_box_num[msknum] = box_num;
    }
  return land_fl;
  }


int b128_box_num( float lat, float lon, float *lat_off, float *lon_off )
/*******************************************************************

   b128_box_num

   purpose: for the 128 / degree land mask file, find the
      index of the degree box from the lat, lon

   Returns type: int - index (0 rel) of the 1 degree box

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             lat              I      latitude of point from
                                                -90. to 90.
      float             lon              I      longitudeof point from
                                                -180 to 180.
      float *           lat_off          O      all positive latitude used
                                                during call to
                                                b128_wd_bit to get mask
                                                at a 128th of a degree point
      float *           lon_off          O      all positive longitude

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       8-Oct-1997      Original development

*******************************************************************/
  {
  int box_num, lat_index, lon_index;
 /*
  *  make positive lat and lon and check them
  */
  *lat_off = 90. + lat;
  *lon_off = 180. + lon;

  if( *lat_off < 0. ) *lat_off = 0.;
  if( *lat_off > 180. ) *lat_off = 180.;
  if( *lon_off < 0. ) *lon_off = 0.;
  if( *lon_off > 360. ) *lon_off = 360.;
 /*
  *  Take care of case of 90 lat and 180 lon properly
  */
  lat_index = ( *lat_off == 180. ) ? 179 : (int)(*lat_off);
  lon_index = ( *lon_off == 360. ) ? 359 : (int)(*lon_off);
 /*
  *  compute the box #
  */
  box_num = lon_index + 360 * lat_index;

  return box_num;
  }



int b128_wd_bit( float lat_off, float lon_off, int *box_wd, int *box_bit )
/*******************************************************************

   b128_box_num

   purpose: for the 128 / degree land mask file, find the
      index of the degree box from the lat, lon

   Returns type: int - nothing now, later, maybe an error condition

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             lat_off          I      all positive latitude
                                                generated in call to
                                                b128_box_num
      float             lon_off          I      all positive longitude
                                                as above
      int *             box_wd           O      word # in array of
                                                128 sq box to find
                                                desired mask value
      int *             box_bit          O      bit # in word of
                                                desired mask value

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       8-Oct-1997      Original development

*******************************************************************/
  {
  int lat_sub, lon_sub;
  int lon_wd;
  double dumb;

 /*
  *  find distance in to the bit from edge of 128 by 128 array
  */
  lat_sub = (int) (modf( lat_off, &dumb ) * 128 );
  if( lat_off == 180. ) lat_sub = 127;

  lon_sub = (int) (modf( lon_off, &dumb ) * 128 );
  if( lon_off == 360. ) lon_sub = 127;
 /*
  *  get the word in longitude direction and bit and then get linear word
  */
  lon_wd = lon_sub / 16;
  *box_bit = lon_sub - lon_wd * 16;
  *box_wd = lon_wd + lat_sub * 8;

  return 0;
  }
