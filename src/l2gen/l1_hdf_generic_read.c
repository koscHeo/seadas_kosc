#include "l1_hdf_generic_read.h"
#include "l12_proto.h"

#define BANDED     2
#define INTERLACED 3

static int32 format  = 0;
static int32 dformat = INTERLACED;
static char **dname;
//static char  dname[NBANDS][32];
static float *data   = NULL;

#define READ_GLBL_ATTR(nam,ptr) {                                           \
  if(SDreadattr(sd_id,SDfindattr(sd_id,(nam)),(VOIDP)(ptr))){               \
    fprintf(stderr,                                                         \
    "-E- %s line %d: Could not get global attribute, %s.\n", \
    __FILE__,__LINE__,(nam));                                       \
  }                                                                         \
}

#define READ_SDS(nam,ptr,s0,s1,s2,e0,e1,e2) {                           \
  int32 start[3];                                                       \
  int32 edge[3];                                                        \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);                             \
  start[0]=(s0); start[1]=(s1); start[2]=(s2);                          \
  if(SDreaddata(SDselect(sd_id, SDnametoindex(sd_id, (nam))),           \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){                            \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n",         \
    __FILE__,__LINE__,(nam));                                           \
  }                                                                     \
}

void get_l1data( int32	sd_id,	
                 int32	scan,
                 int32  npix,
                 int32  nbands,
		 int32  bindx[],
                 int32  format,
                 l1str *l1rec)
{
  //  int16 angles[NBANDS];
    int32 i,j;
    static int firstCall=1;

    READ_SDS("year"     ,l1rec->year, scan, 0,0, 1,   1,1);
    READ_SDS("day"      ,l1rec->day,  scan, 0,0, 1,   1,1);
    READ_SDS("msec"     ,l1rec->msec, scan, 0,0, 1,   1,1);
    READ_SDS("mside"    ,&(l1rec->mside), scan, 0,0, 1,1,1);
    READ_SDS("detnum"   ,&(l1rec->detnum), scan, 0,0, 1,1,1);
    READ_SDS("tilt"     ,&(l1rec->tilt), scan, 0,0, 1,1,1);
    READ_SDS("longitude",l1rec->lon,  scan, 0,0, 1,npix,1);
    READ_SDS("latitude" ,l1rec->lat,  scan, 0,0, 1,npix,1);
    READ_SDS("solz"     ,l1rec->solz, scan, 0,0, 1,npix,1);
    READ_SDS("sola"     ,l1rec->sola, scan, 0,0, 1,npix,1);
    READ_SDS("senz"     ,l1rec->senz, scan, 0,0, 1,npix,1);
    READ_SDS("sena"     ,l1rec->sena, scan, 0,0, 1,npix,1);

    if (firstCall==1) {
        if ( (dname = (char **) malloc(l1rec->nbands*sizeof(char *)))  == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                   __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }

        for (i=0; i<l1rec->nbands; i++)
            if ( (dname[i] = (char *) malloc(32*sizeof(char)) )  == NULL) {
                printf("-E- %s line %d: Error allocating data space.\n",
                       __FILE__,__LINE__);
                exit(EXIT_FAILURE);
            }
        firstCall = 0;
    }

    if (dformat == INTERLACED) {
        READ_SDS("l1b_data" ,l1rec->Lt,   scan, 0,0, 1,npix,nbands);
    } else {
        for (i=0; i<nbands; i++) {
            READ_SDS(dname[i],data,scan,0,0,1,npix,1);
            for (j=0; j<npix; j++)
                l1rec->Lt[j*nbands+bindx[i]] = data[j];
	}
    }

    /* Read L2 flags */
    READ_SDS("l2_flags",l1rec->flags,scan,0,0, 1,npix,1);
    for (j=0; j<npix; j++) {
        l1rec->stlight[j] = ((l1rec->flags[j] & STRAYLIGHT) > 0);
        l1rec->hilt   [j] = ((l1rec->flags[j] & HILT   ) > 0);
        l1rec->navwarn[j] = ((l1rec->flags[j] & NAVWARN) > 0);
        l1rec->navfail[j] = ((l1rec->flags[j] & ATMFAIL) > 0);
    }

}


int openl1_read_hdf_g(filehandle *file)
{
    int32 npix;
    int32 nscan;
    int32 nbands;
    int32 i;
    char  apb[20];
    char  sensor[20];
    int32 sd_id;
    int32 sds_id;

    /* Open the HDF input file */
    sd_id = SDstart(file->name, DFACC_RDONLY);
    if(sd_id == FAIL){
        fprintf(stderr,"-E- %s line %d: SDstart(%s, %d) failed.\n",
        __FILE__,__LINE__,file->name,DFACC_RDONLY);
        return(HDF_FUNCTION_ERROR);
    }

    /* Read some of the level-1A global attributes. */
    READ_GLBL_ATTR("Pixels per Scan Line" , &npix    );
    READ_GLBL_ATTR("Number of Scan Lines" , &nscan   );
    READ_GLBL_ATTR("Number of Bands"      , &nbands  );
    READ_GLBL_ATTR("Sensor Name"          , sensor   );

    /* Determine which kind of radiance format we have */
    if ( sd_select(sd_id,"l1b_data",&sds_id) == 0 )
        dformat = INTERLACED;
    else
        dformat = BANDED;


    file->npix   = npix;
    file->nscan  = nscan;
    file->nbands = nbands;
    file->sd_id  = sd_id;

    return(LIFE_IS_GOOD);
}


int readl1_hdf_g(filehandle *file, int32 recnum, l1str *l1rec)
{
    static int firstCall = 1;
    int32 sd_id = file->sd_id;

    if (firstCall) {

      firstCall = 0;

      if (dformat == BANDED) {

        int32_t i;
        int32_t *wavelen;

        if ( rdsensorinfo(file->sensorID,file->input->evalmask,"Lambda", (void **) &wavelen) != file->nbands ) {
            printf("-E- %s line %d: Error reading sensor table file\n",
                   __FILE__,__LINE__);
            return(1);
        }
     
        for (i=0; i<file->nbands; i++)
            sprintf(dname[i],"Lt_%3d",wavelen[file->bindx[i]]);

        if ((data = (float *) calloc(file->npix, sizeof(float))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                   __FILE__,__LINE__);
            return(1);
        }
      }

    }

    get_l1data(sd_id, recnum, file->npix, file->nbands, (int32 *) file->bindx, format, l1rec);

    l1rec->sensorID = file->sensorID;
    l1rec->npix     = file->npix;

    return(LIFE_IS_GOOD);
}


int closel1_hdf_g(filehandle *file)
{
   if (SDend(file->sd_id)) {
      fprintf(stderr,"-E- %s line %d: SDend(%d) failed for file, %s.\n",
      __FILE__,__LINE__,file->sd_id,file->name);
      return(HDF_FUNCTION_ERROR);
   }

   return(LIFE_IS_GOOD);
}





