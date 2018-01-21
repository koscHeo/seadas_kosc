#include "l12_proto.h"
#include "l12_parms.h"

#define ConfidentCloudy 0
#define ProbablyCloudy  1
#define ProbablyClear   2
#define ConfidentClear  3

int modis_cloud_mask(l1str *l1rec, int32_t ip)
{
    static int firstCall = 1;
    static int32 sd_id;
    static int32 sds_id; 
    static char* cldmask;
    static int32 start[3]; 
    static int32 end  [3]; 
    static int lastScan = -1;

    char  flag;
    int32 status;

    if (firstCall) {

        char  name   [H4_MAX_NC_NAME]  = "";
        char  sdsname[H4_MAX_NC_NAME]  = "";
        char  file   [FILENAME_MAX] = "";
        int32 rank; 
        int32 nt; 
        int32 dims[H4_MAX_VAR_DIMS]; 
        int32 nattrs;

        if (strcmp(l1rec->input->cldfile,"") == 0) {
      	    printf("-E- %s line %d: no cloud flag file provide.\n",
	       	 __FILE__,__LINE__);
            exit(1);
        }
  
        strcpy(file,l1rec->input->cldfile);
        printf("Loading cloud flag file %s\n",file);

        /* Open the file */
        sd_id = SDstart(file, DFACC_RDONLY);
        if (sd_id == FAIL){
            fprintf(stderr,"-E- %s line %d: error opening %s for reading.\n",
            __FILE__,__LINE__,file);
            exit(1);
        }

        strcpy(sdsname,"Cloud_Mask");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);

         if (dims[2] < l1rec->npix || dims[1] < l1rec->nscans) {
            fprintf(stderr,"-E- %s line %d: incompatible dimensions (%d,%d) vs (%d,%d)\n",
            __FILE__,__LINE__,dims[2],dims[1],l1rec->npix,l1rec->nscans);
            exit(1);
	}

        end[0] = 1;
        end[1] = 1;
        end[2] = dims[2];

        if ((cldmask = (char *) malloc(dims[2]*sizeof(char))) == NULL) {
            fprintf(stderr,"-E- %s line %d: error allocating %d bytes for cloud mask\n",
            __FILE__,__LINE__,dims[2]);
            exit(1);
	}

        firstCall = 0;
    }

    if (lastScan != l1rec->iscan) {
        start[0] = 0;
        start[1] = l1rec->iscan;
        start[2] = 0;
        status = SDreaddata(sds_id, start, NULL, end, (VOIDP) cldmask);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS cloud mask at line %d.\n",
                __FILE__,__LINE__,l1rec->iscan);
            exit(1);
        }
        lastScan = l1rec->iscan;
    }


    flag = (cldmask[l1rec->pixnum[ip]] & 6)/2;

    if (flag < ProbablyClear)
        return(1);
    else
        return(0);

}

int modis_cirrus_mask(l1str *l1rec, int32_t ip)
{
    static int firstCall = 1;
    static int32 sd_id;
    static int32 sds_id; 
    static char* cldmask;
    static int32 start[3]; 
    static int32 end  [3]; 
    static int lastScan = -1;

    char  flag;
    int32 status;

    if (firstCall) {

        char  name   [H4_MAX_NC_NAME]  = "";
        char  sdsname[H4_MAX_NC_NAME]  = "";
        char  file   [FILENAME_MAX] = "";
        int32 rank; 
        int32 nt; 
        int32 dims[H4_MAX_VAR_DIMS]; 
        int32 nattrs;

        if (strcmp(l1rec->input->cldfile,"") == 0) {
      	    printf("-E- %s line %d: no cloud flag file provide.\n",
	       	 __FILE__,__LINE__);
            exit(1);
        }
  
        strcpy(file,l1rec->input->cldfile);
        printf("Loading cloud flag file %s\n",file);

        /* Open the file */
        sd_id = SDstart(file, DFACC_RDONLY);
        if (sd_id == FAIL){
            fprintf(stderr,"-E- %s line %d: error opening %s for reading.\n",
            __FILE__,__LINE__,file);
            exit(1);
        }

        strcpy(sdsname,"Cloud_Mask");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);

        if (dims[2] < l1rec->npix || dims[1] < l1rec->nscans) {
            fprintf(stderr,"-E- %s line %d: incompatible dimensions (%d,%d) vs (%d,%d)\n",
            __FILE__,__LINE__,dims[2],dims[1],l1rec->npix,l1rec->nscans);
            exit(1);
	}

        end[0] = 1;
        end[1] = 1;
        end[2] = dims[2];

        if ((cldmask = (char *) malloc(dims[2]*sizeof(char))) == NULL) {
            fprintf(stderr,"-E- %s line %d: error allocating %d bytes for cloud mask\n",
            __FILE__,__LINE__,dims[2]);
            exit(1);
	}

        firstCall = 0;
    }

    if (lastScan != l1rec->iscan) {
        start[0] = 1;
        start[1] = l1rec->iscan;
        start[2] = 0;
        status = SDreaddata(sds_id, start, NULL, end, (VOIDP) cldmask);
        if (status != 0) {
            printf("-E- %s Line %d:  Error reading SDS cloud mask at line %d.\n",
                __FILE__,__LINE__,l1rec->iscan);
            exit(1);
        }
        lastScan = l1rec->iscan;
    }


    flag = (cldmask[l1rec->pixnum[ip]] & 2);

    if (flag == 0)
        return(1);
    else
        return(0);

}
