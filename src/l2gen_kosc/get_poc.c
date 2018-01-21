#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

static float badval =  BAD_FLT;
static float minval =  0.0;
static float maxval = 10.0*1000;


float poc_stramski_443(float *Rrs, float *wave)
{
    static int firstCall = 1;
    static int convert = -1;
    static float a =  203.2;
    static float b = -1.034;
    static int ib1 = -1;
    static int ib2 = -1;

    float poc = badval;
    float Rrs1 = 0.0;
    float Rrs2 = 0.0;

    if (firstCall) {
        firstCall = 0;
        ib1 = bindex_get(443);
        ib2 = bindex_get(545);
        if (ib2 < 0) ib2 = bindex_get(550);
        if (ib2 < 0) ib2 = bindex_get(555);
        if (ib2 < 0) ib2 = bindex_get(560);
        if (ib1 < 0 || ib2 < 0) {
            printf("-E- %s line %d: required bands not available for Stramski POC\n",
            __FILE__,__LINE__);
            exit(1);
	}
    }
    
    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];

    if (Rrs1 >= 0.0 && Rrs2 > 0.0) {
        Rrs2 = conv_rrs_to_555(Rrs2,wave[ib2]);
        poc = a * pow(Rrs1/Rrs2,b);
    }

    return (poc);
}    


float poc_stramski_490(float *Rrs, float *wave)
{
    static int firstCall = 1;
    static int convert_547 = 0;
    static float a =  308.3;
    static float b = -1.639;
    static int ib1 = -1;
    static int ib2 = -1;

    float poc = badval;
    float Rrs1 = 0.0;
    float Rrs2 = 0.0;

    if (firstCall) {
        firstCall = 0;
        ib1 = bindex_get(490);
        ib2 = bindex_get(545);
        if (ib2 < 0) ib2 = bindex_get(550);
        if (ib2 < 0) ib2 = bindex_get(555);
        if (ib2 < 0) ib2 = bindex_get(560);
        if (ib1 < 0 || ib2 < 0) {
            printf("-E- %s line %d: required bands not available for Stramski POC\n",
            __FILE__,__LINE__);
            exit(1);
	}
    }
    
    Rrs1 = Rrs[ib1];
    Rrs2 = Rrs[ib2];

    if (Rrs1 >= 0.0 && Rrs2 > 0.0) {
        Rrs2 = conv_rrs_to_555(Rrs2,wave[ib2]);
        poc = a * pow(Rrs1/Rrs2,b);
    }

    return (poc);
}    


void get_poc(l2str *l2rec, int prodnum, float prod[])
{
    int32_t  ip;

    for (ip=0; ip<l2rec->npix; ip++) {

      prod[ip] = badval;

      switch (prodnum) {

	case CAT_poc_stramski_443:
          prod[ip] = poc_stramski_443(&l2rec->Rrs[ip*l2rec->nbands],l2rec->fwave);
          break;
 
	case CAT_poc_stramski_490:
          prod[ip] = poc_stramski_490(&l2rec->Rrs[ip*l2rec->nbands],l2rec->fwave);
          break;
 
        default:
          printf("Error: %s : Unknown product specifier: %d\n",__FILE__,prodnum);
          exit(1);
          break;
      };

      if (prod[ip] == badval)
        l2rec->flags[ip] |= PRODFAIL;
      else if (prod[ip] < minval || prod[ip] > maxval) 
	l2rec->flags[ip] |= PRODWARN;

    }

}

