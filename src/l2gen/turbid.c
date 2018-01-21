#include "l12_proto.h"
#define TI_BAD -999.


void tindx_shi(l2str *l2rec, int32_t ipix, float *tindx)
{
    static int32_t mask = LAND | CLOUD;

    static int32_t ib748  =   -1;
    static int32_t ib1240 =   -1;
    static int32_t ib2130 =   -1;

    float dr748;
    float dr1240;
    float dr2130;

    int32_t  ip,ip1,ip2,ipb;

    if (ipix < 0) {
        ip1 = 0;
        ip2 = l2rec->npix-1;
    } else {
        ip1 = ipix;
        ip2 = ipix;
    }

    if (ib748 == -1) {
        if ((ib748 = windex(748.0,l2rec->fwave,l2rec->nbands)) < 0) {
	    printf("tindx_shi: incompatible sensor wavelengths (no 748).\n");
            exit(1);
	}
        if ((ib1240 = windex(1240.0,l2rec->fwave,l2rec->nbands)) < 0) {
	    printf("tindx_shi: incompatible sensor wavelengths (no 1240).\n");
            exit(1);
	}
        if ((ib2130 = windex(2130.0,l2rec->fwave,l2rec->nbands)) < 0) {
	    printf("tindx_shi: incompatible sensor wavelengths (no 2130).\n");
            exit(1);
	}
    }

    for (ip=ip1; ip<=ip2; ip++) {

        ipb = ip*l2rec->nbands+ib748;
        dr748  = (l2rec->Lt[ipb]/l2rec->tg_sol[ipb]/l2rec->tg_sen[ipb] - l2rec->Lr[ipb])/l2rec->Fobar[ib748 ];

        ipb = ip*l2rec->nbands+ib1240;
        dr1240 = (l2rec->Lt[ipb]/l2rec->tg_sol[ipb]/l2rec->tg_sen[ipb] - l2rec->Lr[ipb])/l2rec->Fobar[ib1240];

        ipb = ip*l2rec->nbands+ib2130;
	dr2130 = (l2rec->Lt[ipb]/l2rec->tg_sol[ipb]/l2rec->tg_sen[ipb] - l2rec->Lr[ipb])/l2rec->Fobar[ib2130];

        if (((l2rec->flags[ip] & mask) != 0) || dr1240 <= 0.0 || dr2130 <= 0.0)
	    *tindx++ = TI_BAD;
	else
	  /*  *tindx++ = (dr748/dr1240)*exp(-(492./890.)*log(dr1240/dr2130)); */
	    *tindx++ = (dr748/dr1240)*pow((dr1240/dr2130),-(492./890.));
    }

    return;
}



void tindx_morel(l2str *l2rec, int32_t ipix, float *tindx)
{
    static int32_t ib560 = -1;

    float chl;
    float Rrs;
    float Rrs_lim;
    float X;
    int32_t  ip,ip1,ip2;

    if (ipix < 0) {
        ip1 = 0;
        ip2 = l2rec->npix-1;
    } else {
        ip1 = ipix;
        ip2 = ipix;
    }

    if (ib560 == -1) {
        if ((ib560 = windex(560.0,l2rec->fwave,l2rec->nbands)) < 0) {
	    printf("turbid_morel: incompatible sensor wavelengths (no 560).\n");
            exit(1);
	}
	printf("turbid_morel: using %f nm for 560.\n",l2rec->fwave[ib560]);
    }

    for (ip=ip1; ip<=ip2; ip++) {

        chl = l2rec->chl[ip];
        Rrs = l2rec->Rrs[ip*l2rec->nbands+ib560];

        if (l2rec->mask[ip] || chl <= 0.0 || Rrs <= 0.0) {
	    *tindx++ = TI_BAD;
            l2rec->flags[ip] |= PRODFAIL;
        } else if (chl <= 0.2) {
	    *tindx++ = -100;
	} else {
	    X = log10(MIN(chl,10));
	    Rrs_lim = 0.00331 + X*(0.002122 + X*(0.00031587 - X*0.00023145)); 
            *tindx++ = 100.0 * (Rrs - Rrs_lim)/Rrs;
	}
    }

    return;
}

