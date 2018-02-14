/*---------------------------------------------------------------------*/
/* ipar.c -  functions to compute instantaneous PAR                    */
/*---------------------------------------------------------------------*/

#include "l12_proto.h"
#include "l2prod.h"

static float  unitc = 119.625e8;  /* conversion to einsteins/m^2/s */

#define PARW1 400
#define PARW2 700
#define PARWN (PARW2-PARW1)+1

/*---------------------------------------------------------------------*/
/* get_ipar - computes ipar for each pixel in L2 record                */
/*---------------------------------------------------------------------*/
void get_ipar(l2str *l2rec, float ipar[], initstr *initrec)
{
    static float badval = BAD_FLT;
    static int   firstCall = 1;
    static int32_t  nwave;
    static float *wave;
    static float F0vis[PARWN];
    static float nw = 1.334;

    float *tf;
    float *ta;
    float *tw;
    float *ta_tw;
    //    float Ed0m;
    float Ed0p;
    int32_t  ip, ipb, iw, ib;

    if (firstCall) {

        firstCall = 0;
        if ((wave = (float *) calloc(l2rec->nbands,sizeof(float ))) == NULL) {
            printf("-E- %s line %d : error allocating memory for ipar:get_ipar.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        for (iw=0; iw<l2rec->nbands; iw++){
            wave[iw] = l2rec->fwave[iw];  
        }

        // instantaneous solar irradiance at 1-nm intervals
        for (iw=PARW1; iw<=PARW2; iw++){
            ib = iw - PARW1;
            get_f0_thuillier_ext(iw, 1, &F0vis[ib], initrec->f0rec);
            F0vis[ib] *= l2rec->fsol;
        }

        nwave = MIN(windex(900.,wave,l2rec->nbands)+1,l2rec->nbands);
    }
    if ((tf = (float *) calloc(l2rec->nbands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for ipar:get_ipar.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((ta = (float *) calloc(l2rec->nbands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for ipar:get_ipar.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((tw = (float *) calloc(l2rec->nbands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for ipar:get_ipar.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((ta_tw = (float *) calloc(l2rec->nbands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for ipar:get_ipar.\n",
                __FILE__,__LINE__);
        exit(1);
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        ipar[ip] = 0.;

        if (!l2rec->mask[ip]) {

            // note: redefined ipar as above-surface
            // fresnel_sol(wave,nwave,l2rec->solz[ip],l2rec->ws[ip],tf,1);

            for (iw=0; iw<nwave; iw++) {

                ipb = ip*nwave+iw;

                // atmospheric transmittance (sun to surface)
                ta[iw] = l2rec->tg_sol[ipb] * l2rec->t_sol[ipb];

                // water transmittance (1 - (rho_fres + rho_wc))
                // tw[iw] = tf[iw] - l2rec->rhof[ipb];

                // total transmittance per sensor band
                // ta_tw[iw] = ta * tw;
            }

            for (iw=PARW1; iw<=PARW2; iw++){
                ib = iw - PARW1;
                // Ed0m = F0vis[ib] * l2rec->csolz[ip] * linterp(wave,ta_tw,nwave,(float) iw);
                //ipar[ip] += ((float) iw) * Ed0m / unitc;
                Ed0p = F0vis[ib] * l2rec->csolz[ip] * linterp(wave,ta,nwave,(float) iw);
                ipar[ip] += ((float) iw) * Ed0p / unitc;
            }

            if (!finite(ipar[ip])) {
                ipar[ip]= badval;
                l2rec->flags[ip] |= PRODFAIL;
            }

        } else {
            ipar[ip]= badval;
            l2rec->flags[ip] |= PRODFAIL;
        }
    }
    free(tf);
    free(ta);
    free(tw);
    free(ta_tw);
}

/*
	Energy of one photon = hc/lamda.
	So the number of photons generating Ed (watts/m2/nm) is

	Ed / (hc/lamba) = Ed*lamda /(hc)

	= Ed*lamda/ (6.626e-34 J Sec. x 2.998e8 m/Sec()	
	= Ed*lamda/ (1.9865e-25 J m)	

	So using lamda in nm

	(J/Sec/m2/nm) * nm * 1m/1e9nm / (1.9865e-25 J m x 6.022e23 photons / Ein)

	= (Ein/Sec/m2/nm) / 1.19625e8.	
 */
