#ifndef _L1_STRUC_H
#define _L1_STRUC_H

#include <stdint.h>

#include "l12_parms.h"
#include "input_struc.h"
#include "hdf.h"

/* Notice: any changes to this structure may require modifications to the */
/* following routines: alloc_l1.c, cpl1rec.c, l1subpix.c.                 */

typedef struct l1_struct {

    int32_t   sensorID;      
    int32_t   length;     /* number of bytes allocated to data block */
    int32_t   npix;
    int32_t   spix;
    int32_t   epix;
    int32_t   dpix;
    int32_t   nscans;
    int32_t   sscan;
    int32_t   escan;
    int32_t   dscan;
    int32_t   nbands;
    int32_t   nbandsir;
    int32_t   ndets;
    int32_t  *bindx;
    int32_t   iscan;
    int32_t   detnum;
    int32_t   mside;

    /* process control paramters */

    int32_t   landMaskOn;
    int32_t   bathMaskOn;
    int32_t   cloudMaskOn;
    int32_t   glintMaskOn;
    int32_t   hiltMaskOn;
    int32_t   stlightMaskOn;
    int32_t   senzMaskOn;
    int32_t   solzMaskOn;


    /* scan-time-specific data */

    double fsol;

    /* scan attributes */

    float tilt;

    /* All parameters below are scan-length dependent */

    /* sensor band-pass-specific data */


    char   *data;       /* points to start of variable-length data block */

    int32_t   *year;
    int32_t   *day;
    int32_t   *msec;
    int32_t   *nobs;
    float  *lon;
    float  *lat;
    float  *solz;
    float  *sola;
    float  *senz;
    float  *sena;
    float  *Lt;
    float  *Lt_unc;

    float  *Ltir;
    float  *Bt;

    float  *prtemp;	/* per pixel */
    float  *delphi;
    float  *csolz;
    float  *csenz;
    int32_t   *pixnum;
    unsigned char  *slot;           /**< slot number                                */
    float  *alpha;
    float  *scattang;

    float  *ws;
    float  *wd;
    float  *mw;
    float  *zw;
    float  *pr;
    float  *oz;
    float  *wv;
    float  *rh;
    float  *no2_tropo;
    float  *no2_strat;
    float  *no2_frac;
    float  *height;
    float  *elev;
    float  alt;     //altitude of sensor
    short  *ancqc;

    short  *ssttype;	/* per pixel - reference type or climatology */
    int32_t   *flags;
    char   *mask;
    char   *hilt;
    char   *cloud;
    char   *glint;
    char   *land;
    char   *swater;
    char   *ice;
    char   *solzmax;
    char   *senzmax;
    char   *stlight;
    char   *absaer;
    char   *navfail;
    char   *navwarn;
    char   *darkpix;
    char   *filter;
    char   *cirrus;

    float  *t_h2o;
    float  *t_o2;
    float  *tg_sol;
    float  *tg_sen;
    float  *t_sol;
    float  *t_sen;
    float  *rhof;
    float  *tLf;
    float  *Lr;
    float  *L_q;
    float  *L_u;
    float  *polcor;
    float  *dpol;
    float  *TLg;
    float  *rhos;
    float  *glint_coef;
    float  *cloud_albedo;
    float  *aerindex;
    float  *sstref;
    float  *sssref;
    float  *sw_n;
    float  *sw_a;
    float  *sw_bb;
    float  *sw_a_avg;
    float  *sw_bb_avg;
    float  *rho_cirrus;

    /* for MERIS L1 */
    int32_t   *pixdet;         /* detector index of pixel */
    float  *radcor;         /* smile correction */

    /* for MERIS L2 */
    int32_t   n_inprods;
    float **in_prods;
    int32_t  *in_flags;

    /* NBANDS dependent variables */
    int32_t     *iwave ;
    float       *fwave ;
    float       *Fo    ;
    float       *Fobar ;
    float       *Fonom ;
    float       *Tau_r ;
    float       *k_oz  ;
    float       *aw    ;
    float       *bbw   ;

    /* for VIIRS unaggregated and superscan */
    int16 scn_fmt;  /* scan format of data, 0 std, else unaggregated */
    float margin_s;  /* extra scan margin beyond actual samples */

    instr  *input;

} l1str;

#endif




