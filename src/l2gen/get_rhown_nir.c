#include "l12_proto.h"

static const float g1  = 0.0949;
static const float g2  = 0.0794;


/* ---------------------------------------------------------------------- */
/* Convert Rrs[0+] to Rrs[0-]                                             */
/* ---------------------------------------------------------------------- */
float above_to_below(float Rrs)
{    
    return(Rrs / (0.52 + 1.7 * Rrs));
}

/* ---------------------------------------------------------------------- */
/* Convert Rrs[0-] to Rrs[0+]                                             */
/* ---------------------------------------------------------------------- */
float below_to_above(float Rrs)
{    
    return(Rrs * 0.52/(1 - 1.7 * Rrs));
}

/* ------------------------------------------------------------------- */
/* Description:                                                        */
/*	This computes the normalized water-leaving reflectances        */
/*	at the CZCS 670 channel.                                       */
/*                                                                     */
/* Outputs:                                                            */
/*	rhown(670)                                                     */
/*                                                                     */
/* Algorithm Author:                                                   */
/*      Sean Bailey, Futuretech Corporation                            */
/*      based largely on the work of R. Arnone and R. Stumpf           */
/*      using eta relationship from Z. Lee                             */
/* ------------------------------------------------------------------- */

void rhown_red(char *fqfile, float chl, float aw[], float bbw[], float Rrs[],
               float wave[], int32_t nwave, int32_t ib_red, float rhown[])
{
    static int firstCall = 1;
    static int32_t ib5;
    static float eta = 0.5;
    float static chl_min = 0.2;
    float static chl_max = 30.;

    float aw_red, apg_red, bb5;
    float a, bb, chl_in;
    float salbedo;
    float Rrs_red;
    float Rrs5;
    float foq[nwave];

    if (firstCall) {
        ib5 = windex(555,wave,nwave);
        if (fabs(555 - wave[ib5]) > 15){
            printf("%s line %d: can't find reasonable green band\n",__FILE__,__LINE__);
            printf("looking for 555, found %f\n",wave[ib5]);
            exit(EXIT_FAILURE);
        }

        firstCall = 0;
    }

    chl_in = MAX(MIN(chl,chl_max),chl_min);

    aw_red  = aw [ib_red];

    Rrs5 = Rrs[ib5];

    if (Rrs5 <= 0.0)
        Rrs5 = 0.001;

    foqint_morel(fqfile,wave,nwave,0.0,0.0,0.0,chl_in,foq);

    /* Compute total absorption at 670 */
    // NOMAD fit of apg670 to chl
    apg_red = exp(log(chl) * 0.9389 - 3.7589);
    apg_red = MIN(MAX(apg_red,0.0),0.5);
    a = aw_red + apg_red;

    /* Compute backscatter at 550 from Carder/Lee */
    //bb5 = (-0.00182 + 2.058*Rrs5 + bbw5);
    bb5 = (-0.00182 + 2.058*Rrs5);

    /* Translate bb to NIR wavelength */
    bb = bb5 * pow((wave[ib5]/wave[ib_red]),eta) + bbw[ib_red];

    /* Remote-sensing reflectance */
    salbedo = bb / (a + bb) ;
    Rrs_red = foq[ib_red] * salbedo;
    /* Normalized water-leaving reflectance */
    Rrs_red =  below_to_above(Rrs_red);
    rhown[ib_red] = PI*Rrs_red;
}



/* ------------------------------------------------------------------- */
/* Description:                                                        */
/*	This computes the normalized water-leaving reflectances        */
/*	for NIR bands using a bio-optical model and assuming that the  */
/*      NIR absorption is due to the water only.                       */
/*                                                                     */
/* Algorithm Author:                                                   */
/*	Sean Bailey, Futuretech Corporation                            */
/*      based largely on the work of R. Arnone and R. Stumpf           */
/*      using eta relationship from Z. Lee                             */
/* ------------------------------------------------------------------- */

void rhown_nir(char *fqfile,float chl, float aw[], float bbw[], float Rrs[], float wave[],
	       int32_t nwave, float solz, float senz, float phi, 
	       int32_t nir_s, int32_t nir_l, float rhown[])
{
    static int firstCall = 1;
    static int32_t ib2;
    static int32_t ib5;
    static int32_t ib6;

    float a6, aw6, apg6, bbp6;
    float a, bb, eta;
    float static chl_min = 0.2;
    float static chl_max = 30.;
    float foq[nwave];
    float salbedo;
    float Rrs_nir;
    float Rrs2, Rrs5, Rrs6, Rrs6_star;
    int32_t  ib;
    int32_t  dnir = MAX(nir_l-nir_s,1);

    if (firstCall) {
        ib6 = windex(670,wave,nwave);
        if (fabs(670-wave[ib6]) > 50){
            printf("%s line %d: can't find reasonable red band\n",__FILE__,__LINE__);
            printf("looking for 670, found %f\n",wave[ib6]);
            exit(EXIT_FAILURE);
        }

        ib5 = windex(555,wave,nwave);
        if (fabs(555 - wave[ib5]) > 15){
            printf("%s line %d: can't find reasonable green band\n",__FILE__,__LINE__);
            printf("looking for 555, found %f\n",wave[ib5]);
            exit(EXIT_FAILURE);
        }
        ib2 = windex(443,wave,nwave);
        if (fabs(443 - wave[ib2]) > 5){
            printf("%s line %d: can't find reasonable blue band\n",__FILE__,__LINE__);
            printf("looking for 443, found %f\n",wave[ib2]);
            exit(EXIT_FAILURE);
        }
        firstCall = 0;
    }

    Rrs2 = Rrs[ib2];
    Rrs5 = Rrs[ib5];
    Rrs6 = Rrs[ib6];

    if ( Rrs6 <= 0.0) {
        for ( ib=nir_s; ib<=nir_l; ib+=dnir)
            rhown[ib] = 0.0;
        return;
    }
    
    chl = MAX(MIN(chl,chl_max),chl_min);

    // NOMAD fit of apg670 to chl
    apg6 = exp(log(chl) * 0.9389 - 3.7589);
    apg6 = MIN(MAX(apg6,0.0),0.5);

    /* Compute total absorption at 670 */
    aw6  = aw [ib6];
    a6 = aw6 + apg6;

    /* Go below... */   
    Rrs2 = above_to_below(Rrs2);
    Rrs5 = above_to_below(Rrs5);
    Rrs6 = above_to_below(Rrs6);

    //foqint_morel(wave,nwave,0.0,0.0,0.0,chl_in,foq);
    foqint_morel(fqfile, wave,nwave,solz,senz,phi,chl,foq);

    /* Compute the backscatter slope ala Lee */
    eta = 0.0;
    if (Rrs5 > 0.0 && Rrs2 > 0.0) {
        eta = 2.0 * (1. - 1.2 * exp(-0.9*(Rrs2/Rrs5)));
        eta = MIN(MAX(eta,0.0),1.0);
    }

    /* Compute total backscatter at 670 */
    Rrs6_star = Rrs6/foq[ib6];
    bbp6 = (Rrs6_star*a6/(1. - Rrs6_star)) - bbw[ib6];
    
    /* Compute normalized water-leaving reflectance at each NIR wavelength */
    for ( ib=nir_s; ib<=nir_l; ib+=dnir) {

        if (ib == ib6) 
            a = a6;
        else
            a = aw[ib];

        /* Translate bb to NIR wavelength */
	bb = bbp6 * pow((wave[ib6]/wave[ib]),eta) + bbw[ib];

        /* Remote-sensing reflectance */
	salbedo = bb / (a + bb) ;
        Rrs_nir = foq[ib6] * salbedo;
        /* Normalized water-leaving reflectance */
	Rrs_nir =  below_to_above(Rrs_nir);
        rhown[ib] = PI*Rrs_nir;
    }
}


/* ------------------------------------------------------------------- */
/* get_rhown_nir(): calls the appropriate function to compute the      */
/* normalized water-leaving reflectance contribution in the NIR        */
/*                                                                     */
/* B.A. Franz, OBPG 25 January 2006                                    */
/* ------------------------------------------------------------------- */

void get_rhown_eval(char *fqfile, float Rrs[], float wave[], int32_t nir_s, int32_t nir_l,
               int32_t nwave, float aw[], float bbw[], float chl, 
	       float solz, float senz, float phi, float rhown[])
{

    if (wave[nir_l] < MAXWAVE_VIS) {
        rhown_red(fqfile, chl,aw,bbw,Rrs,wave,nwave,nir_l,rhown);
    } else {
        rhown_nir(fqfile,chl,aw,bbw,Rrs,wave,nwave,solz,senz,phi,nir_s,nir_l,rhown);
    }

    return;
}

