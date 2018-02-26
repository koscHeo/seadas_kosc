#include "l12_proto.h"
#include "atrem_corl1.h"

/* ========================================================================== */
/* module atmocor1() - computes pre-aerosol atmospheric components            */
/*                                                                            */
/* Written By: B. A. Franz, SAIC, SeaWiFS Project, February 1998              */
/* Conversion to C, May 2006                                                  */
/* ========================================================================== */


/* ------------------------------------------------------------------- */
/* correction factor to remove oxygen absorption from La(765)          */
/* ------------------------------------------------------------------- */
float oxygen_aer(float airmass)
{
    /* base case: m80_t50_strato  visibility (550nm ,0-2km):25km */
    static float a[] = {-1.0796,9.0481e-2,-6.8452e-3};
    return(1.0+pow(10.0,a[0]+airmass*a[1]+airmass*airmass*a[2]));
}


/* ------------------------------------------------------------------- */
/* correction factor to replace oxygen absorption to Lr(765)           */
/* ------------------------------------------------------------------- */
float oxygen_ray(float airmass)
{
    /* base case is the 1976 Standard atmosphere without aerosols */
    static float a[] = {-1.3491, 0.1155, -7.0218e-3};
    return(1.0/(1.0+pow(10.0,a[0]+airmass*a[1]+airmass*airmass*a[2])));
}	


/* ------------------------------------------------------------------- */
/* main function, loads various quanitities into l1 record for 1 pixel */
/* ------------------------------------------------------------------- */

void atmocor1(l1str *l1rec, int32_t ip)
{
    static float  p0 = STDPR;
    static float  pi = PI;

    int32_t  sensorID = l1rec->sensorID;
    int32_t  nwave    = l1rec->nbands;

    float solz = l1rec->solz[ip];
    float senz = l1rec->senz[ip];
    float raz  = l1rec->delphi[ip];
    float mu0  = l1rec->csolz[ip];
    float mu   = l1rec->csenz[ip];

    float ws   = l1rec->ws[ip];
    float pr   = l1rec->pr[ip];
    float oz   = l1rec->oz[ip];
    float wv   = l1rec->wv[ip];

    float no2_tropo = l1rec->no2_tropo[ip];
    float no2_strat = l1rec->no2_strat[ip];
    float no2_frac = l1rec->no2_frac[ip];

    int32_t  *wave  = l1rec->iwave;
    float *Fo    = l1rec->Fo;
    float *Tau_r = l1rec->Tau_r;

    instr *input = l1rec->input;

    float zero = 0.0;
    int32_t  ib765 = -1;
    float airmass, A;
    float a_o2;
    float glint_coef_q;
    float glint_coef_u;
    int32_t  ib, ipb;
    static float *rhot, *tg_tot;
    static int firstCall=1;
    float scaleRayleigh;

    if (firstCall) {
        if ( (rhot = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to rhot\n");
            exit(FATAL_ERROR);
        }
        if ( (tg_tot = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to tg_tot\n");
            exit(FATAL_ERROR);
        }
        firstCall = 0;
    }

    airmass = 1.0/mu0 + 1.0/mu;


    /* Initialize output values */
    for (ib=0; ib<nwave; ib++) {

        ipb = ip*nwave + ib;

        l1rec->Lr  [ipb] = 0.0;
        l1rec->L_q [ipb] = 0.0;
        l1rec->L_u [ipb] = 0.0;
        l1rec->tLf [ipb] = 0.0;
        l1rec->tg_sol[ipb] = 1.0;
        l1rec->tg_sen[ipb] = 1.0;
        l1rec->t_h2o   [ipb] = 1.0;
        l1rec->t_o2    [ipb] = 1.0;
        l1rec->t_sol   [ipb] = exp(-0.5*pr/p0*Tau_r[ib]/mu0);
        l1rec->t_sen   [ipb] = exp(-0.5*pr/p0*Tau_r[ib]/mu );

        /* Copy TOA radiance to temp var, eliminate missing bands
         * Copy surface reflectance to temp. var for atrem calculation */
        rhot [ib] = pi * l1rec->Lt[ipb]/l1rec->Fo[ib]/mu0;
        //rhot [ib] = l1rec->rhos[ipb];

        if (input->oxaband_opt > 0 && wave[ib] == 765) {
            ib765 = ib;
            l1rec->t_o2[ipb] = 1.0/oxygen_aer(airmass);
        }

        if (sensorID == SEAWIFS) {
            /* this is for rhos only, but effects seawifs land products and cloud   */
            /* will modify get_rhos() to use gaseous transmittance in future update */
            l1rec->t_h2o[ipb] = water_vapor(ib,wv,airmass);
        }

    }
    /* For Atrem:
     * An important refinement we need to make is to separate the ATREM-derived gas
transmittance into separate sun to ground (tg_sol) and ground to sensor (tg_sen)
paths.  Doing this correctly, and being able to account for aircraft altitude,
would seem to require significant modification of the ATREM code.

As a quick and dirty solution, that should be close for spacecraft altitude, we
can effectively just weight by the airmass (slant path distance), I think.

tg_tot = atrem()

mu  = cos(theta)        ; theta = sensor zenith
mu0 = cos(theta0)       ; theta0 = solar zenith
airmass = 1/mu0 + 1/mu  ; total "L-shaped" slant path

A = -ln(tg_tot)/airmass     ; effective optical depth

tg_sol = tg_tot * exp(-A/(1/mu0))
tg_sen = tg_tot * exp(-A/(1/mu ))

It's based on the notion that
tg_tot ~ exp(-A m)
tg_sol ~ exp(-A/mu0)
tg_sen ~ exp(-A/mu).


Currently tg_sol and tg_sen are computed in atmocor1.c.  Airmass terms are there
as well.  If possible, it would be better to call atrem() there, replacing the
current call to gaseous_transmittance(), e.g.:

 //gaseous transmittance

if (gas_opt == 16) {
     tg_tot = atrem()
     for (ib=0; ib<nwave; ib++) {
         ipb = ip*nwave+ib;
     for ib=0 ...
     tg_sol[ipb] = tg_tot[ipb] * exp((1/mu0)/airmass)
     tg_sen[ipb] = tg_tot[ipb] * exp((1/mu )/airmass)
} else {
     ipb = ip*nwave;
gaseous_transmittance(input->gas_opt,sensorID,input->evalmask,nwave,mu0,mu,oz,wv,no2_tropo,no2_strat,
no2_frac,l1rec->tg_sol[ipb],&l1rec->tg_sen[ipb]);
}

You can define ATREM in l12_parms.h.

#define O3_BIT  1
#define CO2_BIT 2
#define NO2_BIT 4
#define H2O_BIT 8
#define ATREM 16


If this is possible, you wouldn't have to mess with any other places where gas
transmittance is applied.

Longer-term, we'd want to determine the rigorous separation of tg_sol and
tg_sen, accounting for aircraft altitude and water vapor profile.

-- bryan
*/
    if ((input->gas_opt & ATREM_BIT) != 0) {

        get_atrem_cor (l1rec->sensorID, l1rec, ip, rhot, tg_tot);
        for (ib=0; ib<nwave; ib++) {
            ipb = ip*nwave + ib;
            A = -log(tg_tot[ib])/airmass;     //effective optical depth
            l1rec->tg_sol[ipb] = exp(-A/mu0);
            l1rec->tg_sen[ipb] = exp(-A/mu );
//            l1rec->tg_sol[ipb] = tg_tot[ib];
//            l1rec->tg_sen[ipb] = 1.0;

        //    printf("RJH:ATREM: %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d \n",ib+1,ip,l1rec->Lt[ipb],Fo[ib],mu0, rhot[ib],tg_tot[ib],debug_atrem.rp94,debug_atrem.r1p14,debug_atrem.cst1,debug_atrem.cst2,debug_atrem.cst3,debug_atrem.cst4,debug_atrem.cst5,debug_atrem.cst6,debug_atrem.jac,debug_atrem.jbc);
        }
    }

        /* gaseous transmittance */

        ipb = ip*nwave;
        gaseous_transmittance(input->gas_opt,sensorID,input->evalmask,nwave,mu0,mu,oz,wv,no2_tropo,no2_strat, no2_frac,
                &l1rec->tg_sol[ipb],&l1rec->tg_sen[ipb]);

    /* white-cap radiances at TOA */

    ipb = ip*nwave;
    whitecaps(sensorID,input->evalmask,nwave,ws,input->wsmax,&l1rec->rhof[ipb]);
    for (ib=0; ib<nwave; ib++) {
        ipb = ip*nwave+ib;
        l1rec->tLf[ipb] = l1rec->rhof[ipb]*l1rec->t_sen[ipb]*l1rec->t_sol[ipb]*Fo[ib]*mu0/pi;
    }

    /* Rayleigh scattering */

    ipb = ip*nwave;

    if (sensorID != AVHRR) {

      rayleigh(sensorID,input->evalmask,nwave,input->pol_opt,solz,senz,raz,Tau_r,Fo,pr,ws,
	       &l1rec->Lr[ipb],&l1rec->L_q[ipb],&l1rec->L_u[ipb]);             
    }

    if (input->oxaband_opt > 0 && ib765 > -1) {
        a_o2 = oxygen_ray(airmass);
        l1rec->Lr [ipb+ib765] *= a_o2;
        l1rec->L_q[ipb+ib765] *= a_o2;
        l1rec->L_u[ipb+ib765] *= a_o2;
    }

    //Scale by the altitude of the sensor, assuming height of atmosphere=100km

    if (l1rec->alt > 0) {
        scaleRayleigh = 1.0 - exp(-l1rec->alt/10); // Assume 10km is e-folding height

        for (ib=0; ib<nwave; ib++) {

            l1rec->Lr [ipb+ib] *= scaleRayleigh;
            l1rec->L_q[ipb+ib] *= scaleRayleigh;
            l1rec->L_u[ipb+ib] *= scaleRayleigh;
        }
    }

    /* glint coefficients and approximate glint radiances */
    /* also add glint to polarization components          */

    /* for avhrr, l1_aci_hdf.c calls avhrrsub5h.f which calls getglint */
    if (sensorID != AVHRR) {

    getglint_iqu_(&senz,&solz,&raz,&ws,&zero,
                  &l1rec->glint_coef[ip],&glint_coef_q,&glint_coef_u);
    //printf("glint_coef=%lf senz=%lf solz=%lf raz=%lf ws=%lf zero=%lf\n",l1rec->glint_coef[ip],senz,solz,raz,ws,zero);
    
    for (ib=0; ib<nwave; ib++) {
        ipb = ip*nwave+ib;
        l1rec->TLg[ipb]  = l1rec->glint_coef[ip] * exp(-(Tau_r[ib]+0.1)*airmass)*Fo[ib];
        l1rec->L_q[ipb] += (glint_coef_q * l1rec->TLg[ipb]);
        l1rec->L_u[ipb] += (glint_coef_u * l1rec->TLg[ipb]);
    }
//    printf("Tau_r=%lf airmass=%lf Fo=%f\n",Tau_r[10],airmass,Fo[10]);
    ipb = ip*nwave + 10;
    }

}


