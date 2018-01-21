/* =========================================================== */
/* Module gaseous_transmittance.c                              */
/*                                                             */
/* Computes sensor-specific transmittance through various      */
/* atmospheric gases.                                          */
/*                                                             */
/* B. Franz, NASA/OBPG, July 2006                              */
/* =========================================================== */

#include "l12_proto.h"

void ozone_transmittance(int32_t sensorID, int32_t evalmask, int32_t nwave, 
               float mu0, float mu, float ozone, 
	       float t_oz_sol[], float t_oz_sen[])
{
    static float *k_oz = NULL;
    float  tau_oz;
    int32_t   iw;
    
    if (k_oz == NULL) {
        rdsensorinfo(sensorID,evalmask,"k_oz",  (void **) &k_oz);
    }
    
    for (iw=0; iw<nwave; iw++) {   
        tau_oz = ozone * k_oz[iw];
        t_oz_sol[iw] = exp( -(tau_oz/mu0) );
        t_oz_sen[iw] = exp( -(tau_oz/mu ) );
    }
}
	    

	
void co2_transmittance(int32_t sensorID, int32_t evalmask, int32_t nwave,
               float mu0, float mu, 
	       float t_co2_sol[], float t_co2_sen[])
{
    static float *t_co2 = NULL;
    int32_t iw;
    
    if (t_co2 == NULL) {
        rdsensorinfo(sensorID,evalmask,"t_co2",  (void **) &t_co2);
    }
    
    for (iw=0; iw<nwave; iw++) {    
        t_co2_sol[iw] = pow(t_co2[iw],1.0/mu0);
        t_co2_sen[iw] = pow(t_co2[iw],1.0/mu );
    }
}

void no2_transmittance(int32_t sensorID, int32_t evalmask, int32_t nwave,
	       float mu0, float mu, 
               float no2_tropo, float no2_strat, float no2_frac,
	       float t_no2_sol[], float t_no2_sen[])
{
    static float *a_no2 = NULL;

    float a_285, a_225;
    float tau_to200;
    float no2_tr200;
    float a0, a1, a2, a3;
    float sec  = 1.0/mu;
    float sec0 = 1.0/mu0;
    int32_t iw;
    
    if (a_no2 == NULL) {
        rdsensorinfo(sensorID,evalmask,"k_no2",  (void **) &a_no2);
    }

    if (no2_tropo > 0.0) 
        /* compute tropo no2 above 200m (Z.Ahmad)    
        no2_tr200 = exp(12.6615 + 0.61676*log(no2_tropo));
           new, location-dependent method */
        no2_tr200 = no2_frac * no2_tropo;
    else
        no2_tr200 = 0.0;


    for (iw=0; iw<nwave; iw++) {    

      if (a_no2[iw] > 0.0) {

        a_285 = a_no2[iw] * (1.0 - 0.003*(285.0-294.0));
        a_225 = a_no2[iw] * (1.0 - 0.003*(225.0-294.0)); 

        tau_to200 = a_285*no2_tr200 + a_225*no2_strat;

        t_no2_sol[iw] = exp(-(tau_to200*sec0));
        t_no2_sen[iw] = exp(-(tau_to200*sec ));

      } else {

        t_no2_sol[iw] = 1.0;
        t_no2_sen[iw] = 1.0;

      }
    }
}


void h2o_transmittance(int32_t sensorID, int32_t evalmask, int32_t nwave,
	       float mu0, float mu, float wv,
	       float t_h2o_sol[], float t_h2o_sen[])
{
    static float *a_h2o = NULL;
    static float *b_h2o = NULL;
    static float *c_h2o = NULL;
    static float *d_h2o = NULL;
    static float *e_h2o = NULL;
    static float *f_h2o = NULL;
    static float *g_h2o = NULL;

    float t_h2o;
    int32_t iw;
    
    if (a_h2o == NULL) {
        rdsensorinfo(sensorID,evalmask,"a_h2o",  (void **) &a_h2o);
        rdsensorinfo(sensorID,evalmask,"b_h2o",  (void **) &b_h2o);
        rdsensorinfo(sensorID,evalmask,"c_h2o",  (void **) &c_h2o);
        rdsensorinfo(sensorID,evalmask,"d_h2o",  (void **) &d_h2o);
        rdsensorinfo(sensorID,evalmask,"e_h2o",  (void **) &e_h2o);
        rdsensorinfo(sensorID,evalmask,"f_h2o",  (void **) &f_h2o);
        rdsensorinfo(sensorID,evalmask,"g_h2o",  (void **) &g_h2o);
    }
    
    for (iw=0; iw<nwave; iw++) {    
        t_h2o = a_h2o[iw] + wv*(b_h2o[iw] + wv*(c_h2o[iw] + wv*(d_h2o[iw] 
              + wv*(e_h2o[iw] + wv*(f_h2o[iw] + wv*g_h2o[iw])))));
        t_h2o_sol[iw] = pow(t_h2o,1.0/mu0);
        t_h2o_sen[iw] = pow(t_h2o,1.0/mu );
    }
}



void gaseous_transmittance(int gasmask, int32_t sensorID, int32_t evalmask, int32_t nwave, float mu0, float mu, 
               float ozone, float wv, float no2_tropo, float no2_strat, float no2_frac,
	       float t_gas_sol[], float t_gas_sen[])
{
    float *t_sol;
    float *t_sen;
    int32_t iw;

    if ((t_sol = (float *) calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for t_sol in gas_trans.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((t_sen = (float *) calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for t_sen in gas_trans.\n",
                __FILE__,__LINE__);
        exit(1);
    }

    if ((gasmask & O3_BIT) != 0) {
        ozone_transmittance(sensorID,evalmask,nwave,mu0,mu,ozone,t_sol,t_sen);
        for (iw=0; iw<nwave; iw++) {    
            t_gas_sol[iw] *= t_sol[iw];
            t_gas_sen[iw] *= t_sen[iw];
        }
    }
    
    if ((gasmask & CO2_BIT) != 0) {
        co2_transmittance(sensorID,evalmask,nwave,mu0,mu,t_sol,t_sen);
        for (iw=0; iw<nwave; iw++) {    
            t_gas_sol[iw] *= t_sol[iw];
            t_gas_sen[iw] *= t_sen[iw];
        }
    }

    if ((gasmask & NO2_BIT) != 0) {
        no2_transmittance( sensorID, evalmask, nwave, mu0, mu, no2_tropo,
          no2_strat, no2_frac, t_sol, t_sen );
        for (iw=0; iw<nwave; iw++) {
            t_gas_sol[iw] *= t_sol[iw];
            t_gas_sen[iw] *= t_sen[iw];
        }
    }
    
    if ((gasmask & H2O_BIT) != 0) {
        h2o_transmittance(sensorID,evalmask,nwave,mu0,mu,wv,t_sol,t_sen);
        for (iw=0; iw<nwave; iw++) {    
            t_gas_sol[iw] *= t_sol[iw];
            t_gas_sen[iw] *= t_sen[iw];
        }
    }
    free(t_sol);
    free(t_sen);
}
