/*
 * get_hab.c
 *
 * Harmful Algal Blooms
 *  Created on: Aug 31, 2015
 *      Author: Rick Healy (richard.healy@nasa.gov)
 */
#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include <sensorInfo.h>
#include "mph_flags.h"
/*
 * Harmful Algal Bloom Indexes and related functions
 * for flagging products and cloud masking
 * R. Healy 9/1/2015 (richard.healy@nasa.gov)
 *
 * Cyanobacteria Index
 * Wynne,T.T.; Stumpf, R.P.; 2013. Spatial and Temporal Patterns in the Seasonal Distribution of
            Toxic Cyanobacteria in Western Lake Erie from 2002–2014,Toxins 2015, 7, 1649-1663; doi:10.3390/toxins7051649

   Maximum Chlorophyll Index
    C.E. Binding ⁎, T.A. Greenberg, R.P. Bukata, The MERIS Maximum Chlorophyll Index; its merits and limitations for inland water
                algal bloom monitoring, JGLR-00579
 *
 */

static unsigned char *flags_habs = NULL;

void get_habs_ci(l2str *l2rec,  l2prodstr *p, float ci[]);
void get_habs_mph(l2str *l2rec,  l2prodstr *p, float mph_chl[]);
void get_habs_mph_flags(l2str *l2rec,  l2prodstr *p, float flags[]);
void get_habs_cldmask_meris(l2str *l2rec, float cld[]);
void get_habs_cldmask_modis(l2str *l2rec, float cld[]);
char get_cloudmask_meris(l1str *l1rec, int32_t ip);
char get_cloudmask_modis(l1str *l1rec, int32_t ip);
char get_cldmask(l1str *lrec,int32_t ip);

void get_habs_ci(l2str *l2rec,  l2prodstr *p, float ci[])
{


    int ib0,ib1,ib2,ib3,firstCall=1;
    float wav0,wav1,wav2,wav3, nonci,citmp,fac;
    int ip,ipb;

    switch (l2rec->sensorID) {
        case HMODISA:
        case HMODIST:
            fac = 1.3;
            break;
        default:
            fac = 1.0;
    }
    switch (p->cat_ix ) {

    case CAT_CI_stumpf:
        // Cyanobacteria Index
        // Wynne, Stumpf algorithm 2013
        wav1 = 665;
        wav2 = 681;
        wav3 = 709;
            break;

        case CAT_CI_cyano:
        case CAT_CI_noncyano:
            // Cyanobacteria Index
            // Wynne, Stumpf algorithm 2013
            wav0 = 620;
            wav1 = 665;
            wav2 = 681;
            wav3 = 709;
            ib0 = bindex_get(wav0);
        break;

    case CAT_MCI_stumpf:

        wav1 = 681;
        wav2 = 709;
        wav3 = 754;
        break;

    default:
        printf("HABS_CI: Hmm, something's really messed up.\n");
        exit(1);
    }

    ib1 = bindex_get(wav1);
    ib2 = bindex_get(wav2);
    ib3 = bindex_get(wav3);

    if (ib1 < 0 || ib2< 0 || ib3 < 0) {
        printf("(M)CI_stumpf: incompatible sensor wavelengths for this algorithm\n");
        exit(1);
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        ipb = l2rec->nbands*ip;

//        if(l2rec->Rrs[ipb+ib1] <= 0.0 || l2rec->Rrs[ipb+ib2] <= 0.0 || l2rec->Rrs[ipb+ib3] <= 0.0) {
        if(l2rec->rhos[ipb+ib1] < 0.0 || l2rec->rhos[ipb+ib2] < 0.0 || l2rec->rhos[ipb+ib3] < 0.0) {
            ci[ip] = BAD_FLT;
            l2rec->flags[ip] |= PRODFAIL;
        } else {
            switch (p->cat_ix ) {

                case CAT_CI_stumpf:
                    ci[ip] =  (l2rec->rhos[ipb+ib3] - l2rec->rhos[ipb+ib1])*(wav2 - wav1)/(wav3 - wav1)
                            - (l2rec->rhos[ipb+ib2] - l2rec->rhos[ipb+ib1]) ;
                    //Scale for geotiff's
//                    ci[ip] = (250/2.5) *
//                             (4 + log10(ci[ip]*fac));
                    break;
                case CAT_CI_cyano:
                case CAT_CI_noncyano:
                    if(l2rec->rhos[ipb+ib0]<0){
                        ci[ip] = BAD_FLT;
                        l2rec->flags[ip] |= PRODFAIL;
                    } else {
                        //0-(%{rho681} - %{rho665} + (%{rho665}-%{rho709})*(681-665)/(709-665))
                        ci[ip] =  fac*((l2rec->rhos[ipb+ib3] - l2rec->rhos[ipb+ib1])*(wav2 - wav1)/(wav3 - wav1)
                                - (l2rec->rhos[ipb+ib2] - l2rec->rhos[ipb+ib1]));
//                        ci[ip] = (250/2.5) *
//                                 (4 + log10(ci[ip])) +0.5;

                        //%{numreal+1} = %{rho665} - %{rho620} + (%{rho620}-%{rho681})*(665-620)/(681-620)
                        nonci  =  l2rec->rhos[ipb+ib1] - l2rec->rhos[ipb+ib0]  + (l2rec->rhos[ipb+ib0] - l2rec->rhos[ipb+ib2])*(wav1 - wav0)/(wav2 - wav0);
                        if (p->cat_ix == CAT_CI_noncyano) {
                            if (nonci >= 0) ci[ip] = 0;
                        } else {
                        if (nonci < 0) ci[ip] = 0;
                        }
                    }
                    break;
               case CAT_MCI_stumpf:
                   ci[ip] =  fac*( l2rec->rhos[ipb+ib2] - l2rec->rhos[ipb+ib1]
                                  - (l2rec->rhos[ipb+ib3] - l2rec->rhos[ipb+ib1])*(wav2 - wav1)/(wav3 - wav1) );
                    //Scale for geotiff's
//                    ci[ip] = (250/4.0) *
//                             (4 + log10(ci[ip]))+0.5;
                    break;
               default:
                     ci[ip] = BAD_FLT;
                    break;
            }

            //if (ci[ip]<0) ci[ip]=0;
            //         %{numreal} = 0-(%{rho681} - %{rho665} + (%{rho665}-%{rho709})*(681-665)/(709-665))
            //         %{numreal} = (%{rho709} - %{rho681} + (%{rho681}-%{rho754})*(709-681)/(754-681))


        }

    }
}

/*
 * Maximum Peak Height of chlorophyll for MERIS
 *
 * R. Healy (9/1/2015) richard.healy@nasa.gov
 *
 * Mark William Matthews , Daniel Odermatt
 * Remote Sensing of Environment 156 (2015) 374–382
 *
 *
 */
void get_habs_mph(l2str *l2rec,  l2prodstr *p, float chl_mph[])
{
    float wav6,wav7,wav8,wav9,wav10,wav14;
    int ip,ipb,ib6,ib7,ib8,ib9,ib10,ib14;
    float Rmax0,Rmax1,wavmax0,wavmax1,ndvi;
    float sipf,sicf,bair,mph0,mph1;
    float *rhos=l2rec->rhos;

    if (l2rec->sensorID != MERIS) {
        printf("MPH not supported for this sensor (%s).\n",
                sensorName[l2rec->sensorID]);
        exit(1);
    }

    if (!flags_habs) {
        flags_habs =  (unsigned char*) calloc(l2rec->npix, sizeof( unsigned char));
    }

    wav6  = 620;
    wav7  = 664;
    wav8  = 681;
    wav9  = 709;
    wav10 = 753;
    wav14 = 885;

    ib6  = bindex_get(wav6);
    ib7  = bindex_get(wav7);
    ib8  = bindex_get(wav8);
    ib9  = bindex_get(wav9);
    ib10 = bindex_get(wav10);
    ib14 = bindex_get(wav14);

    if (ib6 < 0 || ib7< 0 || ib8 < 0 || ib9 < 0 || ib10 < 0 || ib14 < 0) {
        printf("MPH_stumpf: incompatible sensor wavelengths for this algorithm\n");
        exit(1);
    }

    for (ip=0; ip<l2rec->npix; ip++) {
        flags_habs[ip] = 0;

        ipb = l2rec->nbands*ip;

        if(l2rec->Rrs[ipb+ib6] <= 0.0 ) {
//            if(l2rec->Rrs[ipb+ib6] <= 0.0 || l2rec->Rrs[ipb+ib7] <= 0.0 || l2rec->Rrs[ipb+ib8] <= 0.0
//            || l2rec->Rrs[ipb+ib9] <= 0.0) {
            chl_mph[ip] = BAD_FLT;
            l2rec->flags[ip] |= PRODFAIL;
        } else {
            if (rhos[ipb+ib8] > rhos[ipb+ib9]) {
                wavmax0 = wav8;
                Rmax0   = rhos[ipb+ib8];
            }else{
                wavmax0 = wav9;
                Rmax0   = rhos[ipb+ib9];
            }
            if (Rmax0 > rhos[ipb+ib10]) {
                wavmax1 = wavmax0;
                Rmax1   = Rmax0;
            }else{
                wavmax1 = wav10;
                Rmax1   = rhos[ipb+ib10];
            }

            //sun-induced phycocyanin absorption fluorescence
            sipf = rhos[ipb+ib7] - rhos[ipb+ib6] - ( rhos[ipb+ib8]  - rhos[ipb+ib6]) * (664 - 619)/(681 - 619);
            //sun induced chlorophyll fluorescence
            sicf = rhos[ipb+ib8] - rhos[ipb+ib7] - ( rhos[ipb+ib9]  - rhos[ipb+ib7]) * (681 - 664)/(709 - 664);
            //normalised difference vegetation index
            ndvi = (rhos[ipb+ib14] - rhos[ipb+ib7])/(rhos[ipb+ib14] + rhos[ipb+ib7]);
            //backscatter and absorption induced reflectance
            bair = rhos[ipb+ib9] - rhos[ipb+ib7] - ( rhos[ipb+ib14] - rhos[ipb+ib7]) * (709 - 664)/(885 - 664);
            mph0 = Rmax0 - rhos[ipb+ib7] - ( rhos[ipb+ib14] - rhos[ipb+ib7]) * (wavmax0 - 664)/(885 - 664);
            mph1 = Rmax1 - rhos[ipb+ib7] - ( rhos[ipb+ib14] - rhos[ipb+ib7]) * (wavmax1 - 664)/(885 - 664);

            if (wavmax1 != wav10) {

                if (sicf >= 0 || sipf <= 0 || bair <= 0.002) {
                    chl_mph[ip] = 5.24e9*pow(mph0,4) - 1.95e8*pow(mph0,3) + 2.46e6*pow(mph0,2) + 3.02e3*mph0 + 1.97;
                } else {
                    chl_mph[ip] = 22.44 * exp(35.79 * mph1);
                }
            } else {
                if (mph1 >= 0.02 || ndvi >= 0.2) {

                    if (sicf < 0 && sipf > 0) {
                            chl_mph[ip] = 22.44 * exp(35.79 * mph1);

                    } else {
                        chl_mph[ip] = BAD_FLT;
                    }
                } else {
                        chl_mph[ip] = 5.24e9*pow(mph0,4) - 1.95e8*pow(mph0,3) + 2.46e6*pow(mph0,2) + 3.02e3*mph0 + 1.97;
                }
            }
        }
        if (chl_mph[ip] < 0.) chl_mph[ip] = 0.;
    }

}

void get_habs_mph_flags(l2str *l2rec,  l2prodstr *p, float flags[])
{
    //MPH - Maximum Peak Height
    // from "Improved algorithm for routine monitoring of cyanobacteria and
    // eutrophication in inland and near-coastal waters"
    // Matthews and Odermatt, Remote Sensing of Environment (doi:10.1016/j.rse.2014.10.010)
    //

    float wav6,wav7,wav8,wav9,wav10,wav14;
    int ip,ipb,ib6,ib7,ib8,ib9,ib10,ib14;
    float Rmax0,Rmax1,wavmax0,wavmax1,ndvi;
    float sipf,sicf,bair,mph1,chl_mph;
    float *rhos=l2rec->rhos;
    static float thresh=350;

    if (l2rec->sensorID != MERIS) {
        printf("MPH not supported for this sensor (%s).\n",
                sensorName[l2rec->sensorID]);
        exit(1);
    }

    if (!flags_habs) {
       flags_habs   = (unsigned char*) calloc(l2rec->npix,sizeof(unsigned char));
    }

    wav6  = 620;
    wav7  = 664;
    wav8  = 681;
    wav9  = 709;
    wav10 = 753;
    wav14 = 885;

    ib6  = bindex_get(wav6);
    ib7  = bindex_get(wav7);
    ib8  = bindex_get(wav8);
    ib9  = bindex_get(wav9);
    ib10 = bindex_get(wav10);
    ib14 = bindex_get(wav14);

    if (ib6 < 0 || ib7< 0 || ib8 < 0 || ib9 < 0 || ib10 < 0 || ib14 < 0) {
        printf("MPH_stumpf: incompatible sensor wavelengths for this algorithm\n");
        exit(1);
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        flags_habs[ip] = 0;
        ipb = l2rec->nbands*ip;

   //     if(l2rec->rhos[ipb+ib6] <= 0.0 ) {
            if(l2rec->rhos[ipb+ib6] <= 0.0 || 
               l2rec->rhos[ipb+ib7] <= 0.0 || 
               l2rec->rhos[ipb+ib8] <= 0.0 || 
               l2rec->rhos[ipb+ib9] <= 0.0 || 
               (l2rec->flags[ip] & LAND) != 0 || 
               (l2rec->flags[ip] & NAVFAIL) != 0 ) {
            l2rec->flags[ip] |= PRODFAIL;
            flags[ip] = BAD_FLT;
        } else {
            if (rhos[ipb+ib8] > rhos[ipb+ib9]) {
                wavmax0 = wav8;
                Rmax0   = rhos[ipb+ib8];
            }else{
                wavmax0 = wav9;
                Rmax0   = rhos[ipb+ib9];
            }
            if (Rmax0 > rhos[ipb+ib10]) {
                wavmax1 = wavmax0;
                Rmax1   = Rmax0;
            }else{
                wavmax1 = wav10;
                Rmax1   = rhos[ipb+ib10];
            }

            //sun-induced phycocyanin absorption fluorescence
            sipf = rhos[ipb+ib7] - rhos[ipb+ib6] - ( rhos[ipb+ib8]  - rhos[ipb+ib6]) * (664 - 619)/(681 - 619);
            //sun induced chlorophyll fluorescence
            sicf = rhos[ipb+ib8] - rhos[ipb+ib7] - ( rhos[ipb+ib9]  - rhos[ipb+ib7]) * (681 - 664)/(709 - 664);
            //normalised difference vegetation index
            ndvi = (rhos[ipb+ib14] - rhos[ipb+ib7])/(rhos[ipb+ib14] + rhos[ipb+ib7]);
            //backscatter and absorption induced reflectance
            bair = rhos[ipb+ib9] - rhos[ipb+ib7] - ( rhos[ipb+ib14] - rhos[ipb+ib7]) * (709 - 664)/(885 - 664);
            //mph0 = Rmax0 - rhos[ipb+ib7] - ( rhos[ipb+ib14] - rhos[ipb+ib7]) * (wavmax0 - 664)/(885 - 664);
            mph1 = Rmax1 - rhos[ipb+ib7] - ( rhos[ipb+ib14] - rhos[ipb+ib7]) * (wavmax1 - 664)/(885 - 664);

            if (wavmax1 != wav10) {

                if (sicf < 0 && sipf > 0 && bair > 0.002) {
                      flags_habs[ip] |= MPH_CYANO;
                       chl_mph = 22.44 * exp(35.79 * mph1);
                       if (chl_mph > thresh)
                           flags_habs[ip] |= MPH_FLOAT;
                }
            } else {
                if (mph1 >= 0.02 || ndvi >= 0.2) {
                        flags_habs[ip] |= MPH_FLOAT;

                    if (sicf < 0 && sipf > 0) {
                            flags_habs[ip] |= MPH_CYANO;
                            chl_mph = 22.44 * exp(35.79 * mph1);
                            if (chl_mph > thresh)
                                flags_habs[ip] |= MPH_FLOAT;

                    }
                } else {
                        flags_habs[ip] |= MPH_ADJ;
                }
            }

            flags[ip] = flags_habs[ip];
        }

    }

}

void get_habs_cldmask(l2str *l2rec, float cld[])
{
    int ip;

    switch (l2rec->sensorID) {
        case MERIS:
//            for (ip=0;ip < l2rec->npix; ip++)
//               cld[ip] = (float) get_cldmask(l2rec->l1rec,ip);
            get_habs_cldmask_meris(l2rec,  cld);
        break;
        case HMODISA:
        case HMODIST:
            get_habs_cldmask_modis(l2rec,  cld);
        break;
    default:
        printf("HABS cldmsk not supported for this sensor (%s).\n",
                sensorName[l2rec->sensorID]);
        exit(1);
   }
}
void get_habs_cldmask_meris(l2str *l2rec,  float cld[])
{
// Cloud Masking for MERIS

    int ib443,ib490,ib620,ib665,ib681,ib709,ib754,ib885,firstCall=1;
    int ip,ipb;
    float  *rhos = l2rec->rhos, *cloud_albedo = l2rec->cloud_albedo;
    float ftemp,cldtmp;

    if (l2rec->sensorID != MERIS) {
        printf("HABS cloud mask not supported for this sensor (%s).\n",
                sensorName[l2rec->sensorID]);
        exit(1);
    }

    if (!flags_habs) {
       flags_habs   = (unsigned char*) calloc(l2rec->npix,sizeof(unsigned char));
    }

    ib443 = bindex_get(443);
    ib490 = bindex_get(490);
    ib620 = bindex_get(620);
    ib665 = bindex_get(665);
    ib681 = bindex_get(681);
    ib709 = bindex_get(709);
    ib754 = bindex_get(754);
    ib885 = bindex_get(885);

    if (ib443 < 0 || ib490 < 0 || ib620< 0 || ib665 < 0 || ib681 < 0 || ib709 < 0 || ib754 < 0 || ib885 < 0) {
        printf("get_habs_cldmask: incompatible sensor wavelengths for this algorithm\n");
        exit(1);
    }

    for (ip=0; ip<l2rec->npix; ip++) {
        flags_habs[ip] = 0;
        ipb = l2rec->nbands*ip;

        if(rhos[ipb+ib443] >= 0.0 && rhos[ipb+ib620] >= 0.0 && rhos[ipb+ib665] >= 0.0 && rhos[ipb+ib681] >= 0.0 && rhos[ipb+ib709] >= 0.0 && rhos[ipb+ib754] >= 0.0) {
//            %{numreal} = (%{rho620} + %{rho665} + %{rho681}) -3*%{rho443} - (%{rho754}-%{rho443})/(754-443) * (620+665+681- 3*443)
            ftemp = rhos[ipb+ib620] + rhos[ipb+ib665] + rhos[ipb+ib681] - 3*rhos[ipb+ib443] - \
                    (rhos[ipb+ib754] - rhos[ipb+ib443])/(754-443)*(620+665+681 - 3*443);
            cldtmp = cloud_albedo[ip] - 3* ftemp;
            flags_habs[ip] = 0;
            //threshold cld > 0.08
            if (cldtmp > 0.08) {
                flags_habs[ip] = HABS_CLOUD;
            }
            //to deal with scum look at relative of NIR and blue for lower albedos
            if ( (rhos[ipb+ib754]+rhos[ipb+ib709]) > (rhos[ipb+ib443]+rhos[ipb+ib490]) && cloud_albedo[ip] < 0.1) flags_habs[ip] = 0;
            if ( ((rhos[ipb+ib754]+rhos[ipb+ib709]) - (rhos[ipb+ib665]+rhos[ipb+ib490])) > 0.01  && cldtmp < 0.15) flags_habs[ip] = 0;
            if ( (((rhos[ipb+ib754]+rhos[ipb+ib709]) - (rhos[ipb+ib665]+rhos[ipb+ib490])))/cldtmp > 0.1 ) flags_habs[ip] = 0;
            if ( (rhos[ipb+ib665] > 0.1) && (cloud_albedo[ip] > 0.15)) {
                flags_habs[ip] = HABS_CLOUD;
            }

            if (rhos[ipb+ib709 >= 0 && rhos[ipb+ib885] >= 0]) {
                if(rhos[ipb+ib885] > rhos[ipb+ib620] && rhos[ipb+ib885] > rhos[ipb+ib709] && rhos[ipb+ib885] > rhos[ipb+ib754] && rhos[ipb+ib885] > 0.01)
                    flags_habs[ip] |= HABS_NONWTR;
            }
            cld[ip] = flags_habs[ip];
        } else {
            cld[ip] = BAD_FLT;
        }
//        if (cld[ip] < 0  && cld[ip] > BAD_FLT) {
//            printf("RJH: flags_habs <0 %d\n",ip);
//        }
    }


}

void get_habs_cldmask_modis(l2str *l2rec, float cld[])
{
// Cloud Masking for MERIS

    int ib469,ib555,ib645,ib667,ib859,ib1240,ib2130,firstCall=1;
    int ip,ipb;
    float *rhos = l2rec->rhos, *Rrs = l2rec->Rrs, *cloud_albedo = l2rec->cloud_albedo;
    float ftemp,ftemp2,ftemp3,cldtmp;
    float cloudthr=0.027;
    if (l2rec->sensorID != HMODISA && l2rec->sensorID != HMODIST) {
        printf("HABS cloud mask not supported for this sensor (%s).\n",
                sensorName[l2rec->sensorID]);
        exit(1);
    }

    if (!flags_habs) {
       flags_habs   = (unsigned char*) calloc(l2rec->npix,sizeof(unsigned char));
    }

    ib469 = bindex_get(469);
    ib555 = bindex_get(555);
    ib645 = bindex_get(645);
    ib667 = bindex_get(667);
    ib859 = bindex_get(859);
    ib1240 = bindex_get(1240);
    ib2130 = bindex_get(2130);

    if (ib469 < 0 || ib555 < 0 || ib645 < 0 || ib667 < 0 || ib859 < 0 || ib1240 < 0 || ib2130 < 0 ) {
        printf("get_habs_cldmask: incompatible sensor wavelengths for this algorithm\n");
        exit(1);
    }

    for (ip=0; ip<l2rec->npix; ip++) {
        flags_habs[ip] = 0;
        ipb = l2rec->nbands*ip;

        if (l2rec->chl[ip] < 0)
            ftemp = Rrs[ipb+ib667];
        else
            ftemp = Rrs[ipb+ib667]*(0.45 + l2rec->chl[ip]*0.005)/4.3;
//      first correct for turbid water
        if(Rrs[ipb+ib667] < 0.0) ftemp = 0.0;
        ftemp2 = cloud_albedo[ip] - ftemp;

        if (ftemp2 > 0.027) flags_habs[ip] = HABS_CLOUD;
//        if (cloud_albedo[ip] > 0.027 && ftemp2 <= 0.027)
//           printf("RJH: %d %d ftemp=%f rhos667=%f cloud_albedo=%f Rrs=%f Chl=%f\n",ip,l2rec->iscan,ftemp,l2rec->rhos[ipb+ib667],cloud_albedo[ip],Rrs[ipb+ib667],l2rec->chl[ip]);

//        non-water check  1240 is bright relative to 859 and the combination is bright
//        this may hit glint by accident, need to be checked.

        if (rhos[ipb+ib1240]/rhos[ipb+ib859] > 0.5 && (rhos[ipb+ib1240] + rhos[ipb+ib2130]) > 0.10) flags_habs[ip] = HABS_CLOUD;

//        now try to correct for glint
//        region check was thrown out {IF (region = "OM") cloudthr = 0.04} rjh 11/2/2015

        ftemp  = rhos[ipb+ib645] - rhos[ipb+ib555] + (rhos[ipb+ib555] - rhos[ipb+ib859])*(645.0 - 555.0)/(859.0 - 555.0);
        ftemp2 = cloud_albedo[ip] + ftemp;
        if (ftemp2 < cloudthr) flags_habs[ip] = 0;
        if (rhos[ipb+ib859]/rhos[ipb+ib1240] > 4.0) flags_habs[ip] = 0;

//     scum areas

        if ((rhos[ipb+ib859] - rhos[ipb+ib469]) > 0.01 &&  cloud_albedo[ip]< 0.30) flags_habs[ip] = 0;
        if ((rhos[ipb+ib859] - rhos[ipb+ib645]) > 0.01 &&  cloud_albedo[ip]< 0.15) flags_habs[ip] = 0;
        if (rhos[ipb+ib1240] < 0.2)
             ftemp2 = ftemp2 - (rhos[ipb+ib859] - rhos[ipb+ib1240])*fabs(rhos[ipb+ib859] - rhos[ipb+ib1240])/cloudthr;

        ftemp3 = ftemp2;
        if (ftemp2 < cloudthr*2) {
            if ( (rhos[ipb+ib555] - rhos[ipb+ib1240]) > (rhos[ipb+ib469] - rhos[ipb+ib1240]) ) {
                ftemp3 = ftemp2 - (rhos[ipb+ib555] - rhos[ipb+ib1240]);
            }else{
                ftemp3 = ftemp2 - (rhos[ipb+ib469] - rhos[ipb+ib1240]);
            }
        }
//        if ( ftemp3*1000.0 >= 0.5)
//            flags_habs[ip] = 1.0;
//        else
//            flags_habs[ip] = 0.0;

        if ( ftemp3 < cloudthr ) flags_habs[ip] = 0.0;

        if(rhos[ipb+ib555] >= 0 && rhos[ipb+ib1240] >= 0.0 && rhos[ipb+ib1240] > rhos[ipb+ib555])
             flags_habs[ip] |= HABS_NONWTR;

        cld[ip] = flags_habs[ip];

    }


}

char get_cldmask(l1str *l1rec,int32_t ip)
{
    //function for cloud mask by pixel
    switch (l1rec->sensorID) {
        case MERIS:
            return(get_cloudmask_meris(l1rec, ip));
        break;
        case HMODISA:
        case HMODIST:
           return(get_cloudmask_modis(l1rec, ip));
        break;
    default:
        printf("HABS cldmsk not supported for this sensor (%s).\n",
                sensorName[l1rec->sensorID]);
        exit(1);
   }
    return(0);
}

char get_cloudmask_meris(l1str *l1rec, int32_t ip)
{
// Cloud Masking for MERIS

    static int ib443,ib490,ib620,ib665,ib681,ib709,ib754,ib885,firstCall=1;
    int ipb;
    float *rhos = l1rec->rhos, *cloud_albedo = l1rec->cloud_albedo;
    float ftemp,cldtmp;
    char flagcld;

    if (firstCall == 1) {
        ib443 = bindex_get(443);
        ib490 = bindex_get(490);
        ib620 = bindex_get(620);
        ib665 = bindex_get(665);
        ib681 = bindex_get(681);
        ib709 = bindex_get(709);
        ib754 = bindex_get(754);
        ib885 = bindex_get(885);

        if (ib443 < 0 || ib490 < 0 || ib620< 0 || ib665 < 0 || ib681 < 0 || ib709 < 0 || ib754 < 0 || ib885 < 0) {
            printf("get_habs_cldmask: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
        firstCall = 0;
    }
        flagcld = 0;
        ipb = l1rec->nbands*ip;

        if(rhos[ipb+ib443] >= 0.0 && rhos[ipb+ib620] >= 0.0 && rhos[ipb+ib665] >= 0.0 && rhos[ipb+ib681] >= 0.0 && rhos[ipb+ib709] >= 0.0 && rhos[ipb+ib754] >= 0.0) {
            ftemp = rhos[ipb+ib620] + rhos[ipb+ib665] + rhos[ipb+ib681] - 3*rhos[ipb+ib443] - \
                    (rhos[ipb+ib754] - rhos[ipb+ib443])/(754-443)*(620+665+681 - 3*443);
            cldtmp = cloud_albedo[ip] - 3* ftemp;
            flagcld = 0;
            if (cldtmp > 0.08) {
                flagcld = 1;
            }
            if ( (rhos[ipb+ib754]+rhos[ipb+ib709]) > (rhos[ipb+ib443]+rhos[ipb+ib490]) && cloud_albedo[ip] < 0.1) flagcld = 0;
            if ( ((rhos[ipb+ib754]+rhos[ipb+ib709]) - (rhos[ipb+ib665]+rhos[ipb+ib490])) > 0.01  && cldtmp < 0.15) flagcld = 0;
            if ( (((rhos[ipb+ib754]+rhos[ipb+ib709]) - (rhos[ipb+ib665]+rhos[ipb+ib490])))/cldtmp > 0.1 ) flagcld = 0;
            if ( (rhos[ipb+ib665] > 0.1) && (cloud_albedo[ip] > 0.15)) {
                flagcld = 1;
            }

//            if (rhos[ipb+ib709 >= 0 && rhos[ipb+ib885] >= 0]) {
//                if(rhos[ipb+ib885] > rhos[ipb+ib620] && rhos[ipb+ib885] > rhos[ipb+ib709] && rhos[ipb+ib885] > rhos[ipb+ib754] && rhos[ipb+ib885] > 0.01)
//                    flagcld = 0;  // = 2 for NOT Water
//            }
        }
//        if (cld[ip] < 0  && cld[ip] > BAD_FLT) {
//            printf("RJH: flags_habs <0 %d\n",ip);
//        }

        return(flagcld);

}

char get_cloudmask_modis(l1str *l1rec, int32_t ip)
{
// Cloud Masking for MODIS

    int ib469,ib555,ib645,ib667,ib859,ib1240,ib2130,firstCall=1;
    int ipb;
    float *rhos = l1rec->rhos, *cloud_albedo = l1rec->cloud_albedo;
    float ftemp,ftemp2,ftemp3,cldtmp;
    float cloudthr=0.027;
    char flagcld;

    if (firstCall == 1) {
        ib469 = bindex_get(469);
        ib555 = bindex_get(555);
        ib645 = bindex_get(645);
        ib667 = bindex_get(667);
        ib859 = bindex_get(859);
        ib1240 = bindex_get(1240);
        ib2130 = bindex_get(2130);
        if (ib469 < 0 || ib555 < 0 || ib645 < 0 || ib667 < 0 || ib859 < 0 || ib1240 < 0 || ib2130 < 0 ) {
            printf("get_habs_cldmask: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
        firstCall = 0;
    }

        ipb = l1rec->nbands*ip;
        flagcld = 0;
        ftemp = 0; //rhos[ipb+ib667];
//      first correct for turbid water

        if(rhos[ipb+ib667] < 0.0) ftemp = 0.0;
        ftemp2 = cloud_albedo[ip] - ftemp;

        if (ftemp2 > 0.027) flagcld = 1;

//        non-water check  1240 is bright relative to 859 and the combination is bright
//        this may hit glint by accident, need to be checked.

        if (rhos[ipb+ib1240]/rhos[ipb+ib859] > 0.5 && (rhos[ipb+ib1240] + rhos[ipb+ib2130]) > 0.10) flagcld = 1;


//        now try to correct for glint
//        region check was thrown out {IF (region = "OM") cloudthr = 0.04} rjh 11/2/2015

        ftemp  = rhos[ipb+ib645] - rhos[ipb+ib555] + (rhos[ipb+ib555] - rhos[ipb+ib859])*(645.0 - 555.0)/(859.0 - 555.0);
        ftemp2 = cloud_albedo[ip] + ftemp;
        if (ftemp2 < cloudthr) flagcld = 0;
        if (rhos[ipb+ib859]/rhos[ipb+ib1240] > 4.0) flagcld = 0;

//     scum areas

        if ((rhos[ipb+ib859] - rhos[ipb+ib469]) > 0.01 &&  cloud_albedo[ip]< 0.30) flagcld = 0;
        if ((rhos[ipb+ib859] - rhos[ipb+ib645]) > 0.01 &&  cloud_albedo[ip]< 0.15) flagcld = 0;
        if (rhos[ipb+ib1240] < 0.2)
             ftemp2 = ftemp2 - (rhos[ipb+ib859] - rhos[ipb+ib1240])*fabs(rhos[ipb+ib859] - rhos[ipb+ib1240])/cloudthr;

        ftemp3 = ftemp2;
        if (ftemp2 < cloudthr*2) {
            if ( (rhos[ipb+ib555] - rhos[ipb+ib1240]) > (rhos[ipb+ib469] - rhos[ipb+ib1240]) ) {
                ftemp3 = ftemp2 - (rhos[ipb+ib555] - rhos[ipb+ib1240]);
            }else{
                ftemp3 = ftemp2 - (rhos[ipb+ib469] - rhos[ipb+ib1240]);
            }
        }
//        if ( ftemp3*1000.0 >= 0.5)
//            flagcld = 1;
//        else
//            flagcld = 0;

        if ( ftemp3 < cloudthr ) flagcld = 0;

        return(flagcld);

}
