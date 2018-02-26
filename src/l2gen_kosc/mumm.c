/* =================================================================== */
/* MUMM module - turbid water correction for Gordon & Wang atmosphere  */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Ruddick, K., F.Ovidio & M.Rijkeboer (2000). Atmospheric correction  */
/* of SeaWiFS imagery for turbid coastal and inland waters,            */
/* Applied Optics, 39(6), pp897-912.                                   */
/*                                                                     */
/* Written By:                                                         */
/*                                                                     */
/* B. Franz, NASA/OBPG, 03 November 2006 (based on implementation from */
/* http://www.mumm.ac.be/OceanColour/Products/Software/index.php)      */
/*                                                                     */
/* =================================================================== */

#include "l12_proto.h"

/* ------------------------------------------------------------------- */
/* get_rho_mumm(): compute quasi-surface reflectance preferred by MUMM.*/
/* ------------------------------------------------------------------- */

void get_rho_mumm(l2str *l2rec, int32_t ipix, int32_t iw, float *rhom)
{
    int32_t ip, ip1, ip2, ipb;
    float Ltemp;

    if (ipix < 0) {
        ip1 = 0;
        ip2 = l2rec->npix-1;
    } else {
        ip1 = ipix;
        ip2 = ipix;
    }

    for (ip=ip1; ip <= ip2; ip++) {

      if (l2rec->mask[ip]) 
        *rhom++ = 0.0;

      else {

        ipb = ip*l2rec->nbands + iw;

 	Ltemp = ( l2rec->Lt[ipb]
                 / l2rec->tg_sol[ipb]/l2rec->tg_sen[ipb]
                 / l2rec->polcor[ipb]
                 - l2rec->tLf[ipb]
                 - l2rec->Lr[ipb]
		 ) / l2rec->t_o2[ipb]
    	           - l2rec->TLg[ipb];

        *rhom++ = PI*Ltemp/l2rec->Fo[iw]/cos(l2rec->solz[ip]/RADEG);

      } 
    }
}


/* ------------------------------------------------------------------- */
/* get_rhown_mumm(): compute the normalized water-leaving reflectance  */
/* contribution in the NIR using MUMM algorithm.                       */
/* ------------------------------------------------------------------- */

void get_rhown_mumm(l2str *l2rec, int32_t ip, int32_t nir_s, int32_t nir_l, float rhown[])
{
    float alpha   = l2rec->input->mumm_alpha;
    float gamma   = l2rec->input->mumm_gamma;
    float epsilon = l2rec->input->mumm_epsilon;

    float rhom_s, rhom_l;
    float rhoa_s, rhoa_l;

    get_rho_mumm(l2rec,ip,nir_s,&rhom_s);
    get_rho_mumm(l2rec,ip,nir_l,&rhom_l);

    rhoa_l = MAX(MIN((rhom_l*gamma*alpha - rhom_s)/(gamma*alpha - epsilon),rhom_l),0.0);
    rhoa_s = MIN(epsilon*rhoa_l,rhom_s);

    rhown[nir_s] = rhom_s - rhoa_s;
    rhown[nir_l] = rhom_l - rhoa_l;

    return;
}

