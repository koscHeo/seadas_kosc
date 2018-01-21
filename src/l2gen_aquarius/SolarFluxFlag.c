/* based on
function [HFLAG] = Hsun_Flag(x_ECI,y_ECI,z_ECI,xSUN,ySUN,zSUN,Rave)
// 
//         Function Hsun_Flag.m
// INPUT:
// 	  x_ECI,y_ECI,z_ECI - Position vector of satellite in ECI
// 	  xSUN, ySUN, zSUN - Position vector of sun in ECI
//         Rave - Earth Radius at nadir
// OUTPUT: 
//         HFLAG (0, Earth blocakge OR [1, DIRECT and REFLECTED sun possible])
// Functions: 
//         ECI2SRF.m: Matrix for converting ECI to SRF
//
//	  Created by Saji Abraham, 2/17/06
*/

//    Modification history:
// Programmer  Organization  Date      Version  Description of change
// ----------  ------------  ----      -------  ---------------------
// Joel Gales  FutureTech    05/25/11           Fix theta_ind[1] bound check
//

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>


#define PI 3.1415926536

int32_t Hsun_Flag(double orb_pos[3], double sun_vec[3], double *scalt, 
		  int32_t *HFLAG)
{      
  // orb_pos is the Orbital position vector = 
  //                Position vector of satellite in ECI
  // sun_vec is the Position unit vector to the sun in ECI
  // scalt is Spacecraft altitude
  // HFLAG is output, 1: sun above horizon, 0: sun below horizon
  float Rave;		// Earth Radius at nadir
  float scdist;	// Spacecraft distance from center of earth
  float ThetaH, Decr;

  // Conversion of ECI vector into SRF
  scdist = sqrt(orb_pos[0]*orb_pos[0] + 
		orb_pos[1]*orb_pos[1] + 
		orb_pos[2]*orb_pos[2]);
  Rave = scdist - (*scalt);
 
  // Terminator angle relative to sun vector
  ThetaH = acos(Rave/scdist) + PI/2;

  // angle between s/c and sun
  Decr = acos((orb_pos[0]*sun_vec[0] + orb_pos[1]*sun_vec[1] + orb_pos[2]*sun_vec[2])/scdist);

  //    HFLAG: Horizon Flag of sun in ECI'
  //    HFLAG=1: sun appears above the satellite horizon
  //    HFLAG=0: Earth Blockage (sun disappears below the horizon)
  //  printf("%f %f\n", ThetaH, Decr);

  if (ThetaH > Decr) {
    *HFLAG = 1;
  }
  else {
    *HFLAG = 0;
  }
  return 0;
}


/*
 *  SolarFluxFlag.c
 *  
 *
 *  Created by Liang Hong on 12/20/2010.
 *  Copyright 2010 NASA/GSFC/SAIC. All rights reserved.
 *
 *  ref: Solar Flag Algorithm, Revised April 21, 2010
 *
 *	calculates the flag from direct and reflected mean solar flux 
 *
 */
int SolarFluxFlag( uint32_t ibeam, 
		   float RD_Ga4[3][2][260281], 
		   float RD_Ga16[3][2][260281], 
		   float RD_Gm16[3][2][260281], 
		   float SolarFlux, 
		   double orb_pos[3], double *scalt,
		   double Rd[3], double Rr[3], 
		   double Za[3], double Xa[3], double Ya[3], 
		   uint8_t MSFF[4])
{		
  // calc. mean solar flux & solar flare flag per radiometer
  // Rd is the vectors in the direction of direct ray from the sun
  // Rr is the vectors in the direction of reflected ray from the sun
  // Za: direction of antenna boresight ray
  // Xa: vector in the direction of antenna V-port
  // Ya: vector in the direction of antenna H-port
  // MSFF: V direct, V reflect, H direct, H reflect 
	
  float TA4, TA16, TA16m;
  float R[2] = {1, 0.325};	// direct R = 1; reflected R = 0.325
  int ichan, iRay;
  int theta_ind[2], phi_ind[2];
  int32_t HFLAG;

  MSFF[0] = 0;
  MSFF[1] = 0;
  MSFF[2] = 0;
  MSFF[3] = 0;

  // If sun below horizon then set flag values to 0
  Hsun_Flag( &orb_pos[0], &Rd[0], scalt, &HFLAG);
  if ( HFLAG == 0) return 0;

  // calc. Gij[ncyc][6] for V & H pol

  // direct solar
  theta_ind[0] = rint(acos((Za[0]*Rd[0] +  Za[1]*Rd[1] + Za[2]*Rd[2]) /
			   sqrt(Za[0]*Za[0] + Za[1]*Za[1] + Za[2]*Za[2]) /
			   sqrt(Rd[0]*Rd[0] + Rd[1]*Rd[1] + Rd[2]*Rd[2]))
		      * (360/PI));

  if (theta_ind[0] < 0 || theta_ind[0] > 360) {
    printf("theta_ind[0] out of bounds: %d\n", theta_ind[0]);
    return 0;
  }

  phi_ind[0] = rint(atan2((Ya[0]*Rd[0] + Ya[1]*Rd[1] + Ya[2]*Rd[2]) ,
			  (Xa[0]*Rd[0] + Xa[1]*Rd[1] + Xa[2]*Rd[2]))
		    * (360/PI) + 360);

  if (phi_ind[0] < 0 || phi_ind[0] > 720) {
    printf("phi_ind[0] out of bounds: %d\n", phi_ind[0]);
    return 0;
  }

  //printf("Rd:%1d %9.4lf %9.4lf %9.4lf\n", ibeam, Rd[0], Rd[1], Rd[2]);
  //printf("Rr:%1d %9.4lf %9.4lf %9.4lf\n", ibeam, Rr[0], Rr[1], Rr[2]);

  //printf("Xa:%1d %9.4lf %9.4lf %9.4lf\n", ibeam, Xa[0], Xa[1], Xa[2]);
  //printf("Ya:%1d %9.4lf %9.4lf %9.4lf\n", ibeam, Ya[0], Ya[1], Ya[2]);
  //printf("Za:%1d %9.4lf %9.4lf %9.4lf\n", ibeam, Za[0], Za[1], Za[2]);

  //printf("ph:%1d %d\n", ibeam, phi_ind[0]);

  // reflected solar
  theta_ind[1] = rint(acos((Za[0]*Rr[0] + Za[1]*Rr[1] + Za[2]*Rr[2]) /
			   sqrt(Za[0]*Za[0] + Za[1]*Za[1] + Za[2]*Za[2]) /
			   sqrt(Rr[0]*Rr[0] + Rr[1]*Rr[1] + Rr[2]*Rr[2]))
		      * (360/PI));

  if (theta_ind[1] < 0 || theta_ind[1] > 360) {
    printf("theta_ind[1] out of bounds: %d\n", theta_ind[1]);
    return 0;
  }

  phi_ind[1] = rint(atan2((Ya[0]*Rr[0] + Ya[1]*Rr[1] + Ya[2]*Rr[2]) ,
			  (Xa[0]*Rr[0] + Xa[1]*Rr[1] + Xa[2]*Rr[2]))
		    * (360/PI) + 360);
		
  if (phi_ind[1] < 0 || phi_ind[1] > 720) {
    printf("phi_ind[1] out of bounds: %d\n", phi_ind[1]);
    return 0;
  }

  phi_ind[0] = (phi_ind[0] + 360) % 720;
  phi_ind[1] = (phi_ind[1] + 360) % 720;

  // RD[260281 = 721*361]
  for (ichan=0; ichan<2; ichan++) {
    for (iRay=0; iRay<2; iRay++) {
	  
      TA4 = 0.013 * RD_Ga4[ibeam][ichan][theta_ind[iRay]*721+phi_ind[iRay]] 
	* SolarFlux * R[iRay];
      TA16 = 0.013 * RD_Ga16[ibeam][ichan][theta_ind[iRay]*721+phi_ind[iRay]] 
	* SolarFlux * R[iRay];
      TA16m = 0.013 * RD_Gm16[ibeam][ichan][theta_ind[iRay]*721+phi_ind[iRay]] 
	* SolarFlux * R[iRay];
	  
      if (TA16m<0.05 || HFLAG == 0)
	MSFF[ichan*2+iRay] = 0;
      else if ((TA4>0.05 && TA16<0.05 && TA16m<0.05) || 
	       (TA4<0.05 && TA16>0.05) ||
	       (TA4<0.05 && TA16<0.05 && TA16m>0.05))
	MSFF[ichan*2+iRay] = 1;
      else if (TA4>=0.05 && TA16>=0.05)
	MSFF[ichan*2+iRay] = 2;
    }
  }	
  return 0;
}


