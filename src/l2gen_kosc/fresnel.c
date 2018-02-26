/*
 *  theta     Incidence angle (radians)
 *  rel_index Relative index of refraction (n2/n1)
 *  rho       Reflectance
 */

#include <math.h>

float fresnel(float theta,float rel_index)
{
    float rho;
    float alpha;

    if (theta < 0.00001)
        rho = 0.0204078;
    else {
        alpha = asin(sin(theta)/rel_index);
        rho   = ( pow(sin(theta-alpha)/sin(theta+alpha),2.0)
                + pow(tan(theta-alpha)/tan(theta+alpha),2.0))/2.0;
    }

    return(rho);
}


/*
###################################################################
#
# !F90
#
# !Description:
#   The fresnel subroutine calculates the fresnel reflectance given the 
#   relative refractive index between two medium and the incident angle
#     
# !Subroutines and Functions:
#   subroutine FRESNEL (X1,REF,X3)
#   !Input Parameters:
#       X1(real) - Incidence angle in  radians
#       REF(real) - Relative index of refraction (n2/n1); water/air.
#   !Output Parameters:
#       X3(real) - Fresnel Reflectance
#
# !Revision History:
#
# $Id: fresnel.rat,v 1.5 1996/05/08 15:38:04 kay Exp $
#
# $Log: fresnel.rat,v $
# Revision 1.5  1996/05/08  15:38:04  kay
#  update prolog add bang
#
# Revision 1.4  1996/03/11  14:24:29  kay
# update prologs Kay
#
# Revision 1.3  1995/12/04  14:41:11  kay
# add prolog
#
# Revision 1.2  1995/04/27  20:43:19  jim
# Add CVS headers to some files.  Update copyright notices.
#
# Revision 1.1  1994/06/06  19:26:13  jim
# Add module to support new tests.
#
# !Team-Unique Header:
#
# Copyright 1988-1996 by Rosenstiel School of Marine and Atmospheric Science,
# University of Miami, Miami, Florida.
#
#                       All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for non-commercial purposes and without fee is hereby granted,
# provided that the above copyright notice appear in all copies and that both
# that copyright notice and this permission notice appear in supporting
# documentation, and that the names of University of Miami and/or RSMAS not be
# used in advertising or publicity pertaining to distribution of the software
# without specific, written prior permission. 
#
# UNIVERSITY OF MIAMI DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT
# SHALL UNIVERSITY OF MIAMI BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL# DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,# WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING# OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 
#
# !References and Credits:
#
#       Written by:
#       University of Miami
#       Rosensteil School for Marine and Atmospheric Science
#       Division of Meteorology and Physical Oceanography
#       4600 Rickenbacker Cswy
#       Miami,Fl
#
#       contact: SWalsh@rsmas.miami.edu
#
# !Design Notes:
#
# !End############################################################
include "miami"

subroutine FRESNEL (X1,REF,X3)

#  X1  Incidence angle (radians)
#  REF Relative index of refraction (n2/n1)
#  X3  Reflectance
#
#       n1 sin(x1) = n2 sin(x2)
#
#                    tan(x1-x2)**2
#       Refl(par ) = -------------
#                    tan(x1+x2)**2
#
#                    sin(x1-x2)**2
#       Refl(perp) = -------------
#                    sin(x1+x2)**2
#
#  Where:
#       x1  Incidence angle
#       n1  Index refraction of Air
#       x2  Refracted angle
#       n2  Index refraction of Water

single  X1, X2, X3, REF

if (X1 < 0.00001) 
    X3 = .0204078
else {
    X2 = ASIN(SIN(X1)/REF)
    X3 = (SIN(X1-X2)/SIN(X1+X2))**2+(TAN(X1-X2)/TAN(X1+X2))**2
    X3 = X3/2.
}

return
end
*/
