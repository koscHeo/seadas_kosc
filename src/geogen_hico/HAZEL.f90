!!!********************************************************************************
!!!*                                                                              *
!!!*  Name: HAZEL                                                                 *
!!!*  Purpose: Calculates azimuth and elevation                                   *
!!!*  Parameters: none                                                            *
!!!*  Algorithm: This routine was supplied by William Mankin at NCAR, Boulder, CO *
!!!*  Globals used: None.                                                         *
!!!*  Global output:   None.                                                      *
!!!*  Return Codes: None                                                          *
!!!*  Special Considerations: None.                                               *
!!!*                                                                              *
!!!********************************************************************************
module hazel_mod

private
public::hazel

contains

  pure SUBROUTINE HAZEL (H,D,A,E,XLAT)

    implicit none

    real (kind=8), intent(in)::D
    real (kind=8), intent(in),dimension(:)::H,XLAT
    real (kind=8), intent(out),dimension(size(H))::a,e
    real (kind=8),dimension(size(H)) ::sne,sna,csa
    real (kind=8),parameter::pi=3.14159265358979323846d0

!!!C     H = HOUR ANGLE
!!!C     D = DECLINATION
!!!C     A = AZIMUTH
!!!C     E = ELEVATION
!!!C     XLAT = LATITUDE
!!!C     ALL ANGLES IN RADIANS

!!!c      PI = 3.14159265
    SNE = SIN(D)*SIN(XLAT)+COS(D)*COS(XLAT)*COS(H)
    E=ASIN(SNE)
    SNA = COS(D)*SIN(H)
    CSA=(SIN(XLAT)*COS(H)*COS(D)-SIN(D)*COS(XLAT))
    A=ATAN2(SNA,CSA)+PI
  END subroutine hazel

end module hazel_mod
