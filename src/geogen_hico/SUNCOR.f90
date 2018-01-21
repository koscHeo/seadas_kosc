!!!********************************************************************************
!!!*                                                                              *
!!!*  Name: SUNCOR                                                                *
!!!*  Purpose: computes at some reference time TZ the declination of the sun and  *
!!!*           the hour angle-TZ.                                                 *
!!!*  Algorithm: At any time T, the local solar hour angle is HAZ+T(GMT)-XLONG,   *
!!!*           where XLONG is the longitude measured positive west of Greenwich.  *
!!!*           Both time and angles are in radians. To compute azimuth and        *
!!!*           elevation at latitude XLAT, longitude XLONG, and time T, call      *
!!!*           HAZEL(HAZ+T-XLONG,DEC,AZ,EL,XLAT). This routine was copied         *
!!!*           from W. Mankin at NCAR, Boulder, CO.                               *
!!!*  Globals used: none                                                          *
!!!*  Global output: none                                                         *
!!!*  Return Codes: none                                                          *
!!!*  Special Considerations: none                                                *
!!!*                                                                              *
!!!********************************************************************************

module suncor_mod

private
public::suncor

contains

  pure SUBROUTINE SUNCOR(IDAY,MONTH,IYR,TZ,DEC,HAZ)
    use solcor_mod
    use julian_mod
    implicit none
    integer, intent(in)::iday,month,iyr
    real (kind=8), intent(in)::tz
    real (kind=8), intent(out)::dec,haz
    integer::jd
    real (kind=8)::fjd,ras,gsdt,bzero,pzero,solong
    real (kind=8), parameter::pi=3.14159265358979323846d0, two_pi=2*pi

    JD=JULIAN(IDAY,MONTH,IYR)
    !      FJD=0.5+TZ/6.283185307
    FJD=0.5d0 + TZ/two_pi
    CALL SOLCOR(JD,FJD,RAS,DEC,GSDT,BZERO,PZERO,SOLONG)
    HAZ=GSDT-RAS-TZ

    RETURN
  END subroutine suncor
end module suncor_mod

