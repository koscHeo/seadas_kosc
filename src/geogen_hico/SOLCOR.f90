!!!********************************************************************************
!!!*  Name:  SOLCOR                                                               *
!!!*  Purpose: A routine for solar geometry calculations --- copied from          *
!!!*                  W. Mankin at NCAR.                                          *
!!!*  Parameters: none                                                            *
!!!*  Algorithm:  This routine was supplied by W. Mankin at NCAR, Boulder, CO     *
!!!*  Globals used:                                                               *
!!!*  Global output:                                                              *
!!!*  Return Codes: none                                                          *
!!!*  Special Considerations: none                                                *
!!!*                                                                              *
!!!********************************************************************************
module solcor_mod

private
public::solcor

contains

  pure SUBROUTINE SOLCOR(JD,FJD,RAS,DECS,GSDT,BZRO,P,SOLONG)

    implicit none

    real (kind=8), intent(in):: fjd
    integer, intent(in)::jd
    real(kind=8), intent(out):: ras,decs,gsdt,bzro,p,solong

    !The below 4 llines of reals are only used internally and could be made 
    !double precision. DOes this gain us anything?
    real (kind=8)::d,g,xlms,obl,ecc,f,xlts,sndc,csra,omega,thetac,xlmm,frot,ddecs
    real (kind=8):: dgsdt
    real (kind=8), parameter::pi=3.14159265358979323846_8,twopi=2*pi
    real (kind=8), parameter::jyr=365.25_8
    integer::iyr,jdr,irot

!!!c      PI=3.141592654
    D=(JD-2415020)+FJD
    IYR=D/jyr
!    G=-.026601523_8+.01720196977_8*D-1.95d-15*D*D -2._8*PI*IYR
    G=-.026601523_8+.01720196977_8*D-1.95d-15*D*D -twopi*IYR
!    XLMS=4.881627938_8+.017202791266_8*D+3.95d-15*D*D-2._8*PI*IYR
    XLMS=4.881627938_8+.017202791266_8*D+3.95d-15*D*D-twopi*IYR
    OBL=.409319747_8-6.2179d-9*D
    ECC=.01675104_8-1.1444d-9*D
    F=D-jyr*IYR
!    GSDT=1.739935476_8+2._8*PI*F/365.25_8+1.342027d-4*D/365.25_8
    dGSDT=1.739935476_8+twopi*F/jyr+1.342027d-4*D/jyr
!    GSDT=GSDT+6.2831853_8*(FJD-0.5_8)
    GSDT=dGSDT+twopi*(FJD-0.5_8)
    XLTS=XLMS+2._8*ECC*SIN(G)+1.25_8*ECC*ECC*SIN(2._8*G)
    SNDC=SIN(XLTS)*SIN(OBL)
    dDECS=ASIN(SNDC)
    decs=ddecs
    CSRA=COS(XLTS)/COS(dDECS)
    RAS=ACOS(CSRA)
!    IF(SIN(XLTS).LT.0.) RAS=2._8*PI-RAS
    IF(SIN(XLTS).LT.0.) RAS=twopi-RAS
    OMEGA=1.297906_8+6.66992d-7*D
    THETAC=XLTS-OMEGA
    BZRO=ASIN(.126199_8*SIN(THETAC))
    P=-ATAN(COS(XLTS)*TAN(OBL))-ATAN(.127216_8*COS(THETAC))
    XLMM=ATAN2(.992005_8*SIN(THETAC),COS(THETAC))
    JDR=JD-2398220
    IROT=(JDR+FJD)/25.38_8
    FROT=(JDR+FJD)/25.38_8-IROT
!    SOLONG=XLMM-2._8*PI*FROT+PI-3.d-4
    SOLONG=XLMM-twopi*FROT+PI-3.d-4
    IF(SOLONG.LT.0.) then
!       SOLONG=SOLONG+2._8*PI
       SOLONG=SOLONG+twopi
!    else IF(SOLONG.GT.(2._8*PI)) then
    else IF(SOLONG.GT.(twopi)) then
!       SOLONG=SOLONG-2._8*PI
       SOLONG=SOLONG-twopi
    endif

    RETURN
  END subroutine solcor

end module solcor_mod

