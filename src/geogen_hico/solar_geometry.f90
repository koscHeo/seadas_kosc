module sol_geo
  implicit none
  private
  public::solar_geometry

contains

  pure function solar_geometry(IYR,Imn,idy,IH,IM,IS,xlat,xlong) result(sol_zen_az)
    ! 2012-10-31 MJM Made able to handle an array of xlat and xlong, where xlat and xlong
    !                are each vectors
    ! 2002-12-11 MJM Turned into array valued function. Now should only 
    !                necessary solar geometry calculations. This is 
    !                called many nlines times when given appropriate info for
    !                each line.
    ! xlat: positive NORTH
    ! xlon: positive WEST 
    ! 2001-05-16 MJM Solar factors only.
    use suncor_mod
    use hazel_mod
!    use errors_mod,only:fatal_error
    implicit none
!    include 'ASTMATH.CMN'
    integer,intent(in)::iyr,imn,idy,ih,im !year, month,date,hour,minute
    real (kind=8),intent(in)::IS ! seconds
    real (kind=8) ,intent(in),dimension(:)::XLAT, XLONG
    real (kind=8),dimension(2,size(xlat))::sol_zen_az

    real (kind=8), parameter:: Pi      =  4.0D0 * DATAN(1.0D0)
    real (kind=8), parameter:: HalfPi  =  0.5D0*pi
    real (kind=8), parameter:: TwoPi   =  2.0D0*pi
    real (kind=8), parameter:: Rad2Deg = 180.0D0/pi
    real (kind=8), parameter:: Deg2Rad = pi/180.0D0

!    character (len=19),parameter::subroutine_name='solar_geometry'
    real (kind=8)::XH, XM, XS, TT,  DEC, HAZ
    real (kind=8),dimension(size(xlat))::solaz,solzni,XLATR,XLONGR, EL

    XH= IH*1.d0
    XM= IM*1.d0
    XS=IS
    !     TT=6.28318*((XH)/24+XM/1440+XS/86400)
    TT=twopi*((XH)/24+XM/1440+XS/86400)
    
    XLATR = XLAT * Deg2Rad
    XLONGR = XLONG * Deg2Rad

    ! at this TIME there is only ONE astronomical solar position
    CALL SUNCOR(IDY,IMN,IYR,TT,DEC,HAZ)

    CALL HAZEL(HAZ+TT-XLONGR,DEC,SOLAZ,EL,XLATR)

!!!C---Note: DEC, SOLAZ,and EL DEC are in radians

!    IF (EL .LE. 0.) THEN
!       print*,IYR,Imn,idy,IH,IM,IS,xlat,xlong
!       print*,'el=',el
!       call fatal_error(file_name,module_name,subroutine_name,&
!            &'ERROR: Sun is below the horizon!!!'//&
!            &' Check input date, time, latitude and longitude.')
!    ENDIF
    
    !      SOLZNI = 90.0/57.2958 - EL
    SOLZNI = HalfPi - EL
!    sol_zen_az=(/solzni,solaz/)
    sol_zen_az(1,:)=solzni ! radians
    sol_zen_az(2,:)=solaz  !radians

  end function solar_geometry
end module sol_geo
