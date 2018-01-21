!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
module auxtools

  implicit none
  private
  public::frac_to_deg_min_sec,dms_to_frac_deg, &
       deg_min_sec_to_frac_min,deg_min_sec_to_frac_sec,&
       date_to_day_in_year,my_minloc

  integer,parameter,dimension(12):: MD=(/0,31,59,90,120,151,181,212        &
       &     ,243,273,304,334/)  ! only used to calculate iday

  interface dms_to_frac_deg
     module procedure deg_min_sec_to_frac_deg, ideg_imin_sec_to_frac_deg
  end interface
  
contains
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure subroutine frac_to_deg_min_sec(fracd,fracm,fracs,deg,min,sec)
    ! A subroutine to get degrees, min, sec from
    ! fractional pieces. The output are all real, although a 
    ! strong case can be made for degrees and minutes to be
    ! integers.
    ! MJM Unknown date
    implicit none
    real, intent(in)::fracd, fracm, fracs
    real, intent(out)::deg, min, sec
    real:: d_hold, m_hold, s_hold,d_all,m_all,s_all

    d_all=fracd+fracm/60.+fracs/3600.
    d_hold=floor(d_all) !thorough
    m_all=(d_all-d_hold)*60. ! remaining number of arc minutes
    m_hold=floor(m_all)
    s_all=(m_all-m_hold)*60. ! remaining number of arc seconds
    s_hold=floor(s_all)      ! integer portion, not used

    deg=d_hold
    min=m_hold
    sec=s_all
    
  end subroutine frac_to_deg_min_sec
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure subroutine deg_min_sec_to_frac_deg(deg,min,sec,frac_deg)
    ! Convert deg-min-sec angles to fractional degrees
    !MJM 2002-12-06
    implicit none
    real,intent(in)::deg,min,sec
    real,intent(out)::frac_deg
    
    frac_deg= deg + (min + sec/60.)/60.

  end subroutine deg_min_sec_to_frac_deg
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure subroutine ideg_imin_sec_to_frac_deg(deg,min,sec,frac_deg)
    ! Convert deg-min-sec angles to fractional degrees
    !MJM 2002-12-06
    implicit none
    integer,intent(in)::deg,min
    real,intent(in)::sec
    real,intent(out)::frac_deg
    
    frac_deg= real(deg) + (real(min) + sec/60.)/60.

  end subroutine ideg_imin_sec_to_frac_deg
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure subroutine deg_min_sec_to_frac_min(deg,min,sec,frac_min)
    ! Convert deg-min-sec angles to fractional minutes
    !MJM 2002-12-06
    implicit none
    real,intent(in)::deg,min,sec
    real,intent(out)::frac_min
    
    frac_min= 60.*deg + min + sec/60.
    
  end subroutine deg_min_sec_to_frac_min
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure subroutine deg_min_sec_to_frac_sec(deg,min,sec,frac_sec)
    ! Convert deg-min-sec angles to fractional seconds
    !MJM 2002-12-06
    implicit none
    real,intent(in)::deg,min,sec
    real,intent(out)::frac_sec
    
    frac_sec= (60.*deg + min)*60. + sec
    
  end subroutine deg_min_sec_to_frac_sec
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure subroutine date_to_day_in_year(iyr,imn,idy,iday)
    implicit none
    integer, intent(in)::iyr,imn,idy
    integer, intent(out)::iday
    !MJM 2002-12-06
!!!C Calculate the number of days that have passed in this year.  Take leap year
!!!C into account. Passed to solar_irr.f90
    IDAY = MD(IMN) + IDY !first, non-leap year
    
    ! I'm a calendar nut, so fix leap year rule. Not that this will be used 'til 2100,
    ! but it's better to get it correct. Old: lpyr=0 indicates leap year
    ! new: don't assign lpyr variable
    !    LPYR = IYR - (4 * (IYR/4)) ! Strictly a Julian rule, will not work 2100, 
    !           2200, etc.
    !    IF((LPYR.EQ.0).AND.(IDAY.GT.59).AND.(IMN.NE.2)) IDAY = IDAY + 1
    if ((modulo(iyr,4).eq.0.and.modulo(iyr,100).ne.0).or.&
         &(modulo(iyr,400).eq.0)) then
       if ((IDAY.GT.59).AND.(IMN.NE.2)) IDAY = IDAY + 1
    end if
  end subroutine date_to_day_in_year
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure function my_minloc(v) result(loc)
    ! MJM 2004 Apr 22
    ! Function to provide a scalar value as the return value for the 
    ! minimum value index of a vector.
    ! Note: Fortran 95 and Fortran 2003 do not need this, I believe.
    ! Drat. Seems to slow things down just a bit compared to performing this in
    ! in each subroutine.
    implicit none
    real, intent(in),dimension(:)::v
    integer::loc
    integer, dimension(1)::hold

    hold=minloc(v)
    loc=hold(1)

  end function my_minloc
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
end module auxtools
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
