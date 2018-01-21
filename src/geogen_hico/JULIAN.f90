!!!*******************************************************************************
!!!*                                                                             *
!!!*  Name: JULIAN                                                               *
!!!*  Purpose: computes the julian day                                           *
!!!*  Parameters: IDAY - day of the year                                         *
!!!*              MONTH - day of the year                                        *
!!!*  Algorithm: This routine was supplied by William Mankin at NCAR, Boulder, CO*
!!!*  Globals used: none                                                         *
!!!*  Global output:  none                                                       *
!!!*  Return Codes: none                                                         *
!!!*  Special Considerations: none                                               *
!!!*                                                                             *
!!!*******************************************************************************

module julian_mod
  
private
public::julian

contains
  
  pure integer FUNCTION JULIAN(IDAY,MONTH,IYR)
    
    implicit none
    integer, intent(in)::iday,month,iyr

!!!C Local variables
    integer::jyr,i1,i2,i3,leap
    integer,dimension(12),parameter::md=(/0,31,59,90,120,151,181,212            &
         &             ,243,273,304,334/) 

!!!C
!!!c     Next statement commented because 1)taken care of in get_input; 2)
!!!c     it is bad form to modify function arguments;3) we don't allow
!!!c     years <100. We are Y2K compliant!!!!!
!!!c      IF(IYR.LT.100) IYR=IYR+1900 ! we don't allow years<100; see GET_INPUT

    JYR=IYR-1600
    I1=JYR/400           !Number of complete 400 yr spans; each have 97 leap years
    I2=(JYR-400*I1)/100  !Number of 100 year spans, taking away the 400 year spans above
    I3=(JYR-400*I1-100*I2)/4 !Number of four year spans, taking away complete 400 and 100 above
    !    JD(AD1600) normal yr  400 yr  100yr  leap yr this century
    JULIAN=2305447 + 365*JYR + 97*I1 + 24*I2 +I3 !Counts IYR if IYR a leap year
    LEAP=0
    IF(MOD(JYR,4).EQ.0) LEAP=1
    IF(MOD(JYR,100).EQ.0) LEAP=0
    IF(MOD(JYR,400).EQ.0) LEAP=1
    
!!!C     LEAP=1 if iyr is a leap year
    JULIAN=JULIAN+MD(MONTH)+IDAY
    IF(MONTH.LE.2) JULIAN=JULIAN-LEAP !Already counted leap day above, so take out if jan/feb
    
    RETURN
  END function julian
  
end module julian_mod
