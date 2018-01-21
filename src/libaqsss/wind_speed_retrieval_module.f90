module wind_speed_retrieval_module

use static_data_module

implicit none
private
public :: minimum_search 
public :: fd_wspd
public :: fd_est_error
save

! global variable accessed by chi2 functions 
real(4)                                     :: ysigma0_VV, ysigma0_HH, ydew_V, ydew_H, ywin_0, yphir  
real(4)                                     :: ywt_HH, ywt_H, ywt_W
integer(4)                                  :: yirad

real(4)                                     :: ysst, ysss 
! added in June 2013, as dEW is now SST/SSS dependent

type(static_data_struct)                    :: stdata

contains

subroutine fd_wspd(irad, phir, wspd_0, dew, xsigma0_VV, xsigma0_HH, &
     sst, sss, staticdata, &
     wspd_scat, chisq_scat, iflag_wspd_scat)
! scat HH
! December 27, 2012: substitute the NUMREC BRENT with my Golden Section minimum search
  use static_data_module

  implicit none

  integer(4), intent(in)                  :: irad
  real(4), intent(in)                     :: phir     ! relative wind direction
  real(4), intent(in)                     :: wspd_0   ! first guess wind speed (use NCEP)
  real(4), intent(in)                     :: dew(2)
  real(4), intent(in)                     :: xsigma0_VV, xsigma0_HH
  real(4), intent(in)                     :: sst, sss

  type(static_data_struct), intent(in)    :: staticdata

  real(4), intent(out)                    :: wspd_scat, chisq_scat
  integer(4), intent(out)                 :: iflag_wspd_scat

  integer(4), parameter           :: mmax=4
  real(4), dimension(mmax)        :: wa, wb, wc
  real(4), parameter              :: bw=0.5  
  integer(4)                      :: iw, nmax, min_count
  real(4)                         :: w0, w1, w2, y0, y1, y2

  real(4)                         :: dsigma0, ddew, est_err(3), dw

  real(4), dimension(mmax)        :: x_min, y_min
  integer(4), dimension(mmax)     :: m_search

  real(4), parameter              :: very_large_number = 1.0E15
  real(4)                         :: dmin, diff
  integer(4)                      :: imin, jmin

  real(4), parameter              :: eps=1.0E-2
  real(4), parameter              :: wmin=0.0, wmax=30.0
  integer(4),save  :: istart=1

  if(istart.eq.1) then
     istart=0
     stdata = staticdata
  endif

  ! save global variables to be accessed through chi2_HH
  ysigma0_VV = xsigma0_VV
  ysigma0_HH = xsigma0_HH
        
  yphir      = phir
  yirad      = irad 
        
  ywin_0     = wspd_0

  ysst       = sst
  ysss       = sss

  call fd_est_error (irad, wspd_0,  est_err)

  dsigma0 = est_err(1)
  ddew    = est_err(2)
  dw      = est_err(3)
  ywt_HH  = 1.0/(dsigma0**2) ! channel weighted by expected variance
  ywt_W   = 1.0/(dw**2)      ! weight for background field

  if (loc(dew).ne.0) then
     ydew_V     = dew(1)
     ydew_H     = dew(2)
     ywt_H      = 1.0/(ddew**2) ! rad channel H weighted by expected variance 
  endif


!  if (irad.eq.1) write(*,*) 'find minima'
  min_count=0
  nmax = nint(wmax/bw)
  do iw=1,nmax-1
     w0 = wmin + float(iw-1)*bw
     w1 = wmin + float(iw)*bw  
     w2 = wmin + float(iw+1)*bw   
     if (w0<0.0) w0=0.0
     if (loc(dew).ne.0) then
        y0 = chi2_HHH(w0)
        y1 = chi2_HHH(w1)
        y2 = chi2_HHH(w2)
     else
        y0 = chi2_HH(w0)
        y1 = chi2_HH(w1)
        y2 = chi2_HH(w2)
     endif
     if (y1<=y0 .and. y1<y2) then ! found minimum with brackets
        min_count=min_count+1
        wa(min_count)= w0
        wb(min_count)= w1
        wc(min_count)= w2
     endif
  enddo
  
!  if (irad.eq.1) write(*,*) 'min_count: ', min_count 
  if (min_count==0) then ! function is either monotonously rising or falling
     if (loc(dew).ne.0) then
        y0 = chi2_HHH(wmin)
        y1 = chi2_HHH(wmax)
     else
        y0 = chi2_HH(wmin)
        y1 = chi2_HH(wmax)
     endif
     if (y0<y1) then
        x_min(1)=wmin
        y_min(1)=y0
     else
        x_min(1)=wmax
        y_min(1)=y1
     endif
     iflag_wspd_scat=0 ! o.k.
     wspd_scat=x_min(1)
     chisq_scat=y_min(1) 
     return   
  endif
  
  if (min_count>mmax) then
     write(*,*)' chi2 has more than mmax=4 minima'
     write(*,*) irad
     write(*,*) wspd_0, phir
     write(*,*) min_count
     write(*,*) wb(1:min_count)
     iflag_wspd_scat=1 ! no scat wspd
     wspd_scat=missing_val
     chisq_scat=missing_val
     return
  endif

  !  if (irad.eq.1) write(*,*) 'minimum search'
  ! improve on minima using Golden Section             
  ! if more than one minimum choose the one that is closest to NCEP 
  jmin=0
  dmin=very_large_number
  do imin =1,min_count
        
     if (loc(dew).ne.0) then
        call minimum_search(wa(imin),wb(imin),wc(imin),chi2_HHH,eps,x_min(imin),y_min(imin),m_search(imin))
     else
        call minimum_search(wa(imin),wb(imin),wc(imin),chi2_HH, eps,x_min(imin),y_min(imin),m_search(imin))
     endif

     if (m_search(imin)/=0) cycle            
       
     diff  = abs(x_min(imin)-ywin_0)
     if (diff < dmin) then
        dmin=diff
        jmin=imin
     endif
  enddo
        
  if (jmin==0) then
     iflag_wspd_scat=1 ! no scat wspd
     write(*,*) ' Golden Section method did not find minimum'
     wspd_scat=missing_val
     chisq_scat=missing_val
  else 
     iflag_wspd_scat=0 ! o.k.
     wspd_scat=x_min(jmin)
     chisq_scat=y_min(jmin)
  endif
        
  return
end subroutine fd_wspd


subroutine fd_est_error (irad, wspd,  est_err)
  implicit none

  integer(4), intent(in)              ::  irad
  real(4), intent(in)                 ::  wspd

  real(4), dimension(3), intent(out)  ::  est_err
  ! 1= sigma0 2= emiss 3=NCEP wind speed

  real(4), parameter    :: wbin=1.0, wlow = 0.5
  integer(4)            :: iw1, iw2 
  real(4)               :: ww, brief, w1, w2
  real(4), dimension(7) :: d1, d2, d
  integer(4), parameter :: iwmax=22

  est_err=0.0
    
  ww = wspd 
  if (ww < wlow)  ww=wlow
   
  iw1 = floor((ww-wlow)/wbin) + 1
  if (iw1<1)           iw1 = 1
  if (iw1> iwmax-1)    iw1 = iwmax-1 ! linear interpolation above iwmax
  iw2 = iw1 + 1
        
  w1 = wlow + wbin*(iw1-1)
  w2 = w1   + wbin
  brief = (ww-w1)/(w2-w1)
    
  d1(:) = stdata%estimated_error_array(:,iw1)
  d2(:) = stdata%estimated_error_array(:,iw2)
    
  d(:) = d1(:)*(1.0-brief) + d2(:)*brief ! linear interpolation
    
  if (irad==1) then
     est_err(1) = d(1)
     est_err(2) = d(4)
     est_err(3) = d(7)
  else if (irad==2) then
     est_err(1) = d(2)
     est_err(2) = d(5)
     est_err(3) = d(7)       
  else if (irad==3) then
     est_err(1) = d(3)
     est_err(2) = d(6)
     est_err(3) = d(7)
  endif

  return
end subroutine fd_est_error


real(4) function chi2_HH(w)
  use static_data_module

  !     sos evaluation: HH scatterometer + NCEP background wind speed
  !     fixed phir
  !     variable w
  
  implicit none
  real(4)                            :: w

  integer(4)                         :: irad
  real(4), dimension(0:2,npol_sc)    :: bharm  !1=VV, 2=HH, 3=VH, 4=HV
  real(4), dimension(0:2,npol_sc)    :: dbharm  

  real(4)                            :: f0_HH

  real(4) cosd

  irad=yirad
!  if (irad.eq.1) write(*,*) 'In chi2_HH'

  call fd_TM_scat_harmonics (irad, w, stdata, bharm, dbharm)

  f0_HH = bharm(0,2) + bharm(1,2)*cosd(yphir) + bharm(2,2)*cosd(2.0*yphir)
  chi2_HH = ywt_HH*(ysigma0_HH - f0_HH)**2 + ywt_W*((w - ywin_0)**2)

  return
end function chi2_HH


real(4) function chi2_HHH(w)
  use static_data_module
  !     combined radiometer (H-pol) and scatterometer (HH) + NCEP background
  !     sos evaluation
  !     fixed phir
  !     variable w

  implicit none
  real(4)                            :: w
       
  integer(4)                         :: irad
  real(4), dimension(0:2,2)          :: aharm  !1=VV, 2=HH, 3=VH, 4=HV
  real(4), dimension(0:2,npol_sc)    :: bharm  !1=VV, 2=HH, 3=VH, 4=HV
  real(4), dimension(0:2,2)          :: daharm  
  real(4), dimension(0:2,npol_sc)    :: dbharm  

  real(4)                            :: f0_HH, fdew_H

  real(4), dimension(2)              :: dew_1, dew_r, demiss_sst_wspd

  real(4) cosd
  integer, pointer :: nptr => null()
  real(4), parameter :: teff=290. ! used for wind induced emissivity harmonics

  if (stdata%W_STITCH.ne.-1) then
     call fd_TM_emiss_harmonics (irad, w, nptr, nptr, stdata, aharm, daharm)
     dew_1(1:2) = aharm(0,1:2)/teff
     call fd_TM_emiss_harmonics (irad, w, nptr, nptr, stdata, aharm, daharm)
     dew_r(1:2) = aharm(0,1:2)/teff
     call fd_demiss_sst_wspd (irad, ysst, stdata, dew_1, dew_r, demiss_sst_wspd)    
  else
     demiss_sst_wspd = 0.0
  endif

  irad=yirad
  if (stdata%emiss_sst_sss .eq. 1) then
     call fd_TM_emiss_harmonics (irad, w, ysst, ysss, stdata, aharm, daharm)
  else
     call fd_TM_emiss_harmonics (irad, w, nptr, nptr, stdata, aharm, daharm)
  endif

  call fd_TM_scat_harmonics (irad, w, stdata, bharm, dbharm)
	
  f0_HH = bharm(0,2) + bharm(1,2)*cosd(yphir) + bharm(2,2)*cosd(2.0*yphir)
  fdew_H= aharm(0,2) + aharm(1,2)*cosd(yphir) + aharm(2,2)*cosd(2.0*yphir) + demiss_sst_wspd(2)

  chi2_HHH = ywt_HH*(ysigma0_HH - f0_HH)**2 + ywt_H*(ydew_H - fdew_H)**2  + ywt_W*((w - ywin_0)**2)

return
end function chi2_HHH


subroutine minimum_search(a,ystart,b,func,eps,  x_min,y_min,jflag)
  !
  !     written by T. Meissner (Remote Sensing Systems)
  !
  !     NAME:
  !     minimum_search    
  !
  !     USAGE:
  !     CALL minimum_search (a,ystart,b,func,eps,  x_min,y_min,jflag)
  !
  !       
  !         DESCRIPTION:
  !     finds minimum x_min 
  !     of 1 dimensional function func(x) in interval [a,b]
  !         using successive interval partition by Golden Section
  !     needs f(ystart) < f(a) and f(ystart) < f(b), 
  !     i.e. one needs to be certain that a nad b corner the minimum  
  !
  !
  !    INPUT:
  !    NAME       DESCRIPTION                 TYPE          LENGTH                RANGE
  ! ----------------------------------------------------------------------------------------------
  !
  !    a         left interval boundary       REAL(4)  SCALAR
  !    b         right interval boundary      REAL(4)  SCALAR
  !    ystart    first guess for minimum      REAL(4)  SCALAR                [a,b]  
  !    func      external function            REAL(4)  EXTERNAL FUNCTION 
  !    eps       accuracy goal                REAL(4)  SCALAR 
  !
  !
  !    OUPTUT:
  !    NAME       DESCRIPTION                 TYPE                                         RANGE
  ! ----------------------------------------------------------------------------------------------
  !
  !    x_min      minimum                     REAL(4)                        [a,b]
  !    y_min      func(x_min)                 REAL(4)
  !    jflag          error flag                  INTEGER(4)  
  !    = 0        normal return, accuracy goal achieved within nmax iterations
  !    =-1        f(ystart) >= f(a) or f(ystart) >= f(b) 
  !    = 1        accuracy goal NOT achieved within nmax iterations
  !               in this case the routine returns x_min and y_min from the last step
  !       ACCURACY GOAL: 
  !       |Y1 - Y4| < EPS*|Y1+Y4|/2 where Y1 and Y4 are the interval boundaries cornering the minimum 
  !
  !
  !
        
  implicit none

  real(4), intent(in)     :: a,b,ystart,eps
  real(4), external       :: func     
        
  real(4), intent(out)    :: x_min, y_min
  integer(4), intent(out) :: jflag

  real(4)                 :: y1,y2,y3,y4
  real(4)                 :: fy1,fy2,fy3,fy4,fystart
  integer(4)              :: i    

  integer(4), parameter   :: itermax=50 ! max number of iterations in minimum search
  integer(4), parameter   :: nmax = itermax         ! maximum number of iterations

  !      Golden Section Ratio
  real(4), parameter      :: alpha = 0.61803399, beta = 1. - alpha 

  y1 = a
  y4 = b
  fy1 = func(y1)
  fy4 = func(y4)
  fystart = func(ystart)
  if (fystart >= fy1 .or. fystart >= fy4 ) then
     jflag = -1
     return
  endif

  if ( abs(b-ystart) > abs(a-ystart))  then
     y3 = ystart + beta*(b-ystart)        ! determine new point 
     y2 = ystart
  else              
     y3 = ystart  
     y2 = a + alpha*(ystart - a)
  endif

  fy2 = func(y2)
  fy3 = func(y3)


  do i =1,nmax
     
     if ( abs(y4-y1) < eps*(abs(y2)+abs(y3)) )  then
        jflag = 0
        exit
     endif

     if (i == nmax) then
        jflag = 1 ! accuracy not achieved 
        exit
     endif

     if (fy3 < fy2) then
        y1 = y2
        y2 = y3
        y3 = alpha*y2 + beta*y4
        ! y4 stays
        fy2 = fy3
        fy3 = func(y3)
     else
        y4 = y3
        y3 = y2
        y2 = alpha*y3 + beta*y1
        ! y1 stays
        fy3 = fy2
        fy2 = func(y2)
     endif

  enddo
        

  if (fy2 < fy3) then
     x_min = y2
     y_min = fy2
  else
     x_min = y3
     y_min = fy3
  endif


  return                                                   
end subroutine minimum_search

end module wind_speed_retrieval_module
