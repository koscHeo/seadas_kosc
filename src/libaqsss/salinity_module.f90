! module for salinity retrieval
! November 10, 2011

module salinity_module
implicit none
private
public :: fdem0_meissner_wentz_salinity 
public :: dielectric_meissner_wentz
save

real(4), parameter                                      :: freq=1.413
integer(4), parameter :: n_rad=3             ! number of radiometers (i.e., horns)  (inner, middle, outer)	
real(4),    parameter :: missing_val=-9999.  ! indicator for missing value

! for chi2 minimization to determine sss
real(4), parameter                                      :: sss_min=0.0, sss_max=45.0, sss_0=35.0, eps=5.0E-4
integer(4), parameter                           :: itermax=50 ! max number of iterations in minimum search
real(4), dimension(2)                       :: ytb0 ! global variable accessed by chi2
real(4)                                                         :: y_sst, y_tht ! global variable accessed by chi2   
integer(4)                                                      :: y_irad

real(4), parameter, dimension(2,n_rad) :: del_tb0 = reshape( (/0.265,0.220, 0.282,0.209,  0.288,0.205/), (/2,n_rad/) )
! error estimate for tbsur0


public  :: fd_sss


contains


subroutine fd_sss(irad,tht,sst,surtb_mea, sss, chisq_sss, iflag_sss)
use wind_speed_retrieval_module

implicit none

integer(4), intent(in )  :: irad
real(4),    intent(in )  :: tht,sst,surtb_mea(2)
        
real(4),    intent(out)  :: sss, chisq_sss
integer(4), intent(out)  :: iflag_sss

real(4)                                  :: x_min, y_min
integer(4)                               :: m_search


! save global variables to be accesses through chi2
y_irad=irad
y_tht =tht
y_sst =sst
ytb0  =surtb_mea

m_search=0

call minimum_search(sss_min, sss_0, sss_max,chi2,eps,  x_min,y_min,m_search)
if (m_search/=0) then ! accuracy goal in minimum search was not achieved
   iflag_sss=1
   sss=missing_val
   chisq_sss=missing_val
else 
   iflag_sss=0 ! o.k.
   sss=x_min
   chisq_sss=y_min
endif

return
end subroutine fd_sss


real(4) function chi2(sss)
! use V and H pol

!     sos evaluation
!     fixed sst
!     variable sss

implicit none
real(4)                          :: sss
real(4), dimension(2)                    :: ze0, ztb0
real(4)                                                  :: zsst, ztht, surtep

zsst=y_sst
ztht=y_tht
surtep=(zsst+273.15)

call fdem0_meissner_wentz_salinity(freq,ztht,zsst,sss, ze0)
ztb0 = ze0*surtep
chi2 = ((ytb0(1)-ztb0(1))/del_tb0(1,y_irad))**2 + ((ytb0(2)-ztb0(2))/del_tb0(2,y_irad))**2  
     
return
end function chi2


subroutine fdem0_meissner_wentz_salinity(freq,tht,sst,salinity, em0) 
implicit none

real(4), parameter :: f0=17.97510
 
real(4) :: freq,tht,sst,salinity,em0(2)
real(4) :: costht,sinsqtht
real(4) :: e0s,e1s,e2s,n1s,n2s,sig
 
complex(4) :: permit,esqrt,rh,rv
 
complex(4), parameter :: j= (0.,1.)

real(4) cosd

call dielectric_meissner_wentz(sst,salinity,  e0s,e1s,e2s,n1s,n2s,sig)

costht=cosd(tht)
sinsqtht=1.-costht*costht

!     debye law (2 relaxation wavelengths)
permit = (e0s - e1s)/(1.0 - j*(freq/n1s)) + (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s + j*sig*f0/freq
permit = conjg(permit)
        
esqrt=csqrt(permit-sinsqtht)
rh=(costht-esqrt)/(costht+esqrt)
rv=(permit*costht-esqrt)/(permit*costht+esqrt)
em0(1)  =1.-rv*conjg(rv)
em0(2)  =1.-rh*conjg(rh)
 
return
end     subroutine fdem0_meissner_wentz_salinity


!     this routine was carefully editted by fjw to separate out the sst,s dependence
!     i also permultiplied sst**2, sst**3 etc.
!     i numerically verified that the editted version gave exactly the same results are thomas orginal

!     complex dielectric constant: eps
!     t. meissner + f. wentz, ieee tgars, 42(9), 2004, 1836-1849
!
!     input:
!     name   parameter  unit  range
!     sst      sst        [c]   -25 c to 40 c for pure water
!                               -2  c to 34 c for saline water
!     s      salinity   [ppt]  0 to 40
!
!     output:
!     eps    complex dielectric constant
!            negative imaginary part to be consistent with wentz1 convention

subroutine dielectric_meissner_wentz(sst_in,s,   e0s,e1s,e2s,n1s,n2s,sig)
implicit none

real(4), intent(in) :: sst_in,s
real(4), intent(out) :: e0s,e1s,e2s,n1s,n2s,sig
 
real(4), dimension(11), parameter :: &    
     x=(/ 5.7230e+00, 2.2379e-02, -7.1237e-04, 5.0478e+00, -7.0315e-02, 6.0059e-04, 3.6143e+00, &
     2.8841e-02, 1.3652e-01,  1.4825e-03, 2.4166e-04 /)
 
 
real(4), dimension(13), parameter :: &
     z=(/ -3.56417e-03,  4.74868e-06,  1.15574e-05,  2.39357e-03, -3.13530e-05, &
     2.52477e-07, -6.28908e-03,  1.76032e-04, -9.22144e-05, -1.99723e-02, &
     1.81176e-04, -2.04265e-03,  1.57883e-04  /)  ! 2004

real(4), dimension(3), parameter :: a0coef=(/ -0.33330E-02,  4.74868e-06,  0.0e+00/)
real(4), dimension(5), parameter :: b1coef=(/0.23232E-02, -0.79208E-04, 0.36764E-05, -0.35594E-06, 0.89795E-08/)
 
real(4) :: e0,e1,e2,n1,n2
real(4) :: a0,a1,a2,b1,b2
real(4) :: sig35,r15,rtr15,alpha0,alpha1

real(4) sst,sst2,sst3,sst4,s2

sst=sst_in
if(sst.lt.-30.16) sst=-30.16  !protects against n1 and n2 going zero for very cold water

sst2=sst*sst
sst3=sst2*sst
sst4=sst3*sst

s2=s*s
 
!   pure water
e0    = (3.70886e4 - 8.2168e1*sst)/(4.21854e2 + sst) ! stogryn et al.
e1    = x(1) + x(2)*sst + x(3)*sst2
n1    = (45.00 + sst)/(x(4) + x(5)*sst + x(6)*sst2)
e2    = x(7) + x(8)*sst
n2    = (45.00 + sst)/(x(9) + x(10)*sst + x(11)*sst2)
 
!   saline water
!   conductivity [s/m] taken from stogryn et al.
sig35 = 2.903602 + 8.60700e-2*sst + 4.738817e-4*sst2 - 2.9910e-6*sst3 + 4.3047e-9*sst4
r15   = s*(37.5109+5.45216*s+1.4409e-2*s2)/(1004.75+182.283*s+s2)
 
alpha0 = (6.9431+3.2841*s-9.9486e-2*s2)/(84.850+69.024*s+s2)
alpha1 = 49.843 - 0.2276*s + 0.198e-2*s2
rtr15 = 1.0 + (sst-15.0)*alpha0/(alpha1+sst)
 
sig = sig35*r15*rtr15
 
!   permittivity
a0 = exp(a0coef(1)*s + a0coef(2)*s2 + a0coef(3)*s*sst)  
e0s = a0*e0

if(sst.le.30) then
   b1 = 1.0 + s*(b1coef(1) + b1coef(2)*sst + b1coef(3)*sst2 + b1coef(4)*sst3 + b1coef(5)*sst4)
else
   b1 = 1.0 + s*(9.1873715e-04 + 1.5012396e-04*(sst-30))
endif
    
n1s = n1*b1

a1  = exp(z(7)*s + z(8)*s2 + z(9)*s*sst)
e1s = e1*a1

!   b2 = 1.0 + s*(z(10) + z(11)*sst)
b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30))
n2s = n2*b2

a2 = 1.0  + s*(z(12) + z(13)*sst)
e2s = e2*a2
 
return
end subroutine  dielectric_meissner_wentz

end module salinity_module


