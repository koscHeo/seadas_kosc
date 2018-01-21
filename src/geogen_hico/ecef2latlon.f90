module ecef2llh
  implicit none
  private
  public::ecef2latlon
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ecef2latlon(xyz,llh)
    ! Need to check this set of equations; it seems pretty good,
    ! but now there are different equations on Wikipedia. See the papers
    ! I downloaded to the ReferenceFrames folder. 

    ! MJM copied algorithm from web page
    ! http://en.wikipedia.org/wiki/Geodetic_system
    ! COnverted to IDL 2009 Sept. 29.
    ! Converted to F90 2012 Oct 25
    ! 2012/10/31 Converted to Vermeille's algorithm since at least I have a referece. But performance looks
    ! exactly the same. 
    ! 2012/10/31 Changed output order to LON (deg), LAT (DEG), height (m)
    !
    ! ASSUMES WGS 84 ellipsoid
    ! Results are GEODETIC latitude and longitude
    ! These will be different from Geocentric latitude
    !implicit none
    real (kind=8),intent(in),dimension(:,:)::xyz ! ECEF, in meters
    real (kind=8),intent(out),dimension(size(xyz,1),size(xyz,2))::llh ! lon, lat, height
    real (kind=8),parameter::f=1.d0/298.257223563d0
    real (kind=8),parameter::a=6378137d0    
    real (kind=8),parameter::b=a*(1-f),e2= 2*f-f**2,ep2= f*(2-f)/((1-f)**2),ee2= a**2 - b**2
    real (kind=8),parameter:: e=sqrt(e2), a2=a**2, e4=e2**2
!    real (kind=8),dimension(size(xyz,2))::r2,r,ff,G,c,s,p,q,ro,tmp,u,v,zo
    real (kind=8),dimension(size(xyz,2))::p,q,r,s,t,u,v,w,k,d,tmp
    
    include 'astmath.cmn'

!    b = a*(1-f) !% semi-minor axis    
!    e2 = 2*f-f**2 !% first eccentricity squared
!    ep2 = f*(2-f)/((1-f)**2) ! % second eccentricity squared
!    EE2 = a**2 - b**2 !
! ()'s Algorithm????
!!    r2 = XYZ(1,:)**2 + XYZ(2,:)**2!
!!    r = sqrt(r2)  !
!!    FF = 54*b**2* XYZ(3,:)**2 !
!!    G = r2 + (1-e2)*XYZ(3,:)**2 - e2*EE2 !
!!    c = (e2*e2*FF*r2)/(G**3) !
!!    s = ( 1.d0 + c + sqrt(c*c + 2*c) )**(1.d0/3.d0) !
!!    P = FF/(3*(s + 1.d0/s + 1.d0 )**2*G**2)!
!!    Q = sqrt(1.d0 + 2.d0 *e2*e2*P)
!!    ro = -(e2*P*r)/(1+Q) + sqrt((a*a/2)*(1 + 1.d0/Q) - ((1-e2)*P*xyZ(3,:)**2)/(Q*(1+Q)) - P*r2/2)
!!    tmp = (r - e2*ro)**2
!!    U = sqrt( tmp + XYZ(3,:)**2 )
!!    V = sqrt( tmp + (1.d0 - e2)*XYZ(3,:)**2 )
!!    zo = (b**2*XYZ(3,:))/(a*V)

!!    llh(1,:) = atan( (XYZ(3,:) + ep2*zo)/r )*rad2deg
!!    llh(2,:) = atan(XYZ(2,:),XYZ(1,:)) *rad2deg
!!    llh(3,:) = U*( 1 - b**2/(a*V))

! Vermeille's Algorithm (Journal of Geodesy (2002) 76:451-454
!    e=sqrt(e2)
!    a2=a**2
    p=(XYZ(1,:)**2 + XYZ(2,:)**2)/a2
    q=(1-e2)/a2*XYZ(3,:)**2
    r=(p+q-e4)/6
    s=e4*p*q/4/r**3
    t=(1.d0 + s + sqrt(s*(2.d0 + s)))**(1.d0/3.d0)
    u=r*(1.d0 + t + 1.d0/t)
    v=sqrt(u**2+q*e4)
    w=e2*(u+v-q)/2/v
    k=sqrt(u+v+w**2)-w
    D=k*sqrt(xyz(1,:)**2 + xyz(2,:)**2)/(k+e2)
!    
!    llh(1,:)=2*atan(xyz(2,:),(xyz(1,:)+sqrt(xyz(1,:)**2 + xyz(2,:)**2)))
!    or better 
    tmp=sqrt(d**2 + xyz(3,:)**2)

    llh(1,:) = datan2(xyz(2,:),xyz(1,:)) ! longitude in rad
!    llh(2,:) = 2*atan(xyz(3,:),(d+sqrt(d**2+xyz(3,:)**2)))*rad2deg
    llh(2,:) = 2*datan2(xyz(3,:),(d + tmp))  ! latitude in rad
    llh(3,:) = tmp*(k+e2-1)/k   ! height in same units as "a", i.e., meters

  end subroutine ecef2latlon
end module ecef2llh
