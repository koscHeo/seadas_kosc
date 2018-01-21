module intercept
  implicit none
  private
  public::wgs84_intercept
contains
  subroutine wgs84_intercept(rsc,u,out,view_zen,view_az,altitude)
    use ecef2llh
    ! Thanks to beautiful form of equations from Reid Reynolds.
    !implicit none
    real (kind=8),intent(in),dimension(3)::rsc ! spacecraft radius
    real (kind=8),intent(in),dimension(:,:)::u ! unit vector pointing toward earth, a lot of them
    real (kind=8),intent(out),dimension(size(u,1),size(u,2))::out !lon,lat,height (m)
    real (kind=8),intent(out),dimension(size(u,2)),optional::view_zen,view_az
    real (kind=8),intent(in),optional::altitude ! altitude in meters (wrt WGS-84 ellipsoid) of target surface
    ! Parameters
    real (kind=8),parameter::ff=1.d0/298.257223563d0 !WGS-84
    real (kind=8),parameter::re=6378137d0    ! WGS-84
    ! F is really the diagonal elements of a 3x3 matrix, all rest are 0
    ! but we only need the diagonal for the WGS-84 ellipsoid fo rvery simple
    ! mathematics, below. 
    real (kind=8),parameter,dimension(3)::F=[1.d0,1.d0,1.d0/(1.d0-ff)**2]
    !more parameters
    real (kind=8), parameter:: Pi      =  4.0D0 * DATAN(1.0D0)
    real (kind=8), parameter:: HalfPi  =  0.5D0*pi
    real (kind=8), parameter:: TwoPi   =  2.0D0*pi
    real (kind=8), parameter:: Rad2Deg = 180.0D0/pi
    real (kind=8), parameter:: Deg2Rad = pi/180.0D0
    ! Locals
    real (kind=8)::c,alt
    real (kind=8),dimension(size(u,2))::a,b,s,uDotE,uDotN,det
    real (kind=8),dimension(size(u,1),size(u,2))::rout
    real (kind=8),dimension(size(u,1),size(u,2))::snorm
    integer::ns
!    real (kind=8),dimension(3,1)::ssllh

! Start...
    if (present(altitude)) then
       alt=altitude
    else
       alt=0.d0
    endif

! 

! Using F as above (but see the comments - really F as the 3x3 matrix)
! If Rg is a vector from the center to the surface of the earth, then
! the equation of the WGS-84 ellipsoid is transpose(Rg)F*Rg=Re*Re

! The equation of a ray with unit direction (vector) u from the spacecraft 
! (at vector position Rsc) to the earth is the vector equation:
!  Rg = Rsc + s*u ; s is length of ray from the spacecraft to earth; 
! substitue this vector equation into the one above, and use the form 
! of the various vectors and matrices to get the equations below, 
! which yields a simple quadratic equation for s. For us, Rsc is one 
! location, and we have 512 unit direction vectors (elements of u). 
    
! precalculate a few things, using the magic rules of sum () and  matmul() 

    c=sum(F*rsc**2)-re**2 ! scalar, POSITIVE since spacecraft is above earth surface

    b=2.d0*matmul(F*Rsc,u) ! ns elements; negative since we're looking at earth and this is 

    ! essentially a dot prod of spacecraft vector (away from center of earth) and view vector
    ! (towards earth). 
    
    a=matmul(F,u**2)  ! ns elements, positive, since it is like u^2. 

    det= b**2 - 4*a*c

    if (any(det.lt.0)) then
       print*,'ERROR IN WGS84_intercept: invalid answer. no intercept'
       print*,'CHECK INPUT!'
       print*,minval(det),minloc(det)
       print*,maxval(det),maxloc(det)
       stop
    endif

! Note that -b is positive. So the closest root, the one with the smallest s, 
! that is the smallest distance from the spacecraft (the one on the spacecraft 
! side of the earth) is the one with the negative sign before the sqrt(). The OTHER one (positive) 
! is when the ray emerges from the earth.
! Now the distance along the ray, from the spacecraft to the near intercept is:   
    s= (-b -sqrt(det))/(2.d0 * a ) ! ns elements
! Once you know s, rout is the the location in ECEF of the intercept with the Earth,
! traveling a distance s from rsc in the direction(s) u. 
    ns=size(u,2)
    rout= spread(rsc,dim=2,ncopies=ns) + spread(s,dim=1,ncopies=3)*u ! 3xns
    
! Now we need an ECEF->LATLONH conversion
    call ecef2latlon(rout,out) 


! The below essentially use terms from the ECEF -> ENU conversion on the surface of an oblate spheroid,
! our WGS-84 ellipsoid
! The normal to the ellipsoid has the normal snorm
    if (present(view_zen)) then 
       snorm(1,:)=cos(out(1,:))*cos(out(2,:))
       snorm(2,:)=sin(out(1,:))*cos(out(2,:))
       snorm(3,:)=sin(out(2,:))
! The cos(view zenith) is formed by the pointing vector from the ground to the spacecraft
! dotted into the surface normal; this is (for one location) sum(-u(:,i)*snorm(:,i)) 
       view_zen=acos(sum(-u*snorm,dim=1))*rad2deg ! deg from zenith
    endif

    if (present(view_az)) then 
       uDotE=sin(out(1,:))*u(1,:) - cos(out(1,:))*u(2,:)
       uDotN=cos(out(1,:))*sin(out(2,:))*u(1,:) + sin(out(1,:))*sin(out(2,:))*u(2,:) - &
            cos(out(2,:))*u(3,:)
       view_az=datan2(udotE,uDotN)*rad2deg ! deg from N, clockwise
    end if
    out=out*rad2deg

  end subroutine wgs84_intercept
end module intercept
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
