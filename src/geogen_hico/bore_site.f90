module bore_site
  ! MJM 2012/11/14-16 
  implicit none
  private
  public::bore_sight
contains
  pure subroutine bore_sight(dthx,dthy,dthz,rot)
    implicit none
    real (kind=8), intent(in)::dthx,dthy,dthz ! RADIANS
    real (kind=8), intent(out),dimension(3,3)::rot
    real (kind=8)::cx,sx,cy,sy,cz,sz
    ! Really, this just constructs an Euler Angle Sequence (1,2,3) Deibel 8.2 (p24) rotation matrix
    ! This matrix looks like it goes from ISS to HICO; I'll need transpose in caller since I go from HICO to ISS.
    cx=cos(dthx)
    sx=sin(dthx)
    cy=cos(dthy)
    sy=sin(dthy)
    cz=cos(dthz)
    sz=sin(dthz)
    
    rot(:,1)=[cy*cz,cy*sz,-sy]
    rot(:,2)=[sx*sy*cz-cx*sz,sx*sy*sz+cx*cz,cy*sx]
    rot(:,3)=[cx*sy*cz+sx*sz,cx*sy*sz-sx*cz,cy*cx]
    
  end subroutine bore_sight
end module bore_site
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
