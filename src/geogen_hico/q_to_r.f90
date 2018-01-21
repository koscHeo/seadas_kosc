module quat_to_rot
  implicit none
  private
  public:: q_to_r

contains
  subroutine q_to_r( q, a )
!*****************************************************************************80
!
! Very slightly modified my Marcos Montes. 
! MJM Note: scalar is q(1).
!! ROTATION_QUAT2MAT_3D converts rotation from quaternion to matrix form in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries van Dam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Q(4), the quaternion representing the rotation.
!
!    Output, real ( kind = 8 ) A(3,3), the rotation matrix.
!
    implicit none

    integer ( kind = 4 ), parameter :: dim_num = 3
    real ( kind = 8 ),intent(out),dimension(dim_num,dim_num):: a
    real ( kind = 8 ),intent(in):: q(4)
! Locals
    real ( kind = 8 ):: angle
    real ( kind = 8 ):: ca
    real ( kind = 8 ):: cos_phi
    real ( kind = 8 ):: sa
    real ( kind = 8 ):: sin_phi
    real ( kind = 8 ):: v1
    real ( kind = 8 ):: v2
    real ( kind = 8 ):: v3
    
    sin_phi = sqrt ( sum ( q(2:4)**2 ) )
    
    cos_phi = q(1)
    
    angle = 2.0D+00 * atan2 ( sin_phi, cos_phi )
    
    if ( sin_phi == 0.0D+00 ) then
       v1 = 1.0D+00
       v2 = 0.0D+00
       v3 = 0.0D+00
    else
       v1 = q(2) / sin_phi
       v2 = q(3) / sin_phi
       v3 = q(4) / sin_phi
    end if
    
    ca = cos ( angle )
    sa = sin ( angle )
    
    a(1,1) =                    v1 * v1 + ca * ( 1.0D+00 - v1 * v1 )
    a(1,2) = ( 1.0D+00 - ca ) * v1 * v2 - sa * v3
    a(1,3) = ( 1.0D+00 - ca ) * v1 * v3 + sa * v2
    
    a(2,1) = ( 1.0D+00 - ca ) * v2 * v1 + sa * v3
    a(2,2) =                    v2 * v2 + ca * ( 1.0D+00 - v2 * v2 )
    a(2,3) = ( 1.0D+00 - ca ) * v2 * v3 - sa * v1
    
    a(3,1) = ( 1.0D+00 - ca ) * v3 * v1 - sa * v2
    a(3,2) = ( 1.0D+00 - ca ) * v3 * v2 + sa * v1
    a(3,3) =                    v3 * v3 + ca * ( 1.0D+00 - v3 * v3 )
    
    return
  end subroutine q_to_r

end module quat_to_rot
