! This wrapper is needed because the query_disc I/O interface contains 
! assumed-shape arrays which cannot be passed directly by C routines
subroutine query_disc_wrapper( nside, vector, radius, listpix, nlist)
  use pix_tools 

  integer(kind=4), intent(in)                 :: nside
  real(kind=8),    intent(in), dimension(3)   :: vector
  real(kind=8),    intent(in)                 :: radius
  integer(kind=4), intent(out), dimension(0:nlist-1) :: listpix
  integer(kind=4), intent(inout)                :: nlist
  
  call query_disc(nside, vector, radius, listpix, nlist)

end subroutine query_disc_wrapper

