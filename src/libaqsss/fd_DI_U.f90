! Thomas Meissner
! 04/08/2015

subroutine find_DI_U(irad,U,coeffs_dI_U, DI)
use static_data_module
implicit none
integer(4), intent(in)  :: irad
real(4), intent(in)     :: U ! TF(3) 3rd Stokes antenna temperature [K]
real(4), intent(in)     :: coeffs_dI_U(4,n_rad)
real(4), intent(out)    :: DI ! adjustgment in TOI I/2 = (V+H)/2 [K]

! computed from routine analyze_IU_coupling.f90 and fitted in IDL routine plot_DIDQ_U.pro
real(4), dimension(4) :: x, c
real(4)               :: y

x = (/U, U**2, U**3, U**4 /)
c(1:4) = coeffs_dI_U(1:4,irad)

y = dot_product(x,c)

dI = y

return
end subroutine find_di_u
