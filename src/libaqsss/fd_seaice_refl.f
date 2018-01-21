        subroutine fd_seaice_refl(tht, refl)
        implicit none
        
      real(4), intent(in)   ::  tht
        real(4), intent(out)  :: refl(2)

        real(4) costht, sinsqtht
      complex(4) rv,rh, esqrt, permit
 
        real(4) cosd

        permit=(3.,-0.1)        !sea ice 
 
        costht=cosd(tht)
      sinsqtht=1.-costht*costht
        
      esqrt=csqrt(permit-sinsqtht)
      rh=(costht-esqrt)/(costht+esqrt)
      rv=(permit*costht-esqrt)/(permit*costht+esqrt)
      refl(1) =rv*conjg(rv)
      refl(2) =rh*conjg(rh)
 
      return
      end       
        
 
