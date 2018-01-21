        subroutine fd_land_refl(itype,tht,sm, refl)
        implicit none
        
        integer(4) itype
        real(4) tht,sm,refl(2)
        real(4) acoef,alpha,beta
        real(4) costht, sinsqtht
        real(4) cosd
        complex(4) bcoef,rv,rh, esqrt, permit, permit_alpha
 
        data acoef/1.850842/               !from fd_soil_coefs.f assuming freq= 1.413
        data bcoef/(16.06828,-0.7877967)/ ! from fd_soil_coefs.f assuming freq= 1.413
        data alpha/0.65/
        data beta/1.09/

        if(itype.lt.1 .or. itype.gt.3) stop 'itype oob in fd_land_refl, pgm stopped'
 
c     this is the ulaby,moore,fumg model, vol 3 page 2103 that sab used.  also see sab doc on o:\aquarius\orbit_simulator\surface_atmos
      if(itype.eq.1) then  !soil above freezing
        permit_alpha=acoef + bcoef*sm**beta
        permit = permit_alpha**(1./alpha)  !soil
        endif

c     page 2026 and 2027 give the following for fresh-water ice at 1.4 ghz and t=-20:   permit_ice=(3.15,-0.0005)
c     fig e.26 on page 2062 gives real part of snow to be about 1.73 for snow density of 0.4
c     fig e.28 on page 2065 gives imag part of snow to be 0.23*(imag part of ice)=0.23*0.0005  for snow density of 0.4
        if(itype.eq.2)  permit=(1.73,-0.00012)  !snow

        if(itype.eq.3)  permit=(3.,-0.1)        !sea ice 
 
        costht=cosd(tht)
      sinsqtht=1.-costht*costht
        
      esqrt=csqrt(permit-sinsqtht)
      rh=(costht-esqrt)/(costht+esqrt)
      rv=(permit*costht-esqrt)/(permit*costht+esqrt)
      refl(1) =rv*conjg(rv)
      refl(2) =rh*conjg(rh)
 
      return
      end       
        
 
