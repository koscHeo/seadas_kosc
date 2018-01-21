!    TM updatged 11/12/2011: 
!    Include new surface roughness correction model
!    Use wind speed only but no sigma0. 

!     phir.le.-999 special value to indicate isotropic refl is to be returned

      subroutine fd_water_refl_exact(tht,sss,sst,wind,phir,staticdata,refl)
      use static_data_module
      use salinity_module
      implicit none

      real(4),  intent(in)  :: tht,sss,sst,wind,phir
      type(static_data_struct), intent(in) :: staticdata

      real(4),  intent(out) :: refl(2)


      real(4) freq,dems(2),em0(2)
        
      data freq/1.413/

      call fdem0_meissner_wentz_salinity(freq,tht,sst,sss, em0)
      refl=1-em0

!     add wind induced reflective (minus emissivity)
      call fd_dems(tht,wind,phir, sst, staticdata, dems)
      refl=refl - dems
  
      return
      end


      subroutine fd_dems(tht, wind, phir, sst, staticdata, dems)
      use static_data_module

      implicit none

      real(4),  intent(in)  :: tht,wind,phir,sst
      type(static_data_struct), intent(in) :: staticdata

      real(4),  intent(out) :: dems(2)
      real(4)               :: aharm(0:2,2) !1=V, 2=H
      real(4)               :: daharm(0:2,2) !1=V, 2=H

      real(4)               :: dtb_rough(2)

      integer(4) irad,iradsv  
      real(4) thtfx(3),dtht,xmin
      real(4) cosd
        
      real(4), parameter :: teff=290. ! used for wind induced emissivity harmonics
      real(4), parameter :: ysss=35.0
      data thtfx/29.19,37.85,45.48/ !nominal inc. ang. for 3 horns

      integer, pointer :: nptr => null()

c     find horn that has earth inc. angle closest to tht
      xmin=1.d30
      do irad=1,n_rad
         dtht=abs(tht-thtfx(irad))
         if(dtht.lt.xmin) then
            xmin=dtht
            iradsv=irad
         endif
      enddo

      irad=iradsv

      if (staticdata%emiss_sst_sss .eq. 1) then
         call fd_TM_emiss_harmonics (irad, wind, sst, ysss, staticdata, aharm, daharm)
      else
         call fd_TM_emiss_harmonics (irad, wind, nptr, nptr, staticdata, aharm, daharm)
      endif
     
      dtb_rough(1:2) = aharm(0,1:2)
      if (phir .gt. -998) then  ! also include wind direction signal
         dtb_rough(1:2) = dtb_rough(1:2) + aharm(1,1:2)*cosd(phir) + aharm(2,1:2)*cosd(2.0*phir)
      endif      
     
      dems=dtb_rough/teff ! V2.0  08/29/12
      
      return
      end
