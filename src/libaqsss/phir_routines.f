!     08/09/ 2012: I changed w0 from 25.5 to 28.5
! 
!     06/07/2012

!     V6 harmonics form for emissivity
!     V5 harmonics form for sigma0

      subroutine fd_TM_emiss_harmonics(irad,wspd,sst,sss,staticdata,aharm,daharm)
      use static_data_module
      use salinity_module

! harmonic coefficients for wind induced emissivity
      implicit none

      integer(4), intent(in)                            :: irad
      real(4), intent(in)                               :: wspd
      real(4), intent(in), optional                     :: sst ! Celsius
      real(4), intent(in), optional                     :: sss ! psu

      type(static_data_struct), intent(in) :: staticdata

      real(4), dimension(0:2,2), intent(out)    :: aharm !1=V, 2=H
      real(4), dimension(0:2,2), intent(out), optional    :: daharm !1=V, 2=H

      real(4), dimension(2)    :: A0, A1, A2 !1=V, 2=H
      real(4), dimension(2)    :: dA0, dA1, dA2 !1=V, 2=H
        
      real(4)                  :: ww
      integer(4)               :: ipol, iharm
      real(4)                  :: fval, dval

      real(4)                  :: xsst, xsss, xtht
      real(4), dimension(2)    :: xem0, yem0

      real(4)                  :: w0, w1, w2 ! linear extrapolation/cutoff points
      integer(4), parameter    :: npoly=5 
      real(4), parameter       :: freq_aq = 1.41
      real(4), parameter       :: sst_ref=20.0
      real(4), dimension(n_rad):: tht_ref = (/29.19, 37.85, 45.48/)

      A0=0.0
      A1=0.0
      A2=0.0

      do ipol=1,2 

! A0
         iharm=0
         w0 = staticdata%wspd_max_a(iharm,ipol,irad)
         ww = wspd
         if (wspd >= w0) ww=w0  ! extrapolation at w0 
         fval = ww    *staticdata%acoef(iharm,ipol,irad,1) + (ww**2)*staticdata%acoef(iharm,ipol,irad,2) + 
     1        (ww**3)*staticdata%acoef(iharm,ipol,irad,3) + (ww**4)*staticdata%acoef(iharm,ipol,irad,4) + 
     2        (ww**5)*staticdata%acoef(iharm,ipol,irad,5)  
         dval =             staticdata%acoef(iharm,ipol,irad,1) + (2.0*ww)     *staticdata%acoef(iharm,ipol,irad,2) + 
     1        (3.0*(ww**2))*staticdata%acoef(iharm,ipol,irad,3) + (4.0*(ww**3))*staticdata%acoef(iharm,ipol,irad,4) + 
     2        (5.0*(ww**4))*staticdata%acoef(iharm,ipol,irad,5) 
         if (wspd<=w0) then
            A0(ipol) = fval
         else
            A0(ipol) = fval + dval*(wspd-w0)
         endif

         dA0(ipol) = dval

! A1
         iharm=1
         w1 = staticdata%wspd_max_a(iharm,ipol,irad)
         ww = wspd
         if (wspd >= w1) ww=w1  ! cutoff at w1 
         fval = ww    *staticdata%acoef(iharm,ipol,irad,1) + (ww**2)*staticdata%acoef(iharm,ipol,irad,2) + 
     1        (ww**3)*staticdata%acoef(iharm,ipol,irad,3) + (ww**4)*staticdata%acoef(iharm,ipol,irad,4) + 
     2        (ww**5)*staticdata%acoef(iharm,ipol,irad,5)  
         dval =             staticdata%acoef(iharm,ipol,irad,1) + (2.0*ww)     *staticdata%acoef(iharm,ipol,irad,2) + 
     1        (3.0*(ww**2))*staticdata%acoef(iharm,ipol,irad,3) + (4.0*(ww**3))*staticdata%acoef(iharm,ipol,irad,4) + 
     2        (5.0*(ww**4))*staticdata%acoef(iharm,ipol,irad,5) 
         A1(ipol) = fval
         dA1(ipol) = dval

! A2
         iharm=2
         w2 =staticdata%wspd_max_a(iharm,ipol,irad)
         ww = wspd
         if (wspd >= w2) ww=w2  ! cutoff at w2 
         fval = ww    *staticdata%acoef(iharm,ipol,irad,1) + (ww**2)*staticdata%acoef(iharm,ipol,irad,2) + 
     1        (ww**3)*staticdata%acoef(iharm,ipol,irad,3) + (ww**4)*staticdata%acoef(iharm,ipol,irad,4) + 
     2        (ww**5)*staticdata%acoef(iharm,ipol,irad,5)  
         dval =             staticdata%acoef(iharm,ipol,irad,1) + (2.0*ww)     *staticdata%acoef(iharm,ipol,irad,2) + 
     1        (3.0*(ww**2))*staticdata%acoef(iharm,ipol,irad,3) + (4.0*(ww**3))*staticdata%acoef(iharm,ipol,irad,4) + 
     2        (5.0*(ww**4))*staticdata%acoef(iharm,ipol,irad,5) 
         A2(ipol) = fval
         dA2(ipol) = dval
        
      enddo                     !ipol

! The aharm are multiplied by E0(SST)/E0(20)
      if (loc(sst).ne.0 .and. loc(sss).ne.0) then
         xtht=tht_ref(irad)
         xsst=sst
         xsss=sss

         call fdem0_meissner_wentz_salinity(freq_aq,xtht,xsst,   xsss, xem0) 
         call fdem0_meissner_wentz_salinity(freq_aq,xtht,sst_ref,xsss, yem0) 

         A0(1:2) = A0(1:2)*xem0(1:2)/yem0(1:2)
         A1(1:2) = A1(1:2)*xem0(1:2)/yem0(1:2)
         A2(1:2) = A2(1:2)*xem0(1:2)/yem0(1:2)
            
         dA0(1:2) = dA0(1:2)*xem0(1:2)/yem0(1:2)
         dA1(1:2) = dA1(1:2)*xem0(1:2)/yem0(1:2)
         dA2(1:2) = dA2(1:2)*xem0(1:2)/yem0(1:2)
      endif 


      do ipol=1,2
         aharm(0,ipol) = A0(ipol)
         aharm(1,ipol) = A1(ipol)
         aharm(2,ipol) = A2(ipol)
      enddo

!      if (present(dstaticdata%aharm)) then
!         do ipol=1,2
!            dstaticdata%aharm(0,ipol) = dA0(ipol)
!            dstaticdata%aharm(1,ipol) = dA1(ipol)
!            dstaticdata%aharm(2,ipol) = dA2(ipol)
!         enddo
!      endif

      return
      end subroutine fd_TM_emiss_harmonics


      subroutine fd_TM_scat_harmonics (irad, wspd, staticdata, bharm, dbharm)
        ! updated for high wind speeds June 6, 2012
        ! updated for V5 on April 27, 2012
        ! coefficients were fitted using WindSat wspd and NCEP wdir
        ! change of inter/extrapolation points

      use static_data_module  
        ! harmonic coefficients for sigma0
      implicit none

      integer(4), intent(in)                        :: irad
      real(4), intent(in)                           :: wspd
      type(static_data_struct) :: staticdata

      real(4), dimension(0:2,npol_sc), intent(out):: bharm !1=VV, 2=HH, 3=VH, 4=HV
      real(4), dimension(0:2,npol_sc), intent(out), optional  :: dbharm  

      real(4), dimension(npol_sc)                   :: B0, B1, B2 !1=VV, 2=HH, 3=VH, 4=HV
      real(4), dimension(npol_sc)                   :: dB0, dB1, dB2 !1=VV, 2=HH, 3=VH, 4=HV

      real(4)                 ::      ww
      integer(4)              ::      ipol, iharm
      real(4)                 ::      fval, dval

!        real(4), parameter      :: w2=21.0 ! linear extrapolation/cutoff point 
!        real(4), dimension(0:2), parameter      :: w2= (/25.5, 22.5, 22.5/) ! linear extrapolation/cutoff point
      real(4)             :: w2 ! linear extrapolation point
      real(4)             :: wcut ! cutoff point 
      integer(4), parameter :: npoly=5 

      B0=0.0
      B1=0.0
      B2=0.0

      dB0=0.0
      dB1=0.0
      dB2=0.0

      do ipol=1,npol_sc         ! V3 has also VH HV implemented

!     B0
         iharm=0
         w2 = staticdata%wspd_max_b(iharm,ipol,irad)
         wcut=30.0
         if (ipol==3 .or. ipol==4) wcut=35.
         ww = wspd
         if (wspd >= w2) ww=w2  ! linear extrapolation at w2 

         fval = ww    *staticdata%bcoef(iharm,ipol,irad,1) + (ww**2)*staticdata%bcoef(iharm,ipol,irad,2) + 
     1        (ww**3)*staticdata%bcoef(iharm,ipol,irad,3) + (ww**4)*staticdata%bcoef(iharm,ipol,irad,4) + 
     2        (ww**5)*staticdata%bcoef(iharm,ipol,irad,5)   
         dval =             staticdata%bcoef(iharm,ipol,irad,1) + (2.0*ww)     *staticdata%bcoef(iharm,ipol,irad,2) + 
     1        (3.0*(ww**2))*staticdata%bcoef(iharm,ipol,irad,3) + (4.0*(ww**3))*staticdata%bcoef(iharm,ipol,irad,4) + 
     2        (5.0*(ww**4))*staticdata%bcoef(iharm,ipol,irad,5)

         B0(ipol) = fval + dval*(wspd-ww) ! linear extrapolation 
! cutoff at wcut
         if (wspd>wcut) B0(ipol) = fval + dval*(wcut-ww) 
   
! ensure that B0 does not get negative
         if (B0(ipol)<=0.0) B0(ipol)=0.0
         dB0(ipol) = dval
     
! B1
         iharm=1
         w2 = staticdata%wspd_max_b(iharm,ipol,irad)
         wcut=25.0              ! I have lowered that from 30.0 in V9
         ww = wspd
         if (wspd >= w2) ww=w2  ! linear extrapolation at w2  

         fval = ww    *staticdata%bcoef(iharm,ipol,irad,1) + (ww**2)*staticdata%bcoef(iharm,ipol,irad,2) + 
     1        (ww**3)*staticdata%bcoef(iharm,ipol,irad,3) + (ww**4)*staticdata%bcoef(iharm,ipol,irad,4) + 
     2        (ww**5)*staticdata%bcoef(iharm,ipol,irad,5)   
         dval =             staticdata%bcoef(iharm,ipol,irad,1) + (2.0*ww)     *staticdata%bcoef(iharm,ipol,irad,2) + 
     1        (3.0*(ww**2))*staticdata%bcoef(iharm,ipol,irad,3) + (4.0*(ww**3))*staticdata%bcoef(iharm,ipol,irad,4) + 
     2        (5.0*(ww**4))*staticdata%bcoef(iharm,ipol,irad,5)

         B1(ipol) = fval + dval*(wspd-ww) ! linear extrapolation 
!     cutoff at wcut
         if (wspd>wcut) B1(ipol) = fval + dval*(wcut-ww) 
         dB1(ipol) = dval

!     B2
         iharm=2
         w2 = staticdata%wspd_max_b(iharm,ipol,irad)
         wcut=25.0              ! I have lowered that from 30.0 in V9
         ww = wspd
         if (wspd >= w2) ww=w2  ! linear extrapolation at w2 

         fval = ww    *staticdata%bcoef(iharm,ipol,irad,1) + (ww**2)*staticdata%bcoef(iharm,ipol,irad,2) + 
     1        (ww**3)*staticdata%bcoef(iharm,ipol,irad,3) + (ww**4)*staticdata%bcoef(iharm,ipol,irad,4) + 
     2        (ww**5)*staticdata%bcoef(iharm,ipol,irad,5)   
         dval =             staticdata%bcoef(iharm,ipol,irad,1) + (2.0*ww)     *staticdata%bcoef(iharm,ipol,irad,2) + 
     1        (3.0*(ww**2))*staticdata%bcoef(iharm,ipol,irad,3) + (4.0*(ww**3))*staticdata%bcoef(iharm,ipol,irad,4) + 
     2        (5.0*(ww**4))*staticdata%bcoef(iharm,ipol,irad,5)
         B2(ipol) = fval + dval*(wspd-ww) ! linear extrapolation 
!     cutoff at wcut
         if (wspd>wcut) B2(ipol) = fval + dval*(wcut-ww) 
         dB2(ipol) = dval
    
      enddo                     ! ipol


      do ipol=1,npol_sc
         bharm(0,ipol) = B0(ipol)
         bharm(1,ipol) = B1(ipol)
         bharm(2,ipol) = B2(ipol)
      enddo

      if (present(dbharm)) then
         do ipol=1,npol_sc
            dbharm(0,ipol) = dB0(ipol)
            dbharm(1,ipol) = dB1(ipol)
            dbharm(2,ipol) = dB2(ipol)
         enddo
      endif
      
      return
      end subroutine fd_TM_scat_harmonics
