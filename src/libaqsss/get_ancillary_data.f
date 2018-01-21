C       Provide ancillary atmospheric information to retrieval
C       WIND: wind speed PHIW: wind direction
C       TRAN: 
C       TBUP: upwelling radiation by the atmosphere
C       TBDW: downwelling radiation scattered from the surface

      SUBROUTINE GET_ANCILLARY_DATA(LYEAR,IDAYJL,ISECDY,XLAT,XLON,THT, 
     1                              SURTEMP,SURP,CWAT,SWE,SM,subtemp,WIND,PHIW,
     2                              TRAN,TBUP,TBDW, SOLAR_FLUX,
     3                              irad, zang, sclon, fpt_lnd, frc_lnd, dtb,
     4                              fland, gland, fice, gice, waveht, sss_reference,
     5                              tb_landcorr_vpol, tb_landcorr_hpol,     
     6                              ancfile1, ancfile2, ancfile3)
        IMPLICIT NONE

        integer(4), parameter :: n_rad=3  
        integer(4), parameter :: nlon_lnd=2160, nzang_lnd=2881

        INTEGER(4) LYEAR,IDAYJL,ISECDY
        INTEGER(4) LYEARSV,IDAYJLSV,IHOURSV,LYEARA,LYEARB
        INTEGER(4) ILEAP,IHOUR,IHOUR1,LYEAR1,IDAYJL1,IHOUR2,LYEAR2,IDAYJL2
        INTEGER(4) ITIME,ITIME1,ITIME2,IDELTIM1,IDELTIM2
        INTEGER(4) J1,J2,K1,K2,L1,L2
        integer(4) ilat,ilon
        REAL(4) XLAT,XLON,THT,SURTEMP,SURP,CWAT,SM,WIND,PHIW,TRAN,TBUP,TBDW
        REAL(4) UU,VV,WW,subtemp
        REAL(4) A1,A2,B1,B2,C1,C2,D1,D2,BRIEF
        real(4) solar_flux
        real(4) gswe

        integer(4), intent(in)  :: irad
        real(8),    intent(in)  :: zang,sclon

        real(4),    intent(in)  :: fpt_lnd(nlon_lnd,nzang_lnd,n_rad)
        real(4),    intent(in)  :: frc_lnd(nlon_lnd,nzang_lnd,n_rad)
        integer(2), intent(in)  :: dtb(1440,1440,12,2,3) 
        real(4),    intent(out) :: fland, gland, fice, gice, swe, sss_reference
        real(4),    intent(out) :: waveht, tb_landcorr_vpol, tb_landcorr_hpol

        REAL(4) SURP1(360,181),CWAT1(360,181),SM1(360,181),subtemp1(360,181),
     1          UU1(360,181),VV1(360,181),WW1(360,181),SWE1(720,361)
        REAL(4) SURP2(360,181),CWAT2(360,181),SM2(360,181),subtemp2(360,181),
     1          UU2(360,181),VV2(360,181),WW2(360,181),SWE2(720,361)
        REAL(4) SURTEMP1(360*4,180*4), SURTEMP2(360*4,180*4)
        REAL(4) SSS1(360*4,180*4), SSS2(360*4,180*4)
        REAL(4) FICE1(360*2,180*2), FICE2(360*2,180*2)

        REAL(4) TRAN1(360,181,91),TBUP1(360,181,91),TBDW1(360,181,91)
        REAL(4) TRAN2(360,181,91),TBUP2(360,181,91),TBDW2(360,181,91)

        REAL(4)      ::  WAVEHT1(288,157),WAVEHT2(288,157)

        real(4) solarflux1, solarflux2

        real(4) atan2d

        real(4) surtempa, surtempb, sssa, sssb, wt(4)
        real(8) tsum(0:1)
        integer(4) k, j, n, kk

        real(8) xday,xmon

        character(256) ancfile1, ancfile2, ancfile3
        integer ipass

        INTERFACE
           SUBROUTINE READ_DATA(filename,
     1          SURTEMP,SURP,CWAT,SM,subtemp,
     2          UU,VV,WW,SWE,TRAN,TBUP,TBDW,fice, waveht, sss, SOLAR_FLUX)
           USE HDF5             ! This module contains all necessary modules    
           IMPLICIT NONE

           CHARACTER(256) FILENAME
           REAL(4),TARGET :: SURTEMP(360*4,180*4),SURP(360,181),CWAT(360,181),
     1          SM(360,181),UU(360,181),VV(360,181),WW(360,181),SWE(720,361),
     2          subtemp(360,181)
           REAL(4),TARGET :: TRAN(360,181,91),TBUP(360,181,91),TBDW(360,181,91),
     1          fice(360,181), waveht(288,157)
           REAL(4),target :: sss(360*4,180*4)
           REAL(4),target :: SOLAR_FLUX
           END SUBROUTINE READ_DATA
        END INTERFACE

        DATA LYEARSV,IDAYJLSV,IHOURSV/3*-999/
        DATA LYEARA,LYEARB/1975,1973/
        DATA ITIME1,ITIME2/2*-999/

        IF(LYEAR.EQ.0) THEN
           LYEARSV=-999; IDAYJLSV=-999; IHOURSV=-999
           RETURN
        ENDIF

        DATA ipass/1/

        ILEAP=0
        IF(LYEAR.EQ.4*INT(LYEAR/4)) ILEAP=1

        IHOUR=6*INT(ISECDY/21600)

        ! Change .GE. to .GT.
        ! Change stop to exit
        ! JMG 09/02/11
        IF(IHOUR.GT.24) then
           write(*,*) 'ERROR1 IN GCM WIND ROUTINE'
           write(*,*) 'IHOUR: ', ihour
           call exit(1)
        endif

!!!!!!!!! Begin get new y-ancillary file !!!!!!!!!
        IF(LYEAR.EQ.LYEARSV .AND. IDAYJL.EQ.IDAYJLSV .AND. IHOUR.EQ.IHOURSV) 
     1  GOTO 100
        LYEARSV=LYEAR
        IDAYJLSV=IDAYJL
        IHOURSV=IHOUR

        IHOUR1=IHOUR
        LYEAR1=LYEAR
        IDAYJL1=IDAYJL

        ITIME1=31536000*(LYEAR1-LYEARA) + 86400*(IDAYJL1-1) + 3600*IHOUR1 + 
     1         86400*INT((LYEAR1 - LYEARB)/4)
        if (ipass.eq.1) then
           CALL READ_DATA(ancfile1, SURTEMP1,SURP1,CWAT1,SM1,subtemp1,
     1          UU1,VV1,WW1,SWE1,TRAN1,TBUP1,TBDW1,fice1,waveht1,sss1,solarflux1)
        else
           CALL READ_DATA(ancfile2, SURTEMP1,SURP1,CWAT1,SM1,subtemp1,
     1          UU1,VV1,WW1,SWE1,TRAN1,TBUP1,TBDW1,fice1,waveht1,sss1,solarflux1)
        endif

        IHOUR2=IHOUR
        LYEAR2=LYEAR
        IDAYJL2=IDAYJL

        IHOUR2=IHOUR2+6
        IF(IHOUR2.EQ.24) THEN
           IHOUR2=0 
           IDAYJL2=IDAYJL2+1
           IF(IDAYJL2.GT.365+ILEAP) THEN
              IDAYJL2=1
              LYEAR2=LYEAR2+1
           ENDIF
        ENDIF
        ITIME2=31536000*(LYEAR2-LYEARA) + 86400*(IDAYJL2-1) + 3600*IHOUR2 + 
     1         86400*INT((LYEAR2 - LYEARB)/4)
        if (ipass.eq.1) then
           CALL READ_DATA(ancfile2, SURTEMP2,SURP2,CWAT2,SM2,subtemp2,
     1          UU2,VV2,WW2,SWE2,TRAN2,TBUP2,TBDW2,fice2,waveht2,sss2,solarflux2)
        else
           CALL READ_DATA(ancfile3, SURTEMP2,SURP2,CWAT2,SM2,subtemp2,
     1          UU2,VV2,WW2,SWE2,TRAN2,TBUP2,TBDW2,fice2,waveht2,sss2,solarflux2)
        endif

        ipass = ipass + 1

 100    CONTINUE
!!!!!!!!! End get new y-ancillary file !!!!!!!!!

        ITIME=31536000*(LYEAR-LYEARA) + 86400*(IDAYJL-1) + ISECDY +
     1        86400*INT((LYEAR - LYEARB)/4)

        IDELTIM1=ITIME1-ITIME
        IDELTIM2=ITIME2-ITIME
        IF(IDELTIM1.GT.0 .OR. IDELTIM2.LT.0) then
           write(*,*) 'ERROR1 GET_ANCILLARY_DATA, PROGRAM STOPPED'
           call exit(1)
        endif

        IF(IDELTIM2-IDELTIM1.NE.21600) then
           write(*,*) 'ERROR1 GET_ANCILLARY_DATA, PROGRAM STOPPED' 
           call exit(1)
        endif

        A2=(ITIME-ITIME1)/REAL(ITIME2-ITIME1)
        A1=1.-A2
                 
!       return if bad lat/lon  JMG  10/12/11
        IF(ABS(XLAT).GT.90.) then
           write(*,*) 'LAT OOB IN GET_ANCILLARY_DATA, PROGRAM STOPPED'
           return
        endif

        IF(XLON.LT.0..OR.XLON.GT.360.) then
           write(*,*) 'LON OOB IN GET_ANCILLARY_DATA, PROGRAM STOPPED'  
           return
        endif

C       inverse-distance interpolation. B1,B2,C1,C2 are weights and in [0,1].
        BRIEF=90.-XLAT 
        IF(BRIEF.GT.179.999) BRIEF=179.999 
        J1=INT(1+BRIEF)                                                        
        J2=J1+1 
        B1=J1-BRIEF                                                       
        B2=1.-B1
 
        BRIEF=XLON
        IF(BRIEF.GT.359.999) BRIEF=0.
        K1=INT(1+BRIEF)                                                       
        K2=K1+1                                                           
        IF(K2.EQ.361) K2=1                                                
        C1=K1-BRIEF                                                       
        C2=1.-C1  

        BRIEF=THT  
        IF(BRIEF.GT.89.999) BRIEF=89.999
        L1=INT(1+BRIEF) 
        L2=L1+1                                                           
        D1=L1-BRIEF                                                       
        D2=1.-D1 

C     The following could be sped up becuase the same multiplication are done again and again
        
!      SURTEMP=                                                              
!     &    A1*B1*(C1*SURTEMP1(K1,J1)+C2*SURTEMP1(K2,J1))+                 
!     &    A1*B2*(C1*SURTEMP1(K1,J2)+C2*SURTEMP1(K2,J2))+                 
!     &    A2*B1*(C1*SURTEMP2(K1,J1)+C2*SURTEMP2(K2,J1))+                 
!     &    A2*B2*(C1*SURTEMP2(K1,J2)+C2*SURTEMP2(K2,J2))  

      SURP=                                                              
     &    A1*B1*(C1*SURP1(K1,J1)+C2*SURP1(K2,J1))+                 
     &    A1*B2*(C1*SURP1(K1,J2)+C2*SURP1(K2,J2))+                 
     &    A2*B1*(C1*SURP2(K1,J1)+C2*SURP2(K2,J1))+                 
     &    A2*B2*(C1*SURP2(K1,J2)+C2*SURP2(K2,J2))  

      CWAT=                                                              
     &    A1*B1*(C1*CWAT1(K1,J1)+C2*CWAT1(K2,J1))+                 
     &    A1*B2*(C1*CWAT1(K1,J2)+C2*CWAT1(K2,J2))+                 
     &    A2*B1*(C1*CWAT2(K1,J1)+C2*CWAT2(K2,J1))+                 
     &    A2*B2*(C1*CWAT2(K1,J2)+C2*CWAT2(K2,J2))  

      SM=                                                              
     &    A1*B1*(C1*SM1(K1,J1)+C2*SM1(K2,J1))+                 
     &    A1*B2*(C1*SM1(K1,J2)+C2*SM1(K2,J2))+                 
     &    A2*B1*(C1*SM2(K1,J1)+C2*SM2(K2,J1))+                 
     &    A2*B2*(C1*SM2(K1,J2)+C2*SM2(K2,J2))  
        
      subtemp=                                                              
     &    A1*B1*(C1*subtemp1(K1,J1)+C2*subtemp1(K2,J1))+                 
     &    A1*B2*(C1*subtemp1(K1,J2)+C2*subtemp1(K2,J2))+                 
     &    A2*B1*(C1*subtemp2(K1,J1)+C2*subtemp2(K2,J1))+                 
     &    A2*B2*(C1*subtemp2(K1,J2)+C2*subtemp2(K2,J2))  
        
      UU=                                                              
     &    A1*B1*(C1*UU1(K1,J1)+C2*UU1(K2,J1))+                 
     &    A1*B2*(C1*UU1(K1,J2)+C2*UU1(K2,J2))+                 
     &    A2*B1*(C1*UU2(K1,J1)+C2*UU2(K2,J1))+                 
     &    A2*B2*(C1*UU2(K1,J2)+C2*UU2(K2,J2))  
        
      VV=                                                              
     &    A1*B1*(C1*VV1(K1,J1)+C2*VV1(K2,J1))+                 
     &    A1*B2*(C1*VV1(K1,J2)+C2*VV1(K2,J2))+                 
     &    A2*B1*(C1*VV2(K1,J1)+C2*VV2(K2,J1))+                 
     &    A2*B2*(C1*VV2(K1,J2)+C2*VV2(K2,J2))  
         
      WW=                                                              
     &    A1*B1*(C1*WW1(K1,J1)+C2*WW1(K2,J1))+                 
     &    A1*B2*(C1*WW1(K1,J2)+C2*WW1(K2,J2))+                 
     &    A2*B1*(C1*WW2(K1,J1)+C2*WW2(K2,J1))+                 
     &    A2*B2*(C1*WW2(K1,J2)+C2*WW2(K2,J2))                  

!      SWE=                                                              
!     &    A1*B1*(C1*SWE1(2*K1,2*J1)+C2*SWE1(2*K2,2*J1))+                 
!     &    A1*B2*(C1*SWE1(2*K1,2*J2)+C2*SWE1(2*K2,2*J2))+                 
!     &    A2*B1*(C1*SWE2(2*K1,2*J1)+C2*SWE2(2*K2,2*J1))+                 
!     &    A2*B2*(C1*SWE2(2*K1,2*J2)+C2*SWE2(2*K2,2*J2))                  

                                                                       
      TRAN=                                                              
     & D1*(A1*B1*(C1*TRAN1(K1,J1,L1)+C2*TRAN1(K2,J1,L1))+                 
     &     A1*B2*(C1*TRAN1(K1,J2,L1)+C2*TRAN1(K2,J2,L1))+                 
     &     A2*B1*(C1*TRAN2(K1,J1,L1)+C2*TRAN2(K2,J1,L1))+                 
     &     A2*B2*(C1*TRAN2(K1,J2,L1)+C2*TRAN2(K2,J2,L1)))+  
     & D2*(A1*B1*(C1*TRAN1(K1,J1,L2)+C2*TRAN1(K2,J1,L2))+                 
     &     A1*B2*(C1*TRAN1(K1,J2,L2)+C2*TRAN1(K2,J2,L2))+                 
     &     A2*B1*(C1*TRAN2(K1,J1,L2)+C2*TRAN2(K2,J1,L2))+                 
     &     A2*B2*(C1*TRAN2(K1,J2,L2)+C2*TRAN2(K2,J2,L2)))  
                                                                       
      TBUP=                                                              
     & D1*(A1*B1*(C1*TBUP1(K1,J1,L1)+C2*TBUP1(K2,J1,L1))+                 
     &     A1*B2*(C1*TBUP1(K1,J2,L1)+C2*TBUP1(K2,J2,L1))+                 
     &     A2*B1*(C1*TBUP2(K1,J1,L1)+C2*TBUP2(K2,J1,L1))+                 
     &     A2*B2*(C1*TBUP2(K1,J2,L1)+C2*TBUP2(K2,J2,L1)))+  
     & D2*(A1*B1*(C1*TBUP1(K1,J1,L2)+C2*TBUP1(K2,J1,L2))+                 
     &     A1*B2*(C1*TBUP1(K1,J2,L2)+C2*TBUP1(K2,J2,L2))+                 
     &     A2*B1*(C1*TBUP2(K1,J1,L2)+C2*TBUP2(K2,J1,L2))+                 
     &     A2*B2*(C1*TBUP2(K1,J2,L2)+C2*TBUP2(K2,J2,L2)))  
                                                                       
      TBDW=                                                              
     & D1*(A1*B1*(C1*TBDW1(K1,J1,L1)+C2*TBDW1(K2,J1,L1))+                 
     &     A1*B2*(C1*TBDW1(K1,J2,L1)+C2*TBDW1(K2,J2,L1))+                 
     &     A2*B1*(C1*TBDW2(K1,J1,L1)+C2*TBDW2(K2,J1,L1))+                 
     &     A2*B2*(C1*TBDW2(K1,J2,L1)+C2*TBDW2(K2,J2,L1)))+  
     & D2*(A1*B1*(C1*TBDW1(K1,J1,L2)+C2*TBDW1(K2,J1,L2))+                 
     &     A1*B2*(C1*TBDW1(K1,J2,L2)+C2*TBDW1(K2,J2,L2))+                 
     &     A2*B1*(C1*TBDW2(K1,J1,L2)+C2*TBDW2(K2,J1,L2))+                 
     &     A2*B2*(C1*TBDW2(K1,J2,L2)+C2*TBDW2(K2,J2,L2)))  


      WIND=WW
      WIND=1.03*SQRT(UU**2 + VV**2) !SEE O:\amsr_l2\teff_climate\memo2.txt


      IF(UU.EQ.0 .AND. VV.EQ.0) THEN
         PHIW=0
      ELSE
         UU=-UU                 !CONVERT TO METEORLOGICAL SENSE (OUT OF)
         VV=-VV                 !CONVERT TO METEORLOGICAL SENSE (OUT OF)
         PHIW=ATAN2D(UU,VV)     !OUT OF NORTH = 0, OUT OF EAST= 90
         IF(PHIW.LT.0) PHIW=PHIW+360.
      ENDIF

      ! SWH
      IF(ABS(XLAT).GT.78.) then
         waveht = -999.0
      else
         BRIEF=78. - XLAT
         ILAT=NINT(1+BRIEF)
         IF (ILAT<1)   ILAT=1
         IF (ILAT>157) ILAT=157

         BRIEF=XLON/1.25
         ILON=NINT(1+BRIEF) 
         IF (ILON<1)   ILON=1 
         IF (ILON>288) ILON=288
         
         WAVEHT=A1*WAVEHT1(ILON,ILAT) + A2*WAVEHT2(ILON,ILAT)  
      endif

      ! Solar Flux
      solar_flux = a1 * solarflux1 + a2 * solarflux2


      call fd_land( irad, zang,sclon, fpt_lnd, frc_lnd, fland,gland)

      ilat=1 + int(2*(xlat+90))
      ilat = 360 - ilat + 1  ! fix lat reversal JMG 01/11/2011

      ilon=1 + int(2*xlon)
      if(ilat.eq.361) ilat=360
      if(ilon.eq.721) ilon=720

      FICE= A1*FICE1(ilon,ilat) + A2*FICE2(ilon,ilat)

      call fd_gice(a1,a2,fice1,fice2,xlat,xlon,irad, gice)

      call fd_gswe(a1,a2,swe1,swe2,xlat,xlon,irad, gswe)
      swe = gswe

!     sss/sst interpolation
      brief=4*(xlat+89.875)
      j1=int(1+brief)
      j2=j1+1
      b1=j1-brief
      b2=1.-b1          
      if(j1.eq.  0) j1=  1
      if(j2.eq.721) j2=720
 
      brief=4*(xlon-0.125)
      k1=int(1+brief)
      k2=k1+1
      c1=k1-brief
      c2=1.-c1
 
      wt(1)=c1*b1
      wt(2)=c2*b1
      wt(3)=c1*b2
      wt(4)=c2*b2


      sssa=-999; sssb=-999
      tsum=0; n=0
      do j=j1,j2
        do k=k1,k2
           kk=k
           if(kk.eq.   0) kk=1440
           if(kk.eq.1441) kk=   1
           n=n+1
           if(sss1(kk,j).eq.-999.) cycle
           tsum(0)=tsum(0) + wt(n)
           tsum(1)=tsum(1) + wt(n)*sss1(kk,j)
        enddo
      enddo
      if(tsum(0).ne.0) sssa=tsum(1)/tsum(0)

      tsum=0; n=0
      do j=j1,j2
        do k=k1,k2
           kk=k
           if(kk.eq.   0) kk=1440
           if(kk.eq.1441) kk=   1
           n=n+1
           if(sss2(kk,j).le.-999.) cycle
           tsum(0)=tsum(0) + wt(n)
           tsum(1)=tsum(1) + wt(n)*sss2(kk,j)
        enddo
      enddo
      if(tsum(0).ne.0) sssb=tsum(1)/tsum(0)

      if(sssa.lt.-998 .and. sssb.lt.-998) then
        sss_reference=35
        goto 3000
      endif

      if(sssa.lt.-998) sssa=sssb
      if(sssb.lt.-998) sssb=sssa
 
      sss_reference=a1*sssa + a2*sssb

 3000  continue
      ! surface temp
      surtempa=-999; surtempb=-999
      tsum=0; n=0
      do j=j1,j2
        do k=k1,k2
           kk=k
           if(kk.eq.   0) kk=1440
           if(kk.eq.1441) kk=   1
           n=n+1
           if(surtemp1(kk,j).lt.+ (273.15-9.98)) cycle
           tsum(0)=tsum(0) + wt(n)
           tsum(1)=tsum(1) + wt(n)*surtemp1(kk,j)
        enddo
      enddo
      if(tsum(0).ne.0) surtempa=tsum(1)/tsum(0)

      tsum=0; n=0
      do j=j1,j2
        do k=k1,k2
           kk=k
           if(kk.eq.   0) kk=1440
           if(kk.eq.1441) kk=   1
           n=n+1
           if(surtemp2(kk,j).lt.(273.15-9.98)) cycle
           tsum(0)=tsum(0) + wt(n)
           tsum(1)=tsum(1) + wt(n)*surtemp2(kk,j)
        enddo
      enddo
      if(tsum(0).ne.0) surtempb=tsum(1)/tsum(0)

      if(surtempa.lt.-998 .and. surtempb.lt.-998) then
        surtemp=-999
        goto 4000
      endif

      if(surtempa.lt.-998) surtempa=surtempb
      if(surtempb.lt.-998) surtempb=surtempa
 
      surtemp=a1*surtempa + a2*surtempb
 4000  continue



! TB TOA LAND CORR
       xday=idayjl-1+isecdy/86400.d0 

       if(lyear.eq.4*int(lyear/4)) then
          xmon=12.d0*xday/366.d0
       else
          xmon=12.d0*xday/365.d0
       endif

       call fd_tb_toa_land_corr(irad, xmon, zang, sclon, 
     1                          dtb, tb_landcorr_vpol, tb_landcorr_hpol)

       RETURN
       END

                     


      SUBROUTINE READ_DATA(filename,
     1                     SURTEMP,SURP,CWAT,SM,subtemp,
     2                     UU,VV,WW,SWE,TRAN,TBUP,TBDW,fice, waveht, sss, SOLAR_FLUX)
      USE HDF5               ! This module contains all necessary modules    
      IMPLICIT NONE

      CHARACTER(256) FILENAME
      REAL(4),TARGET :: SURTEMP(360*4,180*4),SURP(360,181),CWAT(360,181),
     1                  SM(360,181),UU(360,181),VV(360,181),WW(360,181),
     2                  SWE(720,361),subtemp(360,181)
      REAL(4),TARGET :: TRAN(360,181,91),TBUP(360,181,91),TBDW(360,181,91),
     1                  fice(360,181), waveht(288,157)
      REAL(4),target :: sss(360*4,180*4)
      REAL(4),target :: SOLAR_FLUX

      integer(4) ilat,ilon,klon
      real(4) dummy_lat(1440)

      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims, data_dims_temp
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims_ice, data_dims_gfs, data_dims_swh

      REAL(4), POINTER :: pt

      include "aquarius.fin"

      data_dims(1) = 360
      data_dims(2) = 181
      data_dims(3) = 91

      data_dims_gfs(1) = 720
      data_dims_gfs(2) = 361

      data_dims_temp(1) = 360*4
      data_dims_temp(2) = 180*4

      data_dims_ice(1) = 720
      data_dims_ice(2) = 360

      data_dims_swh(1) = 288
      data_dims_swh(2) = 157

      pt => SURTEMP(1,1)
      call dump_hdf5(filename, 'Surface Temp', data_dims_temp, pt)

!     Reflect about the equator
      do ilat=1,360
         dummy_lat=surtemp(:,ilat)
         surtemp(:,ilat) = surtemp(:,720+1-ilat)
         surtemp(:,720+1-ilat) = dummy_lat
      enddo

      pt => SURP(1,1)
      call dump_hdf5(filename, 'Surface Pressure', data_dims, pt)

      pt => CWAT(1,1)
      call dump_hdf5(filename, 'Cloud Water', data_dims, pt)

      pt => SM(1,1)
      call dump_hdf5(filename, 'Soil Moisture', data_dims, pt)

      pt => subtemp(1,1)
      call dump_hdf5(filename, 'Sub-Surface Temp', data_dims, pt)

      pt => UU(1,1)
      call dump_hdf5(filename, 'UU Wind', data_dims, pt)

      pt => VV(1,1)
      call dump_hdf5(filename, 'VV Wind', data_dims, pt)

      pt => SWE(1,1)
      call dump_hdf5(filename, 'Snow Water', data_dims_gfs, pt)

      pt => FICE(1,1)
      call dump_hdf5(filename, 'FracIce', data_dims_ice, pt)

      pt => waveht(1,1)
      call dump_hdf5(filename, 'Significant Wave Height', data_dims_swh, pt)

      pt => TRAN(1,1,1)
      call dump_hdf5(filename, 'Transmittance', data_dims, pt)

      pt => TBUP(1,1,1)
      call dump_hdf5(filename, 'Upwelling_TB', data_dims, pt)

      pt => TBDW(1,1,1)
      call dump_hdf5(filename, 'Downwelling_TB', data_dims, pt)


      WW=SQRT(UU**2 + VV**2)    !COMPUTE WIND SPEED FORM U AND V COMPONENTS

      pt => SSS(1,1)
      call dump_hdf5(filename, 'SSS_Ref', data_dims_temp, pt)
!     shift to greenwich reference
      do ilat=1,720
         dummy_lat=sss(:,ilat)
         do ilon=1,1440
            klon=ilon+720
            if(klon.gt.1440) klon=klon-1440
            sss(ilon,ilat) = dummy_lat(klon)
         enddo
      enddo

      data_dims(1) = 1
      pt => solar_flux
      call dump_hdf5(filename, 'Solar Flux', data_dims, pt)

      RETURN

      END



      SUBROUTINE FIND_MONTH_DAY(LYEAR,IDAYJL, IMON,IDAY)

        INTEGER*4 IDAYFX(12,0:1)
      DATA IDAYFX/1,32,60,91,121,152,182,213,244,274,305,335,            
     1            1,32,61,92,122,153,183,214,245,275,306,336/           
                                                  
      ILEAP=0
        IF(LYEAR.EQ.4*INT(LYEAR/4)) ILEAP=1

      DO 10 JMON=2,12
        IF(IDAYFX(JMON,ILEAP).GT.IDAYJL) THEN
      IMON=JMON-1
        GO TO 20
        ENDIF
   10 CONTINUE
      IMON=12
   20 CONTINUE

      IDAY=1+IDAYJL-IDAYFX(IMON,ILEAP)
      RETURN
        END


      subroutine fd_gice(a1,a2,fice1,fice2,bore_lat,bore_lon,irad, gicex)
      USE HDF5
      implicit none

      integer(4), parameter :: nlat=360, nlon=720
      real(8), parameter :: pi=3.14159265358979323846d0
      real(8), parameter :: rad=pi/180.d0
      real(8), parameter :: rp=6356.752d3  !earth      polar radius (meters)
      real(8), parameter :: re=6378.137d3 !earth equatorial radius (meters)
      real(8), parameter :: ffac=(rp/re)*(rp/re)
      integer(4), parameter :: n_rad=3 ! number of radiometers (i.e., horns)

      real(4),    intent(in)  :: a1, a2, fice1(nlon,nlat), fice2(nlon,nlat)
      integer(4), intent(in)  :: irad
      real(4),    intent(in)  :: bore_lat,bore_lon
      real(4),    intent(out) :: gicex

      integer(4) ilat,jlat,klat,ilon,jlon,klon,idir,i,ibin,istart
      real(4) d_beta,d_alpha,xlat,xlon,sinlat,coslat,tanlat,axel,
     1        chi,coschi,sinchi,d_s,xlatg,dtht_dthtg,cos10
      real(4) celb(3),xdot,frac_ice
      real(4) cell_area(nlat),cel(3,nlon,nlat)
      real(8) tsum(0:1)
      real(4) cosd, sind, tand, atand
      real(4) value_4

      REAL(4),TARGET :: gain_ice(0:150,n_rad)

      CHARACTER datadir*255, filename*255
      real(4), POINTER :: pt

      INTEGER(HSIZE_T), DIMENSION(2) :: gain_ice_dims

      data istart/1/

      include "aquarius.fin"


      if(istart.eq.1) then
         istart=0

         cos10=cosd(10.)
         d_beta =  180.d0/nlat
         d_alpha=  360.d0/nlon
 
         do klat=1,nlat
            xlatg=-90.d0 + (klat-0.5d0)*d_beta
            value_4 = ffac*tand(xlatg)
            xlat=atand(value_4) !convert from geodetic to geocentric

            sinlat=sind(xlat)
            coslat=cosd(xlat)
            tanlat=sinlat/coslat
            value_4 = re*tanlat/rp
            chi=atand(value_4)
            coschi=cosd(chi)
            sinchi=sind(chi)
            axel=cosd(chi)      !=rearth*sinbeta
            d_s=((re*coschi)**2 + (rp*sinchi)**2) * sqrt((re*sinchi)**2 + 
     1          (rp*coschi)**2)
            dtht_dthtg=(rp*coslat/re)**2 + (re*sinlat/rp)**2
            cell_area(klat)=(cosd(chi)/rp)*dtht_dthtg*d_s*d_beta*d_alpha*rad**2

            do klon=1,nlon
               xlon= (klon-0.5d0)*d_alpha
               cel(1,klon,klat)=cosd(xlon)*coslat
               cel(2,klon,klat)=sind(xlon)*coslat
               cel(3,klon,klat)=sinlat
            enddo               !klon

         enddo                  !klat

         call getenv('OCDATAROOT', datadir)
         pt => gain_ice(0,1)

         gain_ice_dims(1) = 151
         gain_ice_dims(2) = n_rad
         filename = trim(datadir)//'/aquarius/radiometer/gain_ice.h5'
         call dump_hdf5(filename, 'gain_ice', gain_ice_dims, pt)

      endif                     !if istart.eq.1

      if(bore_lat.gt.-40 .and. bore_lat.lt.30) then !no ice possible 
         gicex=0
         return
      endif    

      value_4 = ffac*tand(bore_lat)
      xlat=atand(value_4) !convert from geodetic to geocentric
      celb(1)=cosd(bore_lon)*cosd(xlat)
      celb(2)=sind(bore_lon)*cosd(xlat)
      celb(3)=sind(xlat)

      ilat=1 + int(2*(bore_lat+90))
      ilon=1 + int(2*bore_lon)
      if(ilat.eq.nlat+1) ilat=nlat
      if(ilon.eq.nlon+1) ilon=nlon
      
      tsum=0

      do jlat=ilat-20, ilat+20
         if(jlat.lt.1 .or. jlat.gt.nlat) cycle
         xlat=-90.25 + 0.5*jlat

         jlon=ilon

         xdot=dot_product(celb,cel(:,jlon,jlat))
         if(xdot.lt.cos10) cycle
         call find_ibin(xdot, ibin)

         xlon=-0.25 + 0.5*jlon

         klat = 360 - int(2*(xlat+90))
         klon=1 + int(2*xlon)
         if(klat.eq.361) klat=360
         if(klon.eq.721) klon=720
         frac_ice= A1*FICE1(klon,klat) + A2*FICE2(klon,klat)

         tsum(0) = tsum(0) + cell_area(jlat)*gain_ice(ibin,irad)
         tsum(1) = tsum(1) + cell_area(jlat)*gain_ice(ibin,irad)*frac_ice

         do idir=-1,1,2

            do i=1,nlon-1
               jlon=ilon + idir*i                         
               if(jlon.lt.   1) jlon=jlon+nlon
               if(jlon.gt.nlon) jlon=jlon-nlon

               xdot=dot_product(celb,cel(:,jlon,jlat))
               if(xdot.lt.cos10) exit
               call find_ibin(xdot, ibin)
      
               xlon=-0.25 + 0.5*jlon

               klat = 360 - int(2*(xlat+90))
               klon=1 + int(2*xlon)
               if(klat.eq.361) klat=360
               if(klon.eq.721) klon=720
               frac_ice= A1*FICE1(klon,klat) + A2*FICE2(klon,klat)

               tsum(0) = tsum(0) + cell_area(jlat)*gain_ice(ibin,irad)
               tsum(1) = tsum(1) + cell_area(jlat)*gain_ice(ibin,irad)*frac_ice

            enddo               !i
         enddo                  !idir

      enddo                     !jlat

      if (tsum(1).eq.0) then
         gicex=0
         return
      endif

      if(tsum(0).eq.0) then
         write(*,*) 'tsum(0).eq.0 in find_ice_frac_power'
         call exit(1)
      endif

      gicex=tsum(1)/tsum(0)

      return
      end


      subroutine fd_gswe(a1,a2,fswe1,fswe2,bore_lat,bore_lon,irad, gswex)
      USE HDF5
      implicit none

      integer(4), parameter :: nlat=361, nlon=720
      real(8), parameter :: pi=3.14159265358979323846d0
      real(8), parameter :: rad=pi/180.d0
      real(8), parameter :: rp=6356.752d3  !earth      polar radius (meters)
      real(8), parameter :: re=6378.137d3 !earth equatorial radius (meters)
      real(8), parameter :: ffac=(rp/re)*(rp/re)
      integer(4), parameter :: n_rad=3 ! number of radiometers (i.e., horns)

      real(4),    intent(in)  :: a1, a2, fswe1(nlon,nlat), fswe2(nlon,nlat)
      integer(4), intent(in)  :: irad
      real(4),    intent(in)  :: bore_lat,bore_lon
      real(4),    intent(out) :: gswex

      integer(4) ilat,jlat,klat,ilon,jlon,klon,idir,i,ibin,istart
      real(4) d_beta,d_alpha,xlat,xlon,sinlat,coslat,tanlat,axel,
     1        chi,coschi,sinchi,d_s,xlatg,dtht_dthtg,cos10
      real(4) celb(3),xdot,swe
      real(4) cell_area(nlat),cel(3,nlon,nlat)
      real(8) tsum(0:1)
      real(4) cosd, sind, tand, atand
      real(4) value_4

      REAL(4),TARGET :: gain_ice(0:150,n_rad)

      CHARACTER datadir*255, filename*255
      real(4), POINTER :: pt

      INTEGER(HSIZE_T), DIMENSION(2) :: gain_ice_dims

      data istart/1/

      include "aquarius.fin"


      if(istart.eq.1) then
         istart=0

         cos10=cosd(10.)
         d_beta =  180.d0/nlat
         d_alpha=  360.d0/nlon
 
         do klat=1,nlat
            xlatg=-90.d0 + (klat-0.5d0)*d_beta
            value_4 = ffac*tand(xlatg)
            xlat=atand(value_4) !convert from geodetic to geocentric

            sinlat=sind(xlat)
            coslat=cosd(xlat)
            tanlat=sinlat/coslat
            value_4 = re*tanlat/rp
            chi=atand(value_4)
            coschi=cosd(chi)
            sinchi=sind(chi)
            axel=cosd(chi)      !=rearth*sinbeta
            d_s=((re*coschi)**2 + (rp*sinchi)**2) * sqrt((re*sinchi)**2 + 
     1          (rp*coschi)**2)
            dtht_dthtg=(rp*coslat/re)**2 + (re*sinlat/rp)**2
            cell_area(klat)=(cosd(chi)/rp)*dtht_dthtg*d_s*d_beta*d_alpha*rad**2

            do klon=1,nlon
               xlon= (klon-0.5d0)*d_alpha
               cel(1,klon,klat)=cosd(xlon)*coslat
               cel(2,klon,klat)=sind(xlon)*coslat
               cel(3,klon,klat)=sinlat
            enddo               !klon

         enddo                  !klat

         call getenv('OCDATAROOT', datadir)
         pt => gain_ice(0,1)

         gain_ice_dims(1) = 151
         gain_ice_dims(2) = n_rad
         filename = trim(datadir)//'/aquarius/radiometer/gain_ice.h5'
         call dump_hdf5(filename, 'gain_ice', gain_ice_dims, pt)

      endif                     !if istart.eq.1

      if(bore_lat.gt.-40 .and. bore_lat.lt.30) then !no ice possible 
         gswex=0
         return
      endif    

      value_4 = ffac*tand(bore_lat)
      xlat=atand(value_4) !convert from geodetic to geocentric
      celb(1)=cosd(bore_lon)*cosd(xlat)
      celb(2)=sind(bore_lon)*cosd(xlat)
      celb(3)=sind(xlat)

      ilat=1 + int(2*(bore_lat+90))
      ilon=1 + int(2*bore_lon)
      if(ilat.eq.nlat+1) ilat=nlat
      if(ilon.eq.nlon+1) ilon=nlon
      
      tsum=0

      do jlat=ilat-20, ilat+20
         if(jlat.lt.1 .or. jlat.gt.nlat) cycle
         xlat=-90.25 + 0.5*jlat

         jlon=ilon

         xdot=dot_product(celb,cel(:,jlon,jlat))
         if(xdot.lt.cos10) cycle
         call find_ibin(xdot, ibin)

         xlon=-0.25 + 0.5*jlon

         klat = 360 - int(2*(xlat+90))
         klon=1 + int(2*xlon)
         if(klat.eq.362) klat=361
         if(klon.eq.721) klon=720
         swe= A1*FSWE1(klon,klat) + A2*FSWE2(klon,klat)

         tsum(0) = tsum(0) + cell_area(jlat)*gain_ice(ibin,irad)
         tsum(1) = tsum(1) + cell_area(jlat)*gain_ice(ibin,irad)*swe

         do idir=-1,1,2

            do i=1,nlon-1
               jlon=ilon + idir*i                         
               if(jlon.lt.   1) jlon=jlon+nlon
               if(jlon.gt.nlon) jlon=jlon-nlon

               xdot=dot_product(celb,cel(:,jlon,jlat))
               if(xdot.lt.cos10) exit
               call find_ibin(xdot, ibin)
      
               xlon=-0.25 + 0.5*jlon

               klat = 360 - int(2*(xlat+90))
               klon=1 + int(2*xlon)
               if(klat.eq.362) klat=361
               if(klon.eq.721) klon=720
               swe= A1*FSWE1(klon,klat) + A2*FSWE2(klon,klat)

               tsum(0) = tsum(0) + cell_area(jlat)*gain_ice(ibin,irad)
               tsum(1) = tsum(1) + cell_area(jlat)*gain_ice(ibin,irad)*swe

            enddo               !i
         enddo                  !idir

      enddo                     !jlat

      if (tsum(1).eq.0) then
         gswex=0
         return
      endif

      if(tsum(0).eq.0) then
         write(*,*) 'tsum(0).eq.0 in find_ice_frac_power'
         call exit(1)
      endif

      gswex=tsum(1)/tsum(0)

      return
      end



      subroutine find_ibin(xdot, ibin)
      implicit none

      real(4), intent(in) :: xdot
      integer(4), intent(out)  :: ibin
      real(4) angle
      real(4) acosd

      if(xdot.ge.1) then
         angle=0
      else 
         angle=acosd(xdot)
      endif
      ibin=nint(10*angle)
      if(ibin.gt.150) ibin=150
      return
      end


