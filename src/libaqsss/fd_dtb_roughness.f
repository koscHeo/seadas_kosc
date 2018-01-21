        subroutine fd_dtb_roughness (itype, iwdir, irad, wspd, phir, sigv, sigh, 
     1                               acoef, bcoef, arr_dtbw, dsigma, iflag_dtbw, 
     2                               ioob, dtb_rough)
!      use constants_module
!      use array_module 
        
        ! wind induced emissivity from wind speed + scatterometer sigma0
        implicit none

        integer(4), intent(in)                                  ::  itype 
        ! 0=wspd only  1=wspd+sigma0_VV  2=wspd+sigma0_HH   3=sigma0_VV+sigma0_HH  
        
        integer(4), intent(in)                                  ::  iwdir 
      ! 0=no wdir correction is performed   1= wdir correction is performed

        integer(4), parameter :: n_rad=3
        integer(4), parameter :: nbin_w=60

        integer(4), intent(in)                                  ::  irad
        real(4), intent(in)                                             ::  wspd
        real(4), intent(in)                                             ::  phir
        real(4), intent(in)                                             ::  sigv, sigh  ! sigma0 in real units before wdir signal removed
         real(8), intent(in)                           :: acoef(0:2, 2, 3, 5) 
         real(8), intent(in)                           :: bcoef(0:2, 4, 3, 5) 

        real(8),  intent(in)                           :: arr_dtbw(3,2,n_rad,nbin_w,nbin_w)
        real(4),  intent(in)                           :: dsigma(2,n_rad)
        integer(1), intent(in)                         :: iflag_dtbw(3,2,n_rad,nbin_w,nbin_w)


        integer(4), parameter :: npol_sc=4

        integer(4), intent(out)                                 ::  ioob ! OOB flag  
                                                           ! 0 = can use retrieval    1=underpopulated bin. do not use
        real(4), dimension(2), intent(out)          ::  dtb_rough  !(1=V, 2=H) 
        
        real(4), dimension(0:2,2)                             ::  aharm !1=V,  2=H
        real(4), dimension(0:2,npol_sc)               ::  bharm !1=VV, 2=HH,  3=VH,  4=HV  
 
        real(4) :: dwin
        integer(4)  :: mbin_w 

        data mbin_w/60/
        data dwin/1.0000/

        real(4),    parameter :: missing_val=-9999. ! indicator for missing value

        integer(4)                                                              ::  ibinx1, ibiny1
        integer(4)                                                              ::  ibinx2, ibiny2
        real(4)                                                                 ::  xpoint, ypoint
        real(4)                                                                 ::  xvec0, yvec0
        real(4)                                                                 ::  dx, dy
        real(4)                                                                 ::  x1, y1
        real(4)                                                                 ::  x2, y2
        real(4)                                                                 ::  brief_x, brief_y

        real(4), dimension(2)                                       ::  c1, c2, z1, z2


        real(4), dimension(2)                                       ::  xtbw_scat  !(1=V, 2=H) 

        real(4)                                                                 ::  xsigv, xsigh        ! sigma0 in real units after wdir signal removed
        real(4), dimension(2)                                       ::  dtb_phi, dsigma_phi

        real(4) cosd


      if (mbin_w /= nbin_w) then
                write(*,*) mbin_w, nbin_w
                write(*,*) ' bin mismatch in subroutine initialize_roughness_correction.'
                stop
        endif  

        if (iwdir==1) then
                call fd_TM_emiss_harmonics(irad,wspd, acoef, aharm)
                dtb_phi(1:2)    = aharm(1,1:2)*cosd(phir) + aharm(2,1:2)*cosd(2.0*phir)         
                call fd_TM_scat_harmonics(irad, wspd, bcoef, bharm) 
                dsigma_phi(1:2) = bharm(1,1:2)*cosd(phir) + bharm(2,1:2)*cosd(2.0*phir)
        else if (iwdir==0) then ! isotropic
                dtb_phi=0.0
                dsigma_phi=0.0
        endif


        xsigv = sigv - dsigma_phi(1)   
        xsigh = sigh - dsigma_phi(2)

        ioob=0


        if (itype==0) then  ! use wind speed only
                dtb_rough(:) = aharm(0,:) + dtb_phi(:)
                return
        endif


        if (itype==1) then
                xpoint = wspd
                ypoint = xsigv
                xvec0  = dwin/2.0  
                yvec0  = dsigma(1,irad)/2.0
                dx     = dwin  
                dy     = dsigma(1,irad)
        else if (itype==2) then
                xpoint = wspd
                ypoint = xsigh
                xvec0  = dwin/2.0  
                yvec0  = dsigma(2,irad)/2.0
                dx     = dwin 
                dy     = dsigma(2,irad)
        else if (itype==3) then
                xpoint = xsigv
                ypoint = xsigh
                xvec0  = dsigma(1,irad)/2.0  
                yvec0  = dsigma(2,irad)/2.0
                dx     = dsigma(1,irad)  
                dy     = dsigma(2,irad)
        endif

        
        ibinx1 = floor((xpoint-xvec0)/dx) + 1
        if (ibinx1<1)        ibinx1=1
        if (ibinx1>nbin_w-1) ibinx1=nbin_w-1    
        ibinx2=ibinx1+1
        x1 = xvec0 + (ibinx1-1)*dx
        x2 = x1 + dx
        brief_x = (xpoint-x1)/(x2-x1)

        ibiny1 = floor((ypoint-yvec0)/dy) + 1
        if (ibiny1<1)        ibiny1=1
        if (ibiny1>nbin_w-1) ibiny1=nbin_w-1    
        ibiny2=ibiny1+1
        y1 = yvec0 + (ibiny1-1)*dy
        y2 = y1 + dy
        brief_y = (ypoint-y1)/(y2-y1)

!        if (irad.eq.1) then
!           write(*,*) xpoint,ypoint,sigv,sigh,dsigma_phi
!        endif

        if (iflag_dtbw(itype,1,irad,ibinx1,ibiny1)==0 .and. iflag_dtbw(itype,1,irad,ibinx2,ibiny1)==0) then
                z1 = arr_dtbw(itype, :, irad,ibinx1,ibiny1)
                z2 = arr_dtbw(itype, :, irad,ibinx2,ibiny1)
        else if (iflag_dtbw(itype,1,irad,ibinx1,ibiny1)==0 .and. iflag_dtbw(itype,1,irad,ibinx2,ibiny1)==1) then
                z1 = arr_dtbw(itype, :, irad,ibinx1,ibiny1)
                z2 = z1
        else if (iflag_dtbw(itype,1,irad,ibinx1,ibiny1)==1 .and. iflag_dtbw(itype,1,irad,ibinx2,ibiny1)==0) then
                z1 = arr_dtbw(itype, :, irad,ibinx2,ibiny1)
                z2 = z1
        else
                ioob=1
                dtb_rough=missing_val
!                write(*,*) ibinx1,ibiny1,iflag_dtbw(itype,1,irad,ibinx1,ibiny1),iflag_dtbw(itype,1,irad,ibinx2,ibiny1)
                return
        endif
        c1 = z1*(1.0-brief_x) + z2*brief_x

        if (iflag_dtbw(itype,1,irad,ibinx1,ibiny2)==0 .and. iflag_dtbw(itype,1,irad,ibinx2,ibiny2)==0) then
                z1 = arr_dtbw(itype, :, irad,ibinx1,ibiny2)
                z2 = arr_dtbw(itype, :, irad,ibinx2,ibiny2)
        else if (iflag_dtbw(itype,1,irad,ibinx1,ibiny2)==0 .and. iflag_dtbw(itype,1,irad,ibinx2,ibiny2)==1) then
                z1 = arr_dtbw(itype, :, irad,ibinx1,ibiny2)
                z2 = z1
        else if (iflag_dtbw(itype,1,irad,ibinx1,ibiny2)==1 .and. iflag_dtbw(itype,1,irad,ibinx2,ibiny2)==0) then
                z1 = arr_dtbw(itype, :, irad,ibinx2,ibiny2)
                z2 = z1
        else
                ioob=1
                dtb_rough=missing_val
!                write(*,*) iflag_dtbw(itype,1,irad,ibinx1,ibiny1),iflag_dtbw(itype,1,irad,ibinx2,ibiny1)
                return
        endif
        c2 = z1*(1.0-brief_x) + z2*brief_x

        xtbw_scat = c1*(1.0-brief_y) + c2*brief_y
        if (xtbw_scat(1)<0.0) xtbw_scat(1)=0.0
        if (xtbw_scat(2)<0.0) xtbw_scat(2)=0.0

        dtb_rough = xtbw_scat + dtb_phi


        return
        end subroutine fd_dtb_roughness
