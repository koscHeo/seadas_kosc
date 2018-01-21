! RFI Detection
! T. Meissner, June 2010
! adapted from C. Ruf and S. Misra, April 2010

        subroutine rfi_detection(n_cyc, rfi, gainOff, cellat, cellon, zang, s_acc)
        use l2gen_aquarius_module
        implicit none

        integer(4) n_cyc

        integer(4), dimension(l_ex_data)        :: wm_ex,wd_ex
        real(4), dimension(l_ex_data)           :: td_ex,tm_ex
        integer(4), dimension(l_data)           :: wm_ix,wd_ix
        real(4), dimension(l_data)              :: td_ix,tm_ix

        integer(1), dimension(l_data)           :: iflag_ix
        integer(1), dimension(l_ex_data)        :: iflag_ex, iflag_exx

        real(4), dimension(l_data)              :: sacc_ix
        real(4), dimension(l_ex_data)           :: sacc_ex



        integer(4)                              :: irad, kpol
        integer(4)                              :: ix, i, j, rem, icyc, isubcyc, isacc
        integer(4)                              :: latix, lonix

        integer(4)                              :: nn1, m1, m2
        real(4)                                 :: xsd, xsc1, xsc2, xsc

        real(4), dimension(l_ex_data)           :: tdir_ave, tcln_ave

        integer(4)                              :: iw, jw

        real(4), dimension(4)                   :: gain_vec !VPMH

        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: s_acc
        real(4), dimension(                        n_rad, max_cyc) :: cellat,cellon
        real(8), dimension(                               max_cyc) :: zang
        real(4)  cellon_360

        type(rfi_struct) :: rfi
        type(gain_off_struct) :: gainOff

        if (rfi%idata/=1) stop ' rfi detector parameters have not been read'
        rfi%iflag_rfi = 0 ! default

        do irad=1,n_rad
           do kpol=1,npol       !V/P/M/H

!                       Step 1: Arrange data into 1-dimensional array
              sacc_ix=missing_val
              wm_ix = missing_val
              wd_ix = missing_val
              tm_ix = missing_val
              td_ix = missing_val
              ix = 1
              do  icyc = 1,n_cyc
                 do isubcyc = 1,n_subcyc
                    do isacc = 1,n_Sacc
                                ! Taking latitude/longitude to determine index from which to get Algo parameters
                       cellon_360 = cellon(irad,icyc);
                       if (cellon_360.lt.0) cellon_360 = 360 + cellon_360
                       latix = floor(cellat(irad,icyc)) + 91
                       lonix = floor(cellon_360) + 1                          

                       if (latix<1)   latix=1
                       if (latix>181) latix=181
                       if (lonix<1)   lonix=1
                       if (lonix>360) lonix=360

                       sacc_ix(ix) = s_acc(isacc, isubcyc ,kpol, irad, icyc)
                       wm_ix(ix)   = rfi%wm(latix,lonix)
                       wd_ix(ix)   = rfi%wd(latix,lonix)
                                                        
                       gain_vec    = (/gainOff%gvv(irad,icyc), gainOff%gpp(irad,icyc), 
     &                   gainOff%gmm(irad,icyc), gainOff%ghh(irad,icyc)/)
                       if (abs(gain_vec(kpol)-missing_val)>0.1) then            
!     threshold for inclusion in "clean" mean
                          tm_ix(ix)   = rfi%tm(latix,lonix)*gain_vec(kpol)*rfi%stdta_rad(latix,lonix,kpol,irad) 
                          if (isacc==1 .or. isacc==2) tm_ix(ix)=tm_ix(ix)/sqrt(2.0) ! SA1 and SA2 have 2 obs
                                                        
!     threshold for RFI flag
                          if (zang(icyc).eq.-999) then
                             td_ix(ix)   = rfi%td(latix,lonix)*gain_vec(kpol)*rfi%stdta_rad(latix,lonix,kpol,irad) 
                          else if (zang(icyc).ge.0 .and. zang(icyc).lt.180) then
                             td_ix(ix)   = rfi%td(latix,lonix)*gain_vec(kpol)*rfi%stdta_rad(latix,lonix,kpol,irad)
                          else
                             td_ix(ix)   = rfi%td(latix,lonix)*gain_vec(kpol)*rfi%stdta_rad(latix,lonix,kpol,irad) 
                          endif

                          if (isacc==1 .or. isacc==2) td_ix(ix)=td_ix(ix)/sqrt(2.0) ! SA1 and SA2 have 2 obs
                       endif
                       ix=ix+1
                       
                    enddo       ! isacc
                 enddo          ! isubcyc
              enddo             ! icyc


                                ! Step 2: Expand all S_acc-1&2, to two separate but identical values 
                                !         to keep the same sample time between different samples
              sacc_ex=missing_val
              wm_ex = missing_val
              wd_ex = missing_val
              tm_ex = missing_val
              td_ex = missing_val
              tdir_ave = missing_val
              tcln_ave = missing_val
              j = 1
              do  i = 1,l_data
                 if (j> l_ex_data) stop ' error in j counting in RFI detector Step 2'
                 rem = mod(i,6) 
                 if (rem == 1 .or. rem == 2) then ! copy entry of original array into first 2 slots
                    sacc_ex(j) = sacc_ix(i)
                    wm_ex(j) = wm_ix(i)
                    wd_ex(j) = wd_ix(i)
                    tm_ex(j) = tm_ix(i)
                    td_ex(j) = td_ix(i)
                    j = j + 1
                    sacc_ex(j) = sacc_ix(i) 
                    wm_ex(j) = wm_ix(i)
                    wd_ex(j) = wd_ix(i)
                    tm_ex(j) = tm_ix(i)
                    td_ex(j) = td_ix(i)
                    j = j + 1
                 else if (rem==3 .or. rem==4 .or. rem==5) then ! copy entry of original array into slot 
                    sacc_ex(j) = sacc_ix(i)
                    wm_ex(j) = wm_ix(i)
                    wd_ex(j) = wd_ix(i)
                    tm_ex(j) = tm_ix(i)
                    td_ex(j) = td_ix(i)
                    j = j + 1
                 else if (rem==0) then ! skip, but increase counter
                    j = j + 5
                 end if
              enddo             ! i loop


                                ! Step 3: Initialize iflag to 0. 
                                ! flag will be set to 0
              iflag_ix = 0
              iflag_ex=0
              iflag_exx = 0
        
                        
                                ! Step 4: RFI flagging algorithm
              do i=1,l_ex_data 
                 rem = mod(i,12)
                 if (rem>=8) cycle ! observation that does not view antenna. skip                            
                 if (abs(sacc_ex(i)-missing_val)<0.1) cycle ! invalid observation
                 if (abs(tm_ex(i)-missing_val)<0.1)   cycle ! invalid gain/threshold                    
                 if (abs(td_ex(i)-missing_val)<0.1)   cycle ! invalid gain/threshold            
                 
                                ! iw window size
                 if (mod(wm_ex(i),2)==0) then
                    iw=wm_ex(i)/2
                 else
                    iw=(wm_ex(i)-1)/2
                 endif
                 if (iw<1) stop ' insufficient window WM'

                                ! jw window size
                 jw=wd_ex(i)
                 if (jw<1) stop ' insufficient window WD'

                                ! Step 4a: Calculate "dirty" mean. Exclude sample i. Exclude missing observations.
                 NN1=0                                          
                 XSD=0.0
                 do j=i-iw,i+iw,1
                    if (j<1 .or. j>l_ex_data) cycle
                    if (j==i) cycle ! do not include sample i
                    rem=mod(j,12)
                    if (rem>=8) cycle ! not an antenna observation
                    if (abs(sacc_ex(j)-missing_val)<0.1) cycle ! invalid
                    NN1=NN1+1
                    XSD=XSD+sacc_ex(j)
                 enddo
                 if (NN1 > 0) then
                    XSD = XSD/NN1
                 else
                    iflag_ex(i)=1 ! flag for RFI if no dirty mean can be calculated
                    XSD=missing_val
                    cycle
                 endif
                 tdir_ave(i)=xsd

                                ! Step 4b: Calculate "clean" mean.
                                ! Exclude samples that were previously flagged for RFI.
                                ! Exclude samples in window that differ by more than TM from dirty mean.
                                ! Before sample
                 M1=0                                           
                 XSC1=0.0
                 do j=i-iw,i-1
                    if (j<1 .or. j>l_ex_data) cycle
                    rem=mod(j,12)
                    if (rem>=8) cycle ! not an antenna observation
                    if (abs(sacc_ex(j)-missing_val)<0.1) cycle ! invalid
                    if (iflag_ex(j)==1) cycle !  previously flagged for RFI
                    if (abs(sacc_ex(j)-tdir_ave(i))>TM_ex(i)) cycle ! exceeds threshold and is therefore not included in clean mean
                    M1=M1+1
                    XSC1=XSC1+sacc_ex(j)
                 enddo
                                ! After sample
                 M2=0                                           
                 XSC2=0.0
                 do j=i+1,i+iw
                    if (j<1 .or. j>l_ex_data) cycle
                    rem=mod(j,12)
                    if (rem>=8) cycle ! not an antenna observation
                    if (abs(sacc_ex(j)-missing_val)<0.1) cycle ! invalid
                    if (iflag_ex(j)==1) cycle !  previously flagged for RFI
                    if (abs(sacc_ex(j)-tdir_ave(i))>TM_ex(i)) cycle ! exceeds threshold and is therefore not included in clean mean             
                    M2=M2+1             
                    XSC2=XSC2+sacc_ex(j)
                 enddo
                 if (M1>0.and. M2>0) then
                    XSC = XSC1/(2.0*M1) + XSC2/(2.0*M2)
                 else if (M1>0 .and. M2==0) then
                    XSC = XSC1/M1
                 else if (M1==0 .and. M2>0) then
                    XSC = XSC2/M2
                 else
                    iflag_ex(i)=1 ! flag if no clean mean can be calculated
                    XSC=missing_val
                 endif
                 tcln_ave(i)=xsc
                                
                                ! Step 4c: Flag sample for RFI if it differs by more than TD from clean mean
                 if (abs(tcln_ave(i)-missing_val)>0.1) then
                    if (iflag_ex(i)==0 .and. abs(sacc_ex(i)-tcln_ave(i))>TD_ex(i)) iflag_ex(i)=1
                 endif  
                 
              enddo             ! RFI flagging of sample I


              
                                ! Step 5: Flag 2*wd surrounding samples in window i-wd, i+wd    
                                ! I am restarting the loop over i, because if doing this flagging within the original loop and if 
                                ! the window wd is large then the algorithm would end up flagging everything after an RFI event 
              do i=1,l_ex_data
                 if (iflag_ex(i)==1) then
                    do j=i-jw,i+jw
                       if (j<1 .or. j>l_ex_data) cycle
                       iflag_exx(j)=1
                    enddo
                 endif
              enddo


                                ! Step 6: Decompress extended array into normal array
              j = 1
              do  i = 1,l_data
                 if (j> l_ex_data) stop ' error in j counting in RFI detector Step 2'
                 rem = mod(i,6) 
                 if (rem == 1 .or. rem == 2) then   
                    if (iflag_exx(j)==1) iflag_ix(i)=1
                    j = j + 1
                    if (iflag_exx(j)==1) iflag_ix(i)=1
                    j = j + 1
                 else if (rem==3 .or. rem==4 .or. rem==5) then 
                    if (iflag_exx(j)==1) iflag_ix(i)=1
                    j = j + 1
                 else if (rem==0) then 
                    if (iflag_exx(j)==1) iflag_ix(i)=1
                    j = j + 5
                 end if
              enddo             ! i loop


                                ! Step 7: Re-arrange flag back into multi-dimensional array
              ix = 1
              do icyc = 1,n_cyc
                 do isubcyc = 1,n_subcyc
                    do isacc = 1,n_sacc
                       rfi%iflag_rfi(isacc,isubcyc,kpol,irad,icyc) = iflag_ix(ix)
                       ix = ix+1
                    enddo       ! isacc         
                 enddo          ! isubcyc       
              enddo             ! icyc
              
              
           enddo                ! kpol
        enddo                   ! irad
        
        return
        end subroutine rfi_detection 




        subroutine rfi_detection_CND(n_cyc, rfi, gainOff, cellat, cellon, zang, s_acc)
        use l2gen_aquarius_module
        implicit none

        integer(4) n_cyc

        integer(4), dimension(l_data_c)         :: wm_ix,wd_ix
        real(4), dimension(l_data_c)            :: td_ix,tm_ix
        integer(1), dimension(l_data_c)         :: iflag_ix, iflag_jx
        real(4), dimension(l_data_c)            :: sacc_ix
        real(4)                                                         :: sacc5, sacc6


        integer(4)                                                      :: irad, kpol
        integer(4)                                                      :: ix, i, j, icyc, isubcyc
        integer(4)                                                      :: latix, lonix

        integer(4)                                                      :: nn1, m1, m2
        real(4)                                                         :: xsd, xsc1, xsc2, xsc

        real(4), dimension(l_data_c)            :: tdir_ave, tcln_ave
        
        integer(4)                                                      :: iw, jw

        real(4), dimension(4)                           :: gain_vec !VPMH

        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: s_acc
        real(4), dimension(                        n_rad, max_cyc) :: cellat,cellon
        real(8), dimension(                               max_cyc) :: zang
        real(4)  cellon_360

        type(rfi_struct) :: rfi
        type(gain_off_struct) :: gainOff

        if (rfi%idata/=1) stop ' rfi detector parameters have not been read'
        rfi%iflag_rfi_CND = 0 ! default


        do irad=1,n_rad
           do kpol=1,npol       !V/P/M/H
              
!       Step 1: Arrange data into 1-dimensional array
              sacc_ix=missing_val
              wm_ix = missing_val
              wd_ix = missing_val
              tm_ix = missing_val
              td_ix = missing_val
              tdir_ave=missing_val
              tcln_ave=missing_val
              ix = 1
              do  icyc = 1,n_cyc
                 do isubcyc = 1,n_subcyc
                                ! Taking latitude/longitude to determine index from which to get Algo parameters
                    cellon_360 = cellon(irad,icyc);
                    if (cellon_360.lt.0) cellon_360 = 360 + cellon_360
                    latix = floor(cellat(irad,icyc)) + 91
                    lonix = floor(cellon_360) + 1                                       
                    if (latix<1)   latix=1
                    if (latix>181) latix=181
                    if (lonix<1)   lonix=1
                    if (lonix>360) lonix=360
                    
                    sacc5 = s_acc(5, isubcyc ,kpol, irad, icyc)
                    sacc6 = s_acc(6, isubcyc ,kpol, irad, icyc)
                                                        
                    if (abs(sacc6-missing_val)>0.1 .and. abs(sacc5-missing_val)>0.1) then 
                       sacc_ix(ix) = sacc6-sacc5                                                        
                    endif
                    
                    
                    wm_ix(ix)   = rfi%wm_c(latix,lonix)
                    wd_ix(ix)   = rfi%wd_c(latix,lonix)
                    
                    gain_vec    = (/gainOff%gvv(irad,icyc), gainOff%gpp(irad,icyc), gainOff%gmm(irad,icyc), gainOff%ghh(irad,icyc)/)
                    if (abs(gain_vec(kpol)-missing_val)>0.1) then               
                       tm_ix(ix) = rfi%tm(latix,lonix)*gain_vec(kpol)*
     &                      sqrt(rfi%stdta_rad(latix,lonix,kpol,irad)**2+rfi%stdta_cnd(kpol,irad)**2) 
                                ! threshold for inclusion in "clean" mean. 
                                ! difference between sacc6 and sacc5, which are assumed to be indepenent -> Gaussian error formula
                       
                       if (zang(icyc).eq.-999) then
                          td_ix(ix)   = rfi%td(latix,lonix)*gain_vec(kpol)*
     &                         sqrt(rfi%stdta_rad(latix,lonix,kpol,irad)**2+rfi%stdta_cnd(kpol,irad)**2) 
                       else if (zang(icyc).ge.0 .and. zang(icyc).lt.180) then
                          td_ix(ix)   = rfi%td(latix,lonix)*gain_vec(kpol)*
     &                         sqrt(rfi%stdta_rad(latix,lonix,kpol,irad)**2+rfi%stdta_cnd(kpol,irad)**2) 
                       else
                          td_ix(ix)   = rfi%td(latix,lonix)*gain_vec(kpol)*
     &                         sqrt(rfi%stdta_rad(latix,lonix,kpol,irad)**2+rfi%stdta_cnd(kpol,irad)**2) 
                       endif
                                ! threshold for RFI flag. 

                    endif
                    ix=ix+1
                 enddo          ! isubcyc
              enddo             ! icyc



                                ! Step 2: Initialize iflag to 0. 
                                ! flag will be set to 0
              iflag_ix = 0
              iflag_jx = 0
        
                        
                                ! Step 3: RFI flagging algorithm
              do i=1,l_data_c 
                 if (abs(sacc_ix(i)-missing_val)<0.1) cycle ! invalid observation
                 if (abs(tm_ix(i)-missing_val)<0.1)   cycle ! invalid gain/threshold                    
                 if (abs(td_ix(i)-missing_val)<0.1)   cycle ! invalid gain/threshold            

                                ! iw window size
                 if (mod(wm_ix(i),2)==0) then
                    iw=wm_ix(i)/2
                 else
                    iw=(wm_ix(i)-1)/2
                 endif
                 if (iw<1) stop ' insufficient window WM'

                                ! jw window size
                 jw=wd_ix(i)
                 if (jw<1) stop ' insufficient window WD'

                                ! Step 3a: Calculate "dirty" mean. Exclude sample i. Exclude missing observations.
                 NN1=0                                          
                 XSD=0.0
                 do j=i-iw,i+iw,1
                    if (j<1 .or. j>l_data_c) cycle
                    if (j==i) cycle ! do not include sample i
                    if (abs(sacc_ix(j)-missing_val)<0.1) cycle ! invalid
                    NN1=NN1+1
                    XSD=XSD+sacc_ix(j)
                 enddo
                 if (NN1 > 0) then
                    XSD = XSD/NN1
                 else
                    iflag_ix(i)=1 ! flag for RFI if no dirty mean can be calculated
                    XSD=missing_val
                    cycle
                 endif
                 tdir_ave(i)=xsd

                                ! Step 3b: Calculate "clean" mean.
                                ! Exclude samples that were previously flagged for RFI.
                                ! Exclude samples in window that differ by more than TM_c from dirty mean.
                                ! Before sample
                 M1=0                                           
                 XSC1=0.0
                 do j=i-iw,i-1
                    if (j<1 .or. j>l_data_c) cycle
                    if (abs(sacc_ix(j)-missing_val)<0.1) cycle ! invalid
                    if (iflag_ix(j)==1) cycle !  previously flagged for RFI
                    if (abs(sacc_ix(j)-tdir_ave(i))>TM_ix(i)) cycle ! exceeds threshold and is therefore not included in clean mean
                    M1=M1+1
                    XSC1=XSC1+sacc_ix(j)
                 enddo
                                ! After sample
                 M2=0                                           
                 XSC2=0.0
                 do j=i+1,i+iw
                    if (j<1 .or. j>l_data_c) cycle
                    if (abs(sacc_ix(j)-missing_val)<0.1) cycle ! invalid
                    if (iflag_ix(j)==1) cycle !  previously flagged for RFI
                    if (abs(sacc_ix(j)-tdir_ave(i))>TM_ix(i)) cycle ! exceeds threshold and is therefore not included in clean mean             
                    M2=M2+1             
                    XSC2=XSC2+sacc_ix(j)
                 enddo
                 if (M1>0.and. M2>0) then
                    XSC = XSC1/(2.0*M1) + XSC2/(2.0*M2)
                 else if (M1>0 .and. M2==0) then
                    XSC = XSC1/M1
                 else if (M1==0 .and. M2>0) then
                    XSC = XSC2/M2
                 else
                    iflag_ix(i)=1 ! flag if no clean mean can be calculated
                    XSC=missing_val
                 endif
                 tcln_ave(i)=xsc
                                
                                ! Step 3c: Flag sample for RFI if it differs by more than TD_c from clean mean
                 if (abs(tcln_ave(i)-missing_val)>0.1) then
                    if (iflag_ix(i)==0 .and. abs(sacc_ix(i)-tcln_ave(i))>TD_ix(i)) iflag_ix(i)=1
                 endif  

              enddo             ! RFI flagging of sample I

                
                                ! Step 4: Flag 2*wd_c surrounding samples in window i-wd_c, i+wd_c      
                                ! I am restarting the loop over i, because if doing this flagging within the original loop and if 
                                ! the window wd_c is large then the algorithm would end up flagging everything after an RFI event. 
              do i=1,l_data_c   
                 if (iflag_ix(i)==1) then
                    do j=i-jw,i+jw
                       if (j<1 .or. j>l_data_c) cycle
                       iflag_jx(j)=1
                    enddo
                 endif
              enddo     


                                ! Step 5: Re-arrange flag back into multi-dimensional array
              ix = 1
              do icyc = 1,n_cyc
                 do isubcyc = 1,n_subcyc
                    rfi%iflag_rfi_cnd(isubcyc,kpol,irad,icyc) = iflag_jx(ix)
                    ix = ix+1
                 enddo          ! isubcyc       
              enddo             ! icyc
              

           enddo                ! kpol
        enddo                   ! irad
                

        return
        end subroutine rfi_detection_CND 
