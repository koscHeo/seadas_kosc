!     AQ_calibration module is a code package that contains all routines for turning raw Aquarius radiometer counts 
!     into apparent antenna temperatures.
!     The code is based on the following documents:
!     [1]: J. Piepmeyer, Memorandum from 02/25/2009: Radiometer TAA calibration algorithm for V & H polarizations AQU-MEMO-0013xx
!     [2]: J. Piepmeyer, Memorandum from 03/24/2009: How to calibrate for T_i,U
!     [3]: J. Piepmeyer, Memorandum from 02/20/2009: Radiometer non-linearity correction. 

!    Updates: May 24 2010
!    1. Updated ATBD Algorithm Theoretical Basis Document for T_AA by J. Piepmeyer, April 28, 2010
!    2. Gain Glitch Detector, by J. peng + J. Piepmeyer,  May 2010
!    3. Spreadsheet interactiveversion1.xlsm by J. Piepmeyer, May 14, 2010
!    4. RFI detection code by C. Ruf and S. Misra, April 30, 2010
!    5. Change "stop" to "write(*,*)" for zero counts by J. Gales, Aug 5, 2011

        subroutine count_to_TA(coeff_loss_file, coeff_nl_file,
     1                         n_cyc, t_prt, iTemp, gainOff, rfi, 
     2                         S_acc_raw, L_acc_raw, S_acc, L_acc, 
     3                         cellat, cellon, zang, c_deltaTND, TA_hat_drift_corr,
     4                         TA_hat, TA, TF_hat, TF)
        use l2gen_aquarius_module
        implicit none

        character*256 coeff_loss_file
        character*256 coeff_nl_file

        integer(4) icyc
        integer(4) n_cyc

        real(4), dimension(n_prt,max_cyc) :: t_prt

        type(instrument_temp_struct) :: iTemp
        type(gain_off_struct) :: gainOff
        type(rfi_struct) :: rfi
 
        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc_raw
        real(4), dimension(n_Lacc,           npol, n_rad, max_cyc) :: L_acc_raw
        real(4),    dimension(  n_rad, max_cyc) ::  cellat,cellon
        real(8),    dimension(  n_rad, max_cyc) ::  zang
        real(4),    dimension(  2, n_rad) ::  c_deltaTND
        real(4), dimension( 2, n_rad, max_cyc) :: TA_hat_drift_corr

        real(4), dimension(3, n_rad, max_cyc)  :: TF_hat, TF, TA_hat, TA

        type(prelaunch_struct) :: plc

        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc
        real(4), dimension(n_Lacc,           npol, n_rad, max_cyc) :: L_acc
        integer(4), dimension(               npol, n_rad, max_cyc) :: m_gain, m_off ! dynamical values for gain/offset 

        integer(4) :: NTOT
        integer(4) irad,ipol

#if 0
       open(unit=8,file='s_acc_raw.dat', form='UNFORMATTED')
       write(8) n_cyc
       write(8) s_acc_raw
       close(8)

       open(unit=8,file='l_acc_raw.dat', form='UNFORMATTED')
       write(8) n_cyc
       write(8) l_acc_raw
       close(8)
#endif

!    Step 0: Get prelaunch calibration coeffs
        call get_prelaunch_calibration_coeffs(coeff_loss_file, coeff_nl_file, plc)

!     Step 1: Compute Zone temperatures from PRT values
        call get_zone_temperatures(n_cyc, t_prt, iTemp)

!     Step 2: non-linearity correction
        call nl_correction(n_cyc, iTemp%Tdet, plc, S_acc_raw, L_acc_raw, S_acc, L_acc)

#if 0
       open(unit=8,file='s_acc.dat', form='UNFORMATTED')
       write(8) s_acc
       close(8)

       open(unit=8,file='l_acc.dat', form='UNFORMATTED')
       write(8) l_acc
       close(8)
#endif

!     Step 3: gain glitch detector
        call gain_glitch_detector(n_cyc, rfi, L_acc, m_gain, m_off, ntot)

!     Step 4: Determine Dicke load reference temperature and CND and ND injection temperatures + offset temperatures
        do icyc =1,n_cyc
        iTemp%TCND(:,:,icyc) = plc%TCND_0(:,:) + plc%dTCND_dT4(:,:)*(iTemp%T4(:,:,icyc) - plc%Tref_CND) + 
     1                         plc%dTCND_dT(:,:)*(iTemp%TCND_P(:,:,icyc) - plc%Tref_CND)  
     
        iTemp%TDL_offset(:,:,icyc) = plc%TDL_offset_0(:,:) + plc%dTDL_offset_dT(:,:) *
     1                               (iTemp%TDL(:,:,icyc) - plc%Tref_DL_offset)

        iTemp%T0(:,:,icyc) = iTemp%TDL(:,:,icyc) + iTemp%TDL_offset(:,:,icyc)

        iTemp%TND(:,:,icyc) = plc%TND_0(:,:) + plc%dTND_dT*(iTemp%TND_P(:,:,icyc) -  plc%Tref_ND)
        iTemp%TND(:,:,icyc) = iTemp%TND(:,:,icyc) * (1 - c_deltaTND)

        iTemp%TND_offset(:,:,icyc) = plc%TND_offset_0(:,:) + plc%DTND_offset_dT(:,:) *
     1                               (iTemp%TND_P(:,:,icyc) - plc%Tref_ND_offset)
 
        enddo

!       open(unit=8,file='T0.dat', form='UNFORMATTED'); write(8) iTemp%T0; close(8)
!       open(unit=8,file='TCND.dat', form='UNFORMATTED'); write(8) iTemp%TCND; close(8)
!       open(unit=8,file='TDL.dat', form='UNFORMATTED'); write(8) iTemp%TDL; close(8)
!       open(unit=8,file='TND.dat', form='UNFORMATTED'); write(8) iTemp%TND; close(8)
!       open(unit=8,file='TND_offset.dat', form='UNFORMATTED'); write(8) iTemp%TND_offset; close(8)
!       open(unit=8,file='TDL_offset.dat', form='UNFORMATTED'); write(8) iTemp%TDL_offset; close(8)

!     Step 5: gains and offsets
        call find_cal_VH(n_cyc, m_gain, m_off, iTemp, gainOff%gvv, gainOff%ghh, gainOff%ov, gainOff%oh, L_acc)

        call find_cal_PM(n_cyc, m_gain, m_off, iTemp%T0, iTemp%TND, iTemp%TDL_offset, iTemp%TND_offset, 
     1                          gainOff%gpv, gainOff%gph, gainOff%op, gainOff%om, 
     2                          gainOff%gmv, gainOff%gmh, gainOff%gpp, gainOff%gmm, L_acc)


!     Step 6: RFI detector
        call get_rfi_auxiliary_data(rfi)

!       Sets rfi%iflag_rfi
        call rfi_detection(n_cyc, rfi, gainOff, cellat, cellon, zang, s_acc)

!       Sets rfi%iflag_rfi_cnd
        call rfi_detection_CND(n_cyc, rfi, gainOff, cellat, zang, cellon, s_acc)

!     Step 7: calculate gains for U channel

        call find_cal_hU(n_cyc, plc%CL5, plc%CLMM, plc%dL5_dT5, itemp%T5, iTemp%TCND, plc%Tref_5, gainOff, S_acc) ! without RFI filter 
        call find_cal_gU(n_cyc, plc%CL5, plc%CLMM, plc%dL5_dT5, itemp%T5, iTemp%TCND, plc%Tref_5, gainOff, rfi, S_acc) ! with RFI filter

!       open(unit=8,file='go.dat', form='UNFORMATTED')
!       write(8) gainOff%gvv, gainOff%ghh, gainOff%ov, gainOff%oh, gainOff%op, gainOff%om, gainOff%gpv, gainOff%gph, 
!     1           gainOff%gmv, gainOff%gmh, gainOff%hpU, gainOff%hmU, gainOff%gpU, gainOff%gmU
!       close(8)

!     Step 8: Count to TA before front end loss corrections. averageing over subcycles
!             Sets rfi%num_rfi
        call find_TA_hat(n_cyc, gainOff, rfi, TA_hat, TF_hat, S_acc)

!     Apply TA_hat drift correction (if not applicable TA_hat_drift_corr values have been set to 0)
        do icyc=1,n_cyc
           do irad=1,n_rad
              do ipol=1,2
                 if (abs(TA_hat(ipol,irad,icyc)-missing_val)<0.1) cycle

! Remove for run 5  JMG  11/12/14
!                 TA_hat_drift_corr(ipol, irad, icyc) = TA_hat_drift_corr(ipol, irad, icyc) * iTemp%TND(ipol, irad, icyc)
                 TA_hat(ipol, irad, icyc) = 
     1                TA_hat(ipol, irad, icyc) + TA_hat_drift_corr(ipol, irad, icyc)
              enddo
           enddo
        enddo

        if (c_deltaTND(1,1) .ne. 0.0) then
#if 0
           open(unit=8,file='rfi.dat', form='UNFORMATTED')
           write(8) rfi%iflag_rfi, rfi%iflag_rfi_cnd, rfi%num_rfi
           close(8)

           open(unit=9,file='rfi_glitch.dat', form='UNFORMATTED')
           write(9) rfi%iflag_glitch
           close(9)
#endif
        endif


!     Step 9: Front end loss corrections
        call fe_loss_corr(n_cyc, plc, iTemp, TA_hat, TF_hat, TF, TA)

!       open(unit=8,file='ta.dat', form='UNFORMATTED')
!       write(8) ta_hat, ta
!       close(8)

        return
        end subroutine count_to_TA



!     Step1: performs mapping between PRT and temperature in the various zones which are fed into fthe forward model
        subroutine get_zone_temperatures(n_cyc, t_prt, iTemp)
        use l2gen_aquarius_module
        implicit none

        integer(4) icyc
        integer(4) n_cyc

        real(4),    dimension(n_prt,max_cyc) :: t_prt

        type(instrument_temp_struct) :: iTemp

        do icyc=1,n_cyc

!     Reflector temperature (zone 1),  average over 8 thermistors,  same for each radiometer 
        iTemp%T1(1,icyc) = sum(t_prt(34:41,icyc))/8.0 + 273.15
        iTemp%T1(2,icyc) = sum(t_prt(34:41,icyc))/8.0 + 273.15
        iTemp%T1(3,icyc) = sum(t_prt(34:41,icyc))/8.0 + 273.15

!     zone 2A
        iTemp%T2A(1,icyc) = t_prt(1,icyc) + 273.15
        iTemp%T2A(2,icyc) = t_prt(2,icyc) + 273.15
        iTemp%T2A(3,icyc) = t_prt(3,icyc) + 273.15

!       zone 3
        iTemp%T3(1,1,icyc) = t_prt(5,icyc) + 273.15 !RAD1 V
        iTemp%T3(2,1,icyc) = t_prt(4,icyc) + 273.15 !RAD1 H
        iTemp%T3(1,2,icyc) = t_prt(7,icyc) + 273.15 !RAD2 V
        iTemp%T3(2,2,icyc) = t_prt(6,icyc) + 273.15 !RAD2 H
        iTemp%T3(1,3,icyc) = t_prt(9,icyc) + 273.15 !RAD3 V
        iTemp%T3(2,3,icyc) = t_prt(8,icyc) + 273.15 !RAD3 H

!     zone 2B = arithmetic average of 2A and 2B ?
        iTemp%T2B(1,icyc) = (iTemp%T2A(1,icyc) + (iTemp%T3(1,1,icyc)+iTemp%T3(2,1,icyc))/2.0) /2.0
        iTemp%T2B(2,icyc) = (iTemp%T2A(2,icyc) + (iTemp%T3(1,2,icyc)+iTemp%T3(2,2,icyc))/2.0) /2.0
        iTemp%T2B(3,icyc) = (iTemp%T2A(3,icyc) + (iTemp%T3(1,3,icyc)+iTemp%T3(2,3,icyc))/2.0) /2.0

!     zone 4
        iTemp%T4(1,1,icyc) = t_prt(11,icyc) + 273.15 ! RAD1 V
        iTemp%T4(2,1,icyc) = t_prt(10,icyc) + 273.15 ! RAD1 H
        iTemp%T4(1,2,icyc) = t_prt(13,icyc) + 273.15 ! RAD2 V
        iTemp%T4(2,2,icyc) = t_prt(12,icyc) + 273.15 ! RAD2 H
        iTemp%T4(1,3,icyc) = t_prt(15,icyc) + 273.15 ! RAD3 V
        iTemp%T4(2,3,icyc) = t_prt(14,icyc) + 273.15 ! RAD3 H

 !    zone 5
        iTemp%T5(1,1,icyc) = t_prt(17,icyc) + 273.15 ! RAD1 V
        iTemp%T5(2,1,icyc) = t_prt(16,icyc) + 273.15 ! RAD1 H
        iTemp%T5(1,2,icyc) = t_prt(19,icyc) + 273.15 ! RAD2 V
        iTemp%T5(2,2,icyc) = t_prt(18,icyc) + 273.15 ! RAD2 H
        iTemp%T5(1,3,icyc) = t_prt(21,icyc) + 273.15 ! RAD3 V
        iTemp%T5(2,3,icyc) = t_prt(20,icyc) + 273.15 ! RAD3 H

!     Dicke
        iTemp%TDL(1,1,icyc) = (t_prt(24,icyc) + t_prt(25,icyc))/2.0 + 273.15 ! RAD1 V
        iTemp%TDL(2,1,icyc) = (t_prt(22,icyc) + t_prt(23,icyc))/2.0 + 273.15 ! RAD1 H
        iTemp%TDL(1,2,icyc) = (t_prt(28,icyc) + t_prt(29,icyc))/2.0 + 273.15 ! RAD2 V
        iTemp%TDL(2,2,icyc) = (t_prt(26,icyc) + t_prt(27,icyc))/2.0 + 273.15 ! RAD2 H
        iTemp%TDL(1,3,icyc) = (t_prt(32,icyc) + t_prt(33,icyc))/2.0 + 273.15 ! RAD1 V
        iTemp%TDL(2,3,icyc) = (t_prt(30,icyc) + t_prt(31,icyc))/2.0 + 273.15 ! RAD1 H

!     ND 
        iTemp%TND_P(1,1,icyc) = (t_prt(72,icyc) + t_prt(74,icyc))/2.0 + 273.15 ! RAD1 V
        iTemp%TND_P(2,1,icyc) = (t_prt(73,icyc) + t_prt(75,icyc))/2.0 + 273.15 ! RAD1 H
        iTemp%TND_P(1,2,icyc) = (t_prt(77,icyc) + t_prt(79,icyc))/2.0 + 273.15 ! RAD2 V
        iTemp%TND_P(2,2,icyc) = (t_prt(79,icyc) + t_prt(80,icyc))/2.0 + 273.15 ! RAD2 H
        iTemp%TND_P(1,3,icyc) = (t_prt(82,icyc) + t_prt(84,icyc))/2.0 + 273.15 ! RAD3 V
        iTemp%TND_P(2,3,icyc) = (t_prt(83,icyc) + t_prt(85,icyc))/2.0 + 273.15 ! RAD3 H

!     CND 
        iTemp%TCND_P(1,1,icyc) = t_prt(71,icyc) + 273.15 !RAD1
        iTemp%TCND_P(2,1,icyc) = t_prt(71,icyc) + 273.15 !RAD1
        iTemp%TCND_P(1,2,icyc) = t_prt(76,icyc) + 273.15 !RAD2
        iTemp%TCND_P(2,2,icyc) = t_prt(76,icyc) + 273.15 !RAD2
        iTemp%TCND_P(1,3,icyc) = t_prt(81,icyc) + 273.15 !RAD3
        iTemp%TCND_P(2,3,icyc) = t_prt(81,icyc) + 273.15 !RAD3


!     Detector
        iTemp%Tdet(1:4,1,icyc) = t_prt(42:45,icyc) + 273.15 ! RAD1 V/P/M/H
        iTemp%Tdet(1:4,2,icyc) = t_prt(46:49,icyc) + 273.15 ! RAD2 V/P/M/H
        iTemp%Tdet(1:4,3,icyc) = t_prt(50:53,icyc) + 273.15 ! RAD3 V/P/M/H

        enddo  !icyc

        return
        end subroutine get_zone_temperatures



!     Step 2: non-linearity correction
!     non linearity correction
!     the raw counts are divided (normalized) to a single step accumulation 
!     SA1 and SA2 are divided by 2, LA1-4 are dvided by 10, LA5-8 are divided by 2.

        subroutine nl_correction(n_cyc, Tdet, plc, S_acc_raw, L_acc_raw, S_acc, L_acc)
        use l2gen_aquarius_module
        implicit none

        type(prelaunch_struct) :: plc

!     Short accummulation radiometer raw counts 
!     6 records for each of the 12 subcycles. 
!     S1 and S2 are double accumulatons and will be divided by 2 during pocessing.
!     Polarization order (last dimension) is  1=V, 2=P, 3=M, 4=H
        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc_raw      ! raw SA counts
        real(4),    dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc  ! non-linearity corrected SA counts 

!     Long accumulation radiometer raw counts
!     8 records per cycle. 
!     L1-L4 are to be divided by 10. L5-L8 are to be divided by 2. The divisions will be performed during processing. 
        real(4), dimension(n_Lacc,           npol, n_rad, max_cyc) :: L_acc_raw      ! raw LA counts
        real(4),    dimension(n_Lacc,           npol, n_rad, max_cyc) :: L_acc ! non-linearity corrected LA counts

        real(4), dimension(4,n_rad,max_cyc) :: Tdet

        integer(4) icyc,irad,ipol, is, iacc
        integer(4) n_cyc

        real(4) vave, c2, c3, delta_t
        integer(4), parameter :: inl=1 ! 0= no non-lin  1= with non-lin 


        S_acc=missing_val
        L_acc=missing_val

        do icyc=1,n_cyc
           do irad=1,n_rad
              do ipol=1,4

                 if (inl.ne.1) then
                    c2=0.0
                    c3=0.0
                 else
                    if (abs(Tdet(ipol,irad,icyc)-missing_val)<0.1) cycle ! invalid detector temperature
                    delta_t = Tdet(ipol,irad,icyc) - plc%Tref_nl(ipol,irad)
!                   write(*,*) Tdet(ipol,irad,icyc), plc%Tref_nl(ipol,irad)
                    c2  = plc%Dnl20(ipol,irad) + plc%Dnl21(ipol,irad)*delta_t + plc%Dnl22(ipol,irad)*(delta_t**2)
                    c3  = plc%Dnl30(ipol,irad) + plc%Dnl31(ipol,irad)*delta_t + plc%Dnl32(ipol,irad)*(delta_t**2)
                 endif

                 do is=1,n_subcyc
                    do iacc=1,n_Sacc
                       if (S_acc_raw(iacc,is,ipol,irad,icyc).eq.missing_count) cycle
                       vave=S_acc_raw(iacc,is,ipol,irad,icyc)
                       S_acc(iacc,is,ipol,irad,icyc) = vave + c2*(vave**2) + c3*(vave**3)       
                    enddo       !iacc
                 enddo          !is
        
                 do iacc=1,n_Lacc
                    if (L_acc_raw(iacc,ipol,irad,icyc).eq.missing_count) cycle
                    vave=L_acc_raw(iacc,ipol,irad,icyc) !/2.0
                    L_acc(iacc,ipol,irad,icyc) = vave + c2*(vave**2) + c3*(vave**3)
                 enddo          !iacc
        
              enddo             !ipol
           enddo                !irad
        enddo                   !icyc

                
        return
        end     subroutine nl_correction 



!     Step 3: Gain Glitch Detector
        subroutine gain_glitch_detector(n_cyc, rfi, L_acc, m_gain, m_off, ntot)
        use l2gen_aquarius_module
        implicit none

        integer(4) n_cyc
        integer(4) :: NTOT
        integer(1), dimension(max_cyc) :: jflag
        real(4), dimension(n_cyc) :: Y1, Y2, Z
        real(4) :: xsum

        integer(4) :: icyc, irad, kpol, ibox, i1, i2, iwin
        integer(4) :: nwin1, nwin2

        real(4), dimension(n_Lacc,           npol, n_rad, max_cyc) :: L_acc
        integer(4), dimension(               npol, n_rad, max_cyc) :: m_gain, m_off ! dynamical values for gain/offset 

        type(rfi_struct) :: rfi

        rfi%iflag_glitch=0

        if (mod(N1,2)==0) then
                nwin1=(N1/2)-1
        else
                nwin1=(N1-1)/2
        endif 

        if (mod(N2,2)==0) then
                nwin2=(N2/2)-1
        else
                nwin2=(N2-1)/2
        endif
        
         

        do irad=1,n_rad
           do kpol=1,npol       !V/P/M/H

              jflag=0           ! reset flag

!             boxcar filter
              do icyc=1,n_cyc
                 ntot=0
                 xsum=0.0
                 do ibox=icyc-nwin1, icyc+nwin1
                    if (ibox<1 .or. ibox>n_cyc) cycle
                    if (abs(L_acc(1,kpol,irad,ibox)-missing_val)<0.1) cycle
                    ntot=ntot+1
                    xsum = xsum + L_acc(1,kpol,irad,ibox)
                 enddo                           
                 if (ntot>0) then
                    Y1(icyc)=xsum/ntot
                 else
                    Y1(icyc)=missing_val
                 endif
              enddo             !icyc   


!             differential filter
              do icyc=1,n_cyc
                 i2 = icyc+nwin2
                 i1 = icyc-nwin2         
                 if  (i2>n_cyc .or. i2<1 .or. i1>n_cyc .or. i1<1) then
                    Y2(icyc)=missing_val
                    cycle
                 else if (abs(Y1(i2)-missing_val)<0.1 .or. abs(Y1(i1)-missing_val)<0.1 ) then
                    Y2(icyc) = missing_val
                 else
                    Y2(icyc) = Y1(i2) - Y1(i1)
                 endif
              enddo             ! icyc

!             threshold check
              do icyc=1,n_cyc
                 if (abs(Y2(icyc)-missing_val)<0.1) then
                    Z(icyc)=missing_val
                 else
                    Z(icyc) = abs(Y2(icyc))/sigma_glitch(kpol,irad)
                 endif
                 if (abs(Z(icyc)-missing_val)<0.1) cycle
                 if (Z(icyc)>threshold_glitch(kpol,irad)) jflag(icyc)=1
              enddo             ! icyc

!            flag N2/2 samples to the left and to the right
              do icyc=1,n_cyc
                 if (jflag(icyc)==1) then
                    i2 = icyc+nwin2
                    i1 = icyc-nwin2     
                    do iwin=i1,i2
                       if (iwin>=1 .and. iwin<=n_cyc) rfi%iflag_glitch(kpol,irad,iwin)=1
                    enddo 
                 endif          ! sample icyc flagged
              enddo             ! icyc            
                
!            set dynamical gain/offset averaging windows
              do icyc=1,n_cyc
                 if (rfi%iflag_glitch(kpol,irad,icyc)==0) then
                    m_gain(kpol,irad,icyc)   = n_gain
                    m_off(kpol,irad,icyc)    = n_off
                 else
                    m_gain(kpol,irad,icyc)   = n_red
                    m_off(kpol,irad,icyc)    = n_red
                 endif
              enddo             ! icyc
                        
           enddo                ! kpol
        enddo                   ! irad                                    
                                
        return
        end subroutine gain_glitch_detector


!     Step 5a: calculates averaged gain and offset for V/H
        subroutine find_cal_VH(n_cyc, m_gain, m_off, iTemp, gvv, ghh, ov, oh, L_acc)
        use l2gen_aquarius_module
        implicit none

        integer(4) icyc, jcyc, irad, ipol, kpol,nsum
        integer(4) n_cyc

        real(4) t_0, t_1, c_0, c_1, d_t 
        real(8) xsum

        integer(4) ivalid(2,max_cyc)
        real(8)     xgain(2,max_cyc), xoff(2,max_cyc)

        real(4), dimension(n_rad,max_cyc) :: gvv, ghh, ov, oh

!       real(4), dimension(2,n_rad,max_cyc) :: T0, TND

        real(4), dimension(n_Lacc, npol, n_rad, max_cyc) :: L_acc ! non-linearity corrected LA counts

        integer(4), dimension(     npol, n_rad, max_cyc) :: m_gain, m_off ! dynamical values for gain/offset 

        type(instrument_temp_struct) :: iTemp

!       write (*,*) 'gvv ', loc(gvv)
!       write (*,*) 'iTemp ', loc(iTemp%t1(1,1)), loc(iTemp)

        do irad=1,n_rad ! radiometer loop

           do ipol=1,2
              if (ipol==1) kpol=1
              if (ipol==2) kpol=4
                        
! first pass: calculate gain and offset for each cycle
              do icyc=1,n_cyc

                 ivalid(ipol,icyc)=0

                 t_0 = iTemp%T0(ipol, irad,  icyc)
                 d_t = iTemp%TND(ipol, irad,  icyc)

                 if(abs(t_0-missing_val).lt.0.1) cycle 
                 if(abs(d_t-missing_val).lt.0.1) cycle
                                
                 t_1 = t_0 + d_t
                                 
                 if (ipol==1) then  
                    call avg2(L_acc(1, kpol,irad,icyc), L_acc(4, kpol,irad,icyc), missing_val, c_0) ! v pol w/o ND (cold)
                    if(abs(c_0-missing_val).lt.0.1) cycle
                              
                    call avg2(L_acc(2, kpol,irad,icyc), L_acc(3, kpol,irad,icyc), missing_val, c_1) ! v pol with ND (hot)
                    if(abs(c_1-missing_val).lt.0.1) cycle
                 endif          ! vpol

                 if (ipol==2) then
                    call avg2(L_acc(1, kpol,irad,icyc), L_acc(2, kpol,irad,icyc), missing_val, c_0) ! h pol w/o ND (cold)
                    if(abs(c_0-missing_val).lt.0.1) cycle

                    call avg2(L_acc(3, kpol,irad,icyc), L_acc(4, kpol,irad,icyc), missing_val, c_1) ! h pol with ND (hot)
                    if(abs(c_1-missing_val).lt.0.1) cycle
                 endif          ! hpol

                 ivalid(ipol,icyc) =  1                          
                 xgain(ipol,icyc) = (c_1-c_0)/d_t
                 xoff(ipol,icyc) =  c_0-xgain(ipol,icyc)*t_0

              enddo             ! icyc, first pass through cycles
           enddo                ! ipol loop
                        
! second pass through cycles: averaging
           do icyc=1,n_cyc 
              xsum=0.d0
              nsum=0
              do jcyc=icyc-m_gain(1,irad,icyc)/2, icyc+m_gain(1,irad,icyc)/2
                 if (jcyc<1 .or. jcyc>n_cyc) cycle
                 if (ivalid(1,jcyc)==0)      cycle
                 nsum=nsum+1
                 xsum=xsum+xgain(1,jcyc)
              enddo

              if (nsum.ne.0) then
                 gvv(irad,icyc)=xsum/nsum
              else
                 gvv(irad,icyc)=missing_val
              endif                           

              xsum=0.d0
              nsum=0
              do jcyc=icyc-m_gain(4,irad,icyc)/2,icyc+m_gain(4,irad,icyc)/2
                 if (jcyc<1 .or. jcyc>n_cyc) cycle
                 if (ivalid(2,jcyc)==0) cycle
                 nsum=nsum+1
                 xsum=xsum+xgain(2,jcyc)
              enddo
              if (nsum.ne.0) then
                 ghh(irad,icyc)=xsum/nsum
              else
                 ghh(irad,icyc)=missing_val
              endif

              xsum=0.d0
              nsum=0
              do jcyc=icyc-m_off(1,irad,icyc)/2, icyc+m_off(1,irad,icyc)/2
                 if (jcyc<1 .or. jcyc>n_cyc) cycle
                 if (ivalid(1,jcyc)==0) cycle
                 nsum=nsum+1
                 xsum=xsum+xoff(1,jcyc)
              enddo
              if (nsum.ne.0) then
                 ov(irad,icyc)=xsum/nsum
              else
                 ov(irad,icyc)=missing_val
              endif

              xsum=0.d0
              nsum=0
              do jcyc=icyc-m_off(1,irad,icyc)/2, icyc+m_off(1,irad,icyc)/2
                 if (jcyc<1 .or. jcyc>n_cyc) cycle
                 if (ivalid(2,jcyc)==0) cycle
                 nsum=nsum+1
                 xsum=xsum+xoff(2,jcyc)
              enddo
              if (nsum.ne.0) then
                 oh(irad,icyc)=xsum/nsum
              else
                 oh(irad,icyc)=missing_val
              endif

           enddo                ! second pass through cycles
           
        enddo                   ! radiometer loop

        return
        end subroutine find_cal_VH



!     Step 5b: calculates averaged gain and offset coefficients for P/M channels
        subroutine find_cal_PM(n_cyc, m_gain, m_off, T0, TND, TDL_offset, TND_offset, gpv, gph, op, om, gmv, gmh, gpp, gmm, L_acc)
        use l2gen_aquarius_module
        implicit none

        integer(4) icyc, jcyc, irad, ipol, nsum
        integer(4) ivalid(2,max_cyc)
        integer(4) n_cyc

        real(4) t0v, t0v_offset, tndv, tndv_offset, t0h, t0h_offset, tndh, tndh_offset   
        real(4) gvec(2, 3,max_cyc)                                                              

        real(8) xsum(2), osum, det
        real(8) amat(4,3), amat_t(3,4), xmat(3,3), xmat_inv(3,3), xmat_tes(3,3), bvec(4),cvec(3)

        real(4), dimension(2,n_rad,max_cyc) :: TND
        real(4), dimension(2,n_rad,max_cyc) :: T0

        real(4), dimension(2,n_rad,max_cyc) :: TDL_offset
        real(4), dimension(2,n_rad,max_cyc) :: TND_offset

        real(4), dimension(n_rad,max_cyc) :: gpv, gph, op, om, gmv, gmh, gpp, gmm

        real(4), dimension(n_Lacc,           npol, n_rad, max_cyc) :: L_acc ! non-linearity corrected LA counts
        integer(4), dimension(               npol, n_rad, max_cyc) :: m_gain, m_off ! dynamical values for gain/offset 

        do irad=1,n_rad ! radiometer loop
                        
                ! first pass: calculate gain and offset for each cycle
                do icyc=1,n_cyc

                        ivalid(1:2,icyc)=0

                        t0v =  T0(1, irad, icyc)
                        t0h =  T0(2, irad, icyc)
                        t0v_offset =  TDL_offset(1, irad, icyc)
                        t0h_offset =  TDL_offset(2, irad, icyc)
                        tndv = TND(1, irad, icyc)
                        tndh = TND(2, irad, icyc)
                        tndv_offset = TND_offset(1, irad, icyc)
                        tndh_offset = TND_offset(2, irad, icyc)

                        if(abs(t0v-missing_val)<0.1)    cycle 
                        if(abs(t0h-missing_val)<0.1)    cycle 
                        if(abs(t0v_offset-missing_val)<0.1)     cycle
                        if(abs(t0h_offset-missing_val)<0.1)     cycle
                        if(abs(tndv-missing_val)<0.1)   cycle
                        if(abs(tndh-missing_val)<0.1)   cycle                   
                        if(abs(tndv_offset-missing_val)<0.1)    cycle
                        if(abs(tndh_offset-missing_val)<0.1)    cycle   
                                                        
!             LA1
                        amat(1,1) = t0v - t0v_offset
                        amat(1,2) = t0h - t0h_offset
                        amat(1,3) = 1.d0

!             LA2            
                        amat(2,1) = t0v               + tndv - tndv_offset 
                        amat(2,2) = t0h + t0h_offset  
                        amat(2,3) = 1.d0

!             LA3
                        amat(3,1) = t0v               + tndv + tndv_offset
                        amat(3,2) = t0h               + tndh - tndh_offset     
                        amat(3,3) = 1.d0

!             LA4
                        amat(4,1) = t0v + t0v_offset   
                        amat(4,2) = t0h               + tndh + tndh_offset           
                        amat(4,3) = 1.d0

                        amat_t=transpose(amat)
                        xmat = matmul(amat_t,amat)

                        call invert_3by3(xmat,  xmat_inv,det)
                        xmat_tes=matmul(xmat,xmat_inv)
                        xmat_tes(1,1)=xmat_tes(1,1)-1
                        xmat_tes(2,2)=xmat_tes(2,2)-1
                        xmat_tes(3,3)=xmat_tes(3,3)-1
                        if(abs(det).lt.1.d-10 .or. maxval(abs(xmat_tes)).gt.1.d-6) write(*,*) 'error in invert_3by3, pgm stopped'

                        ! P
                ipol=2
                        if(abs(L_acc(1, ipol,irad,icyc)-missing_val)>0.1 .and. abs(L_acc(2, ipol,irad,icyc)-missing_val)>0.1 .and. 
     &             abs(L_acc(3, ipol,irad,icyc)-missing_val)>0.1 .and. abs(L_acc(4, ipol,irad,icyc)-missing_val)>0.1 ) then
                        bvec(1:4) = L_acc(1:4,ipol,irad,icyc)
                        cvec(1:3) = matmul(amat_t,bvec)
                        gvec(1,1:3,icyc) = matmul(xmat_inv,cvec)
                        ivalid(1,icyc)=1
                        endif                   

                        ! M
                ipol=3
                        if(abs(L_acc(1, ipol,irad,icyc)-missing_val)>0.1 .and. abs(L_acc(2, ipol,irad,icyc)-missing_val)>0.1 .and. 
     &             abs(L_acc(3, ipol,irad,icyc)-missing_val)>0.1 .and. abs(L_acc(4, ipol,irad,icyc)-missing_val)>0.1 ) then
                        bvec(1:4) = L_acc(1:4,ipol,irad,icyc)
                        cvec(1:3) = matmul(amat_t,bvec)
                        gvec(2,1:3,icyc) = matmul(xmat_inv,cvec)
                        ivalid(2,icyc)=1
                        endif

                enddo ! first pass through cycles
                
                ! second pass through cycles: averaging
                do icyc=1,n_cyc 

                        ! gpv, gph
                        xsum=0.d0
                        nsum=0
                        do jcyc=icyc-m_gain(2,irad,icyc)/2, icyc+m_gain(2,irad,icyc)/2
!                       do jcyc=icyc-n_gain/2,icyc+n_gain/2
                                if (jcyc<1 .or. jcyc>n_cyc) cycle
                                if (ivalid(1,jcyc)==0) cycle
                                nsum=nsum+1
                                xsum=xsum+gvec(1,1:2,jcyc)
                        enddo
                        if (nsum.ne.0) then
                                gpv(irad,icyc)=xsum(1)/nsum
                                gph(irad,icyc)=xsum(2)/nsum
                                gpp(irad,icyc)=gpv(irad,icyc)+gph(irad,icyc) ! only used for calculating sigma_count in RFI detector
                        else
                                gpv(irad,icyc)=missing_val
                                gph(irad,icyc)=missing_val
                                gpp(irad,icyc)=missing_val
                        endif                           

                        ! op
                        osum=0.d0
                        nsum=0
 
                        do jcyc=icyc-m_off(2,irad,icyc)/2, icyc+m_off(2,irad,icyc)/2
!                       do jcyc=icyc-n_off/2,icyc+n_off/2
                                if (jcyc<1 .or. jcyc>n_cyc) cycle
                                if (ivalid(1,jcyc)==0) cycle
                                nsum=nsum+1
                                osum=osum+gvec(1,3,jcyc)
                        enddo
                        if (nsum.ne.0) then
                                op(irad,icyc)=osum/nsum
                        else
                                op(irad,icyc)=missing_val
                        endif

                        ! gmv, gmh
                        xsum=0.d0
                        nsum=0
                        do jcyc=icyc-m_gain(3,irad,icyc)/2, icyc+m_gain(3,irad,icyc)/2
!                       do jcyc=icyc-n_gain/2,icyc+n_gain/2
                                if (jcyc<1 .or. jcyc>n_cyc) cycle
                                if (ivalid(2,jcyc)==0) cycle
                                nsum=nsum+1
                                xsum=xsum+gvec(2,1:2,jcyc)
                        enddo
                        if (nsum.ne.0) then
                                gmv(irad,icyc)=xsum(1)/nsum
                                gmh(irad,icyc)=xsum(2)/nsum
                                gmm(irad,icyc)=gmv(irad,icyc)+gmh(irad,icyc) ! only used for calculating sigma_count in RFI detector
                        else
                                gmv(irad,icyc)=missing_val
                                gmh(irad,icyc)=missing_val
                                gmm(irad,icyc)=missing_val
                        endif                           

                        ! om
                        osum=0.d0
                        nsum=0
                        do jcyc=icyc-m_off(3,irad,icyc)/2, icyc+m_off(3,irad,icyc)/2
!                       do jcyc=icyc-n_off/2,icyc+n_off/2
                                if (jcyc<1 .or. jcyc>n_cyc) cycle
                                if (ivalid(2,jcyc)==0) cycle
                                nsum=nsum+1
                                osum=osum+gvec(2,3,jcyc)
                        enddo
                        if (nsum.ne.0) then
                                om(irad,icyc)=osum/nsum
                        else
                                om(irad,icyc)=missing_val
                        endif

                enddo ! second pass through cycles

        enddo ! radiometer loop

        return
        end subroutine find_cal_PM




!     Step 7a: calculates averaged gains hpU, hmU without RFI filter
        subroutine find_cal_hU(n_cyc, CL5, CLMM, dL5_dT5, T5, TCND, Tref_5, gainOff, S_acc) 
        use l2gen_aquarius_module
        implicit none

        real(4), dimension(2,n_rad)     ::     CL5, CLMM, dL5_dT5
        real(4), dimension(2,n_rad,max_cyc) :: T5
        real(4), dimension(2,n_rad,max_cyc) :: TCND

        real(4)         ::      Tref_5 ! ref. temp for FE losses and ND

        integer(4) n_cyc
        integer(4)      icyc, jcyc, irad, ipol, kpol, is, iacc,nsum 
        integer(4)      ivalid(2,max_cyc)

        real(4) XTCND_U, XTCND_U_prime
        real(4) L5(2),XTCND_prime(2), LMM(2)
        real(4) xvcnd(2), xv(2), ygv(2), ygh(2)
        real(4) xgU(2,max_cyc)

        real(8) xsum
        real(4) cosd

        type(gain_off_struct) :: gainOff
        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc     ! non-linearity corrected SA counts 

        do irad=1,n_rad ! radiometer loop
                        
! first pass: calculate gain for each cycle
           do icyc=1,n_cyc

              ivalid(1:2,icyc)=0

              if(abs( TCND(1, irad, icyc)-missing_val)<0.1)   cycle 
              if(abs( TCND(2, irad, icyc)-missing_val)<0.1)   cycle 
              if(abs(   T5(1, irad, icyc)-missing_val)<0.1)   cycle 
              if(abs(   T5(2, irad, icyc)-missing_val)<0.1)   cycle 


              XTCND_U = 2.0*cosd(delta_phi(irad))*sqrt(TCND(1,irad,icyc)*TCND(2,irad,icyc))
              L5(1:2)=  CL5(1:2,irad) + dL5_dT5(1:2,irad)*(T5(1:2,irad,icyc) - Tref_5)
              LMM(1:2)= CLMM(1:2,irad)
              if (L5(1)< 1.0E-15 .or. L5(2) < 1.0E-15) 
     &             write(*,*) ' loss factor 0 or negative. pgm stopped in find_cal_U' 
              if (LMM(1)< 1.0E-15 .or. LMM(2) < 1.0E-15) 
     &             write(*,*) ' impedence mismatch 0 or negative. pgm stopped in find_cal_U' 
              XTCND_U_prime = XTCND_U/sqrt(L5(1)*L5(2) * LMM(1)*LMM(2))
              if (abs(XTCND_U_prime) < 1.0E-15) write(*,*) ' TCND prime is zero. pgm stopped in find_cal_U' 
              XTCND_prime(1:2) = TCND(1:2,irad,icyc)/(L5(1:2)*LMM(1:2))
                        
              do ipol=1,2       ! index for P/M
                        
                 if(ipol==1) kpol=2 !P
                 if(ipol==2) kpol=3 !M

                 if(ipol==1 .and. abs(gainOff%gpv(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                 if(ipol==1 .and. abs(gainOff%gph(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                 if(ipol==2 .and. abs(gainOff%gmv(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                 if(ipol==2 .and. abs(gainOff%gmh(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                        
                 if(ipol==1) then
                    ygv(ipol)=gainOff%gpv(irad,icyc)
                    ygh(ipol)=gainOff%gph(irad,icyc)
                 else
                    ygv(ipol)=gainOff%gmv(irad,icyc)
                    ygh(ipol)=gainOff%gmh(irad,icyc)
                 endif                   
                        
! average voltage of antena + CND view over subcycles
                 nsum=0
                 xsum=0.d0
                 do is=1,n_subcyc
                    if (abs(S_acc(6, is,kpol,irad,icyc)-missing_val)<0.1) cycle ! subcycle loop
                    nsum=nsum + 1
                    xsum=xsum + S_acc(6, is,kpol,irad,icyc)
                 enddo
                 if (nsum.eq.0) cycle ! pol loop 
                 xvcnd(ipol) = xsum/nsum
                        
! average voltage of antenna view over subcycles and accumulations
                 nsum=0
                 xsum=0.d0
                 do is=1,n_subcyc
                    do iacc=1,5
                       if (abs(S_acc(iacc,is,kpol,irad,icyc)-missing_val)<0.1) cycle ! subcycle loop
                       nsum=nsum + 1
                       xsum=xsum + S_acc(iacc,is,kpol,irad,icyc)
                    enddo
                 enddo
                 if (nsum.eq.0) cycle ! pol loop 
                 xv(ipol) = xsum/nsum

                 ivalid(ipol,icyc) = 1
                 xgU(ipol,icyc)=(xvcnd(ipol)-xv(ipol))/XTCND_U_prime-(XTCND_prime(1)*
     &                ygv(ipol)+XTCND_prime(2)*ygh(ipol))/XTCND_U_prime
                                                                
              enddo             ! ipol loop P/M
                        
           enddo                ! icyc, first pass through cycles 



! second pass through cycles: averaging
           do ipol=1,2          ! P/M
              do icyc=1,n_cyc 

                 xsum=0.d0
                 nsum=0

                 do jcyc=icyc-n_gain/2,icyc+n_gain/2
                    if (jcyc<1 .or. jcyc>n_cyc) cycle
                    if (ivalid(ipol,jcyc)==0) cycle
                    nsum=nsum+1
                    xsum=xsum+xgU(ipol,jcyc)
                 enddo

                 if (ipol==1) then
                    if (nsum.ne.0) then
                       gainOff%hpU(irad,icyc)=xsum/nsum
                    else
                       gainOff%hpU(irad,icyc)=missing_val
                    endif           
                 endif

                 if (ipol==2) then
                    if (nsum.ne.0) then
                       gainOff%hmU(irad,icyc)=xsum/nsum
                    else
                       gainOff%hmU(irad,icyc)=missing_val
                    endif           
                 endif

              enddo             ! cycles
           enddo                ! pol P/M

        enddo                   ! radiometer loop

        return
        end subroutine find_cal_hU



!     Step 7b: calculates averaged gains gpU, gmU with RFI filter
!     need to run RFI detector first
        subroutine find_cal_gU(n_cyc, CL5, CLMM, dL5_dT5, T5, TCND, Tref_5, gainOff, rfi, S_acc)
        use l2gen_aquarius_module
        implicit none

        real(4), dimension(2,n_rad)     ::     CL5, CLMM, dL5_dT5
        real(4), dimension(2,n_rad,max_cyc) :: T5
        real(4), dimension(2,n_rad,max_cyc) :: TCND

        real(4)         ::      Tref_5 ! ref. temp for FE losses and ND

        integer(4) n_cyc
        integer(4)      icyc, jcyc, irad, ipol, kpol, is, iacc,nsum 
        integer(4)      ivalid(2,max_cyc)

        real(4) XTCND_U, XTCND_U_prime
        real(4) L5(2),XTCND_prime(2), LMM(2)
        real(4) xvcnd(2), xv(2), ygv(2), ygh(2)
        real(4) xgU(2,max_cyc)

        real(8) xsum
        real(4) cosd

        type(gain_off_struct) :: gainOff
        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc     ! non-linearity corrected SA counts 

        type(rfi_struct) :: rfi

        do irad=1,n_rad         ! radiometer loop
                        
! first pass: calculate gain for each cycle
           do icyc=1,n_cyc

              ivalid(1:2,icyc)=0

              if(abs( TCND(1, irad, icyc)-missing_val)<0.1)   cycle 
              if(abs( TCND(2, irad, icyc)-missing_val)<0.1)   cycle 
              if(abs(   T5(1, irad, icyc)-missing_val)<0.1)   cycle 
              if(abs(   T5(2, irad, icyc)-missing_val)<0.1)   cycle 


              XTCND_U = 2.0*cosd(delta_phi(irad))*sqrt(TCND(1,irad,icyc)*TCND(2,irad,icyc))
              L5(1:2)=  CL5(1:2,irad) + dL5_dT5(1:2,irad)*(T5(1:2,irad,icyc) - Tref_5)
              LMM(1:2)= CLMM(1:2,irad)
              if (L5(1)< 1.0E-15 .or. L5(2) < 1.0E-15) write(*,*) ' loss factor 0 or negative. pgm stopped in find_cal_U' 
              if (LMM(1)< 1.0E-15 .or. LMM(2) < 1.0E-15) 
     &             write(*,*) ' impedence mismatch 0 or negative. pgm stopped in find_cal_U' 
              XTCND_U_prime = XTCND_U/sqrt(L5(1)*L5(2) * LMM(1)*LMM(2))
              if (abs(XTCND_U_prime) < 1.0E-15) then
                 write(*,*) 'icyc: ', icyc, 'irad: ', irad
                 write(*,*) ' TCND prime is zero. pgm stopped in find_cal_gU' 
              endif
              XTCND_prime(1:2) = TCND(1:2,irad,icyc)/(L5(1:2)*LMM(1:2))
                        
              do ipol=1,2       ! index for P/M
                        
                 if(ipol==1) kpol=2 !P
                 if(ipol==2) kpol=3 !M

                 if(ipol==1 .and. abs(gainOff%gpv(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                 if(ipol==1 .and. abs(gainOff%gph(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                 if(ipol==2 .and. abs(gainOff%gmv(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                 if(ipol==2 .and. abs(gainOff%gmh(irad,icyc)-missing_val)<0.1) cycle ! pol loop  
                        
                 if(ipol==1) then
                    ygv(ipol)=gainOff%gpv(irad,icyc)
                    ygh(ipol)=gainOff%gph(irad,icyc)
                 else
                    ygv(ipol)=gainOff%gmv(irad,icyc)
                    ygh(ipol)=gainOff%gmh(irad,icyc)
                 endif                   
                        
! average voltage of antenna + CND view over subcycles
                 nsum=0
                 xsum=0.d0
                 do is=1,n_subcyc
                    if (abs(S_acc(6, is,kpol,irad,icyc)-missing_val)<0.1) cycle ! subcycle loop
                    if (rfi%iflag_rfi(6, is,kpol,irad,icyc)==1)  cycle ! RFI in CND observation
                    if (rfi%iflag_rfi_CND(is,kpol,irad,icyc)==1) cycle ! RFI in CND observation

                    nsum=nsum + 1
                    xsum=xsum + S_acc(6, is,kpol,irad,icyc)
                 enddo
                 if (nsum.eq.0) cycle ! pol loop 
                 xvcnd(ipol) = xsum/nsum
                        
! average voltage of antenna view over subcycles and accumulations
                 nsum=0
                 xsum=0.d0
                 do is=1,n_subcyc
                    do iacc=1,5
                       if (abs(S_acc(iacc,is,kpol,irad,icyc)-missing_val)<0.1) cycle ! subcycle loop
                       if (rfi%iflag_rfi(iacc,is,kpol,irad,icyc)==1) cycle ! RFI 
                       nsum=nsum + 1
                       xsum=xsum + S_acc(iacc,is,kpol,irad,icyc)
                    enddo
                 enddo
                 if (nsum.eq.0) cycle ! pol loop 
                 xv(ipol) = xsum/nsum

                 ivalid(ipol,icyc) = 1
                 xgU(ipol,icyc)=(xvcnd(ipol)-xv(ipol))/XTCND_U_prime-(XTCND_prime(1)*
     &                ygv(ipol)+XTCND_prime(2)*ygh(ipol))/XTCND_U_prime
                                                                
              enddo             ! ipol loop P/M
                        
           enddo                ! icyc, first pass through cycles 



! second pass through cycles: averaging
           do ipol=1,2          ! P/M
              do icyc=1,n_cyc 

                 xsum=0.d0
                 nsum=0

                 do jcyc=icyc-n_gain/2,icyc+n_gain/2
                    if (jcyc<1 .or. jcyc>n_cyc) cycle
                    if (ivalid(ipol,jcyc)==0) cycle
                    nsum=nsum+1
                    xsum=xsum+xgU(ipol,jcyc)
                 enddo

                 if (ipol==1) then
                    if (nsum.ne.0) then
                       gainOff%gpU(irad,icyc)=xsum/nsum
                    else
                       gainOff%gpU(irad,icyc)=missing_val
                    endif           
                 endif

                 if (ipol==2) then
                    if (nsum.ne.0) then
                       gainOff%gmU(irad,icyc)=xsum/nsum
                    else
                       gainOff%gmU(irad,icyc)=missing_val
                    endif           
                 endif

              enddo             ! cycles
           enddo                ! pol P/M

        enddo                   ! radiometer loop

        return
        end subroutine find_cal_gU



!     Step 8: count to TA before front end loss corrections
        subroutine find_TA_hat(n_cyc, gainOff, rfi, TA_hat, TF_hat, S_acc)
        use l2gen_aquarius_module
        implicit none

        integer(4), parameter :: ntot=n_subcyc*5 !total number of earth counts per cycle

        integer(4) icyc,irad,ipol, is, iacc, nsum, ibad, icase
        integer(4) n_cyc
        real(4) vave(4) ! counts averaged over one full cycle
        real(4) V_vec(2), O_vec(2), GU_Vec(2), Z_vec(2)
        real(4) tax(2),G_mat(2,2),xpU,xmU       
        real(8) xsum    

        type(gain_off_struct) :: gainOff
        type(rfi_struct) :: rfi

        real(4), dimension(3,n_rad,max_cyc)  :: TF_hat, TA_hat  

        real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc     ! non-linearity corrected SA counts 

        TA_hat=missing_val      
        TF_hat=missing_val      ! Added 12/2011 JMG to fix negative Tf bug

        do icase=1,2            !1 is no rfi filtering, 2 is yes rfi filtering

           do icyc=1,n_cyc
              do irad=1,n_rad

                 if(icase.eq.1) then     
                    xpU=gainOff%hpU(irad,icyc)
                    xmU=gainOff%hmU(irad,icyc)
                 else
                    xpU=gainOff%gpU(irad,icyc)
                    xmU=gainoff%gmU(irad,icyc)
                 endif


                 if (abs(gainOff%gvv(irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%ghh(irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%gpv(irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%gph(irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%gmv(irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%gmh(irad,icyc)-missing_val)<0.1) cycle
                 if (abs(xpU           -missing_val)<0.1)         cycle
                 if (abs(xmU           -missing_val)<0.1)         cycle
                 if (abs(gainOff%ov( irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%oh( irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%op( irad,icyc)-missing_val)<0.1) cycle
                 if (abs(gainOff%om( irad,icyc)-missing_val)<0.1) cycle

                 if (icase.eq.2) then
                    if (abs(gainOff%gpU(irad,icyc)-missing_val)<0.1) cycle
                    if (abs(gainOff%gmU(irad,icyc)-missing_val)<0.1) cycle
                 endif

                 if (abs(xmU) < 1.0E-15 .or. abs(xpU) < 1.0E-15) write(*,*) ' xpU or xmU zero. pgm find_TA_hat stopped.'

                 ibad=0

                 do ipol=1,4
!     sum counts over cycle
                    nsum=0
                    xsum=0.0
                    rfi%num_rfi(ipol,irad,icyc)=0
                    do is=1,n_subcyc
                       do iacc=1,5
                          if(abs(S_acc(iacc,is,ipol,irad,icyc)-missing_val)<0.1)  cycle
                   
                          if(icase.eq.2 .and. rfi%iflag_rfi(iacc,is,ipol,irad,icyc)==1)  then
                             rfi%num_rfi(ipol,irad,icyc)=rfi%num_rfi(ipol,irad,icyc)+1
                             cycle ! RFI
                          endif 

                          nsum=nsum + 1
                          xsum=xsum + S_acc(iacc,is,ipol,irad,icyc)
                       enddo
                    enddo

                    if (nsum.ne.0) then
                       vave(ipol)=xsum/nsum
                    else
                       vave(ipol)=missing_val
                       ibad=1
                    endif

                 enddo          ! pol
        
                 if(ibad.eq.1) cycle !if any of the four pol count values are missing, dont compute ta  

                 tax(1) = (vave(1)-gainOff%ov(irad,icyc))/gainOff%gvv(irad,icyc)
                 tax(2) = (vave(4)-gainOff%oh(irad,icyc))/gainOff%ghh(irad,icyc)

                 V_vec = (/vave(2), vave(3)/)
                 O_vec = (/gainOff%op(irad,icyc), gainOff%om(irad,icyc)/)
                 G_mat(1,1:2) = (/gainOff%gpv(irad,icyc), gainOff%gph(irad,icyc)/)
                 G_mat(2,1:2) = (/gainOff%gmv(irad,icyc), gainOff%gmh(irad,icyc)/)

!                 GU_vec(1) = gainOff%gpU(irad,icyc)/(gainOff%gmU(irad,icyc)**2 + gainOff%gpU(irad,icyc)**2)
!                 GU_vec(2) = gainOff%gmU(irad,icyc)/(gainOff%gmU(irad,icyc)**2 + gainOff%gpU(irad,icyc)**2)

                 GU_vec(1) = xpU/(xmU**2 + xpU**2)
                 GU_vec(2) = xmU/(xmU**2 + xpU**2)

                 Z_vec = V_vec - O_vec - matmul(G_mat,tax)

                 if(icase.eq.1) then
                    TA_hat(1:2, irad, icyc) = tax
                    TA_hat(  3, irad, icyc) = dot_product(GU_vec, Z_vec)
                 else
                    TF_hat(1:2, irad, icyc) = tax
                    TF_hat(  3, irad, icyc) = dot_product(GU_vec, Z_vec)
                 endif

              enddo             !icyc
           enddo                !irad

        enddo                   !icase

        return
        end subroutine find_TA_hat 


 !     Step 9: front end loss corrections: TA_hat into TA
        subroutine fe_loss_corr(n_cyc, plc, iTemp, TA_hat, TF_hat, TF, TA)
        use l2gen_aquarius_module
        implicit none

        integer(4)      icyc,irad,ipol,ibad,icase
        integer(4) n_cyc
        real(4)         L1
        real(4)         L2A(2), L2B(2), L3(2), L4(2), L5(2), LMM(2)
        real(4)         tax(3), TA_p(2), TA_fp(2), TA_fpp(2), TA_opp(2), TA_pp(2), TAA(2)  
        real(4)         TAA_U, TA_U, lfac

        type(prelaunch_struct) :: plc
        type(instrument_temp_struct) :: iTemp

        real(4), dimension(3,n_rad,max_cyc)  :: TA_hat, TF_hat, TA, TF

        real(4) cosd

        TA=missing_val ! default                
        TF=missing_val ! default

        do icyc=1,n_cyc
           do irad=1,n_rad
           
              if(abs(iTemp%T1( irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T2A(irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T2B(irad,icyc)-missing_val)<0.1) cycle
           
              if(abs(iTemp%T3(1,irad,icyc)-missing_val)<0.1 .or. abs(iTemp%T3(2,irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T4(1,irad,icyc)-missing_val)<0.1 .or. abs(iTemp%T4(2,irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T5(1,irad,icyc)-missing_val)<0.1 .or. abs(iTemp%T5(2,irad,icyc)-missing_val)<0.1) cycle


              L1        =       plc%CL1(1,irad)
              L2A(1:2)  =       plc%CL2A(:,irad) 
              L2B(1:2)  =       plc%CL2B(:,irad) 
              L3(1:2)   =       plc%CL3(:,irad)  
              L4(1:2)   =       plc%CL4(:,irad)  + plc%dL4_dT4(:,irad)*(iTemp%T4(1:2,irad,icyc) - plc%Tref_4)
              L5(1:2)   =       plc%CL5(:,irad)  + plc%dL5_dT5(:,irad)*(iTemp%T5(1:2,irad,icyc) - plc%Tref_5)
              LMM(1:2)  =       plc%CLMM(1:2,irad)

              do icase=1,2      !1 is no rfi filtering, 2 is yes rfi filtering

                 if(icase.eq.1) then 
                    tax=TA_hat(:,irad,icyc)
                 else
                    tax=TF_hat(:,irad,icyc)
                 endif

!       find TA with no rfi filtering
                 ibad=0
                 do ipol =1,3
                    if(abs(tax(ipol)-missing_val)<0.1) ibad=1
                 enddo 
              
                 if(ibad.eq.1) cycle
                        
!       V/H
                 TA_p(1:2)   = LMM(1:2)*tax(1:2)- (LMM(1:2)-1.0)*iTemp%TND_P(1:2,irad,icyc) ! impedance mismatch correction
                 TA_fp( 1:2) = L5(1:2)*TA_p(1:2)   - (L5(1:2)-1.0)*iTemp%T5(1:2,irad,icyc) ! Zone 5 correction
                 TA_fpp(1:2) = L4(1:2)*TA_fp(1:2)  - (L4(1:2)-1.0)*iTemp%T4(1:2,irad,icyc) ! Zone 4 correction
                 TA_opp(1:2) = L3(1:2)*TA_fpp(1:2) - (L3(1:2)-1.0)*iTemp%T3(1:2,irad,icyc) ! Zone 3 correction
                 TA_pp( 1:2) = L2A(1:2)*L2B(1:2)*TA_opp(1:2) - L2A(1:2)*(L2B(1:2)-1.0)*iTemp%T2B(irad,icyc) -
     1                        (L2A(1:2)-1.0)*iTemp%T2A(irad,icyc) ! Zone 2A/2B correction
                 TAA(   1:2) = L1*TA_pp(1:2) - (L1-1.0)*iTemp%T1(irad,icyc) ! Zone 1 correction

!       3rd Stokes U
                 TA_U = tax(3)
                 lfac = L1 * sqrt(L2A(1)*L2A(2)) * sqrt(L2B(1)*L2B(2)) * sqrt(L3(1)*L3(2)) * sqrt(L4(1)*L4(2)) * 
     1                       sqrt(L5(1)*L5(2)) * sqrt(LMM(1)*LMM(2))
                 TAA_U = (lfac/cosd(delta_phi_f(irad))) * TA_U
        
                 if(icase.eq.1) then 
                    TA(1:2,irad,icyc) = TAA(1:2)
                    TA(  3,irad,icyc) = TAA_U
                 else
                    TF(1:2,irad,icyc) = TAA(1:2)
                    TF(  3,irad,icyc) = TAA_U
                 endif

              enddo             !icase

           enddo                !irad
        enddo                   !icyc

        return
        end subroutine fe_loss_corr


#if 0
        subroutine fe_loss_fac(n_cyc, plc, iTemp, lossFac)
        use l2gen_aquarius_module
        implicit none

        integer(4)      icyc,irad
        integer(4) n_cyc
        real(4)         L1
        real(4)         L2A(2), L2B(2), L3(2), L4(2), L5(2), LMM(2), lossFac(2,n_rad,n_cyc)

        type(prelaunch_struct) :: plc
        type(instrument_temp_struct) :: iTemp

        do icyc=1,n_cyc
           do irad=1,n_rad
           
              if(abs(iTemp%T1( irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T2A(irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T2B(irad,icyc)-missing_val)<0.1) cycle
           
              if(abs(iTemp%T3(1,irad,icyc)-missing_val)<0.1 .or. abs(iTemp%T3(2,irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T4(1,irad,icyc)-missing_val)<0.1 .or. abs(iTemp%T4(2,irad,icyc)-missing_val)<0.1) cycle
              if(abs(iTemp%T5(1,irad,icyc)-missing_val)<0.1 .or. abs(iTemp%T5(2,irad,icyc)-missing_val)<0.1) cycle


              L1        =       plc%CL1(1,irad)
              L2A(1:2)  =       plc%CL2A(:,irad) 
              L2B(1:2)  =       plc%CL2B(:,irad) 
              L3(1:2)   =       plc%CL3(:,irad)  
              L4(1:2)   =       plc%CL4(:,irad)  + plc%dL4_dT4(:,irad)*(iTemp%T4(1:2,irad,icyc) - plc%Tref_4)
              L5(1:2)   =       plc%CL5(:,irad)  + plc%dL5_dT5(:,irad)*(iTemp%T5(1:2,irad,icyc) - plc%Tref_5)
              LMM(1:2)  =       plc%CLMM(1:2,irad)

              lossFac(1:2,irad,icyc) = LMM * L5 * L4 * L3 * L2A * L2B * L1
           enddo                !irad
        enddo                   !icyc

        return
        end subroutine fe_loss_fac
#endif


!     average two numbers x1 and x2, but check to see if one or both are missing values.
        subroutine avg2(x1,x2,xm, avg)
        implicit none

        real(4), intent(in)  :: x1,x2,xm
        real(4), intent(out) :: avg

        if(abs(x1-xm).gt.0.1 .and. abs(x2-xm).gt.0.1) then !both good,usual case 
           avg =(x1 + x2)/2. 
           return
        endif
         
        if(abs(x1-xm).gt.0.1) then !x1 is there, x2 is missing 
           avg = x1
           return
        endif

        if(abs(x2-xm).gt.0.1) then !x2 is there, x1 is missing
           avg = x2
           return
        endif

        avg=xm                  !both missing, return value indicating missing
        return
        end subroutine avg2


