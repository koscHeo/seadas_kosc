        module l2gen_aquarius_module
        implicit none

        integer(4), parameter :: max_cyc= 5000 ! max number of block in l1a file = (5872 sec/orbit)*(1.2 overlap)/(1.44 sec/block) + slop
        integer(4), parameter :: n_rad=3 ! number of radiometers (i.e., horns)  (inner, middle, outer)  
        integer(4), parameter :: n_prt=85 ! number of thermistors

        integer(4), parameter :: n_subcyc=12 ! number of sub cycles per cycle
        integer(4), parameter :: npol=4 !number of polarizations (v, h, plus, minus)    
        integer(4), parameter :: n_Sacc=6      ! number of S(short) accum (antenna) 
                                                                                   ! There are 8 accumulation for each subcycle
                                                                                   ! Accum 1+2 and accum 3+4 are summed together
                                                                                   ! This leaves a total of 8-2=6 SA
                                                                                   ! Each of SA1 and SA2 are multiplied by 2 

        integer(4), parameter :: n_Lacc=8          ! number of L(long)  accum (calibration)
                                                                                   ! Laccum #5 contains also antenna measurements 
                                                                                   ! There are 4 accum for each subcycle
                                                                                   ! The accum for subcycles 1-10  are summed together
                                                                                   ! The accum for subcycles 11-12 are summed together
                                                                                   ! This leaves 4 LA (each of them multiplied by 10) 
                                                                                   ! and 4 LA (each of them multiplied by 2) 

        integer(2), parameter :: missing_count=-1   ! indicator for missing or bad count
        real(4),    parameter :: missing_val=-9999.     ! indicator for missing value

        integer(4), parameter :: n_gain=42     ! number of cycles used for gain   averaging (= 1 min / 1.44 s)
        integer(4), parameter :: n_off =206    ! number of cycles used for offset averaging (= 5 min / 1.44 s)
        integer(4), parameter :: n_red =14 ! reduced number of average samples (= 20 s / 1.44 s). sed if gain glitch flag =1 

        real(4),    parameter :: delta_phi(n_rad)=(/1.4, 6.6, 4.5/), delta_phi_f(n_rad)=(/1.2, 0.9, 0.9/)! phase imbalances of CND network / OMT subsystem      

!       integer(4), parameter :: ntot=n_subcyc*5  !total number of earth counts per cycle

        real(4),    parameter :: rfi_ta=3.0     ! (kelvin) rfi threshold, if dif between count and medium exceeds this it is flagged

        integer(4), parameter :: nomega=1441, nzang=1441
        integer(4), parameter :: nlon_lnd=2160, nzang_lnd=2881

!     ==================================================================================================
!     ==================================constants for glitch detector ==================================
!     ==================================================================================================

        ! sigma for computing Z in glitch detector
      real(4), parameter, dimension(npol, n_rad) :: sigma_glitch = reshape( 
     1  (/ 0.074, 0.048, 0.060, 0.069,    0.075, 0.048, 0.047, 0.067,    0.060, 0.047, 0.051, 0.073/), 
     2  (/npol, n_rad/) )
 
      ! glitch detection threshold
        real(4), parameter, dimension(npol, n_rad) :: threshold_glitch = 8.0 

      ! glitch detection windows
        integer(4), parameter   :: N1=41, N2=69


!     nominal gain/offsets
!     email from J. Piepmeyer 2/11/2009
        real(4), parameter, dimension(n_rad)   :: gvv_nom= (/1.03, 0.94, 0.93/) 
        real(4), parameter, dimension(n_rad)   :: ghh_nom= (/0.95, 0.91, 1.02/)
        real(4), parameter, dimension(n_rad)   :: gpp_nom= (/0.68, 0.62, 0.60/) 
        real(4), parameter, dimension(n_rad)   :: gmm_nom= (/0.75, 0.68, 0.68/) 
        real(4), parameter, dimension(n_rad)   :: gpv_nom=gpp_nom/2.0, gph_nom=gpp_nom/2.0
        real(4), parameter, dimension(n_rad)   :: gmv_nom=gmm_nom/2.0, gmh_nom=gmm_nom/2.0
        real(4), parameter, dimension(n_rad)   :: gpU_nom=0.40, gmU_nom=-0.40
        real(4), parameter, dimension(n_rad,4) :: goffset_nom=500.0

!     ==================================================================================================
!     ==================================constants for rfi detection ====================================
!     ==================================================================================================

        integer(4), parameter   :: l_data     = n_Sacc*n_subcyc*max_cyc ! 360000
        integer(4), parameter   :: l_ex_data  = 12    *n_subcyc*max_cyc ! 720000
        integer(4), parameter   :: l_data_c   =        n_subcyc*max_cyc ! 6000

! nominal ocean temperature v-pol = 140K, h-pol = 80K AFTER FRONT END LOSSES
!       Band dependent values Paolo De Matthaeis  09/01/12
!        real(4), parameter, dimension(npol) :: ta_ocean_nom = (/170.0, 150.0, 150.0, 130.0/) !VPMH
!        real(4), parameter, dimension(npol*n_rad) :: ta_ocean_nom = (/   ! VPMH Band 123
!     1                      190.0343, 182.7769, 187.0934, 
!     2                      186.6665, 192.0265, 197.7660,
!     3                      181.4487, 185.3383, 188.2318,
!     4                      177.7489, 180.5001, 184.1999 /)

! nominal cnd temperature using 510 K
        real(4), parameter, dimension(npol) :: ta_cnd_nom   = 510.0 

!     NEDT parameters
        real(4), parameter :: receiver_noise=74.6  !Kelvin
        real(4), parameter :: bw_tau=474.34  !sqrt(bandwidth*integration time)=sqrt( 25 MHz * 0.009 sec)


!       RFI aux parameters
        character(len=200), parameter :: rfi_aux_WM='aquarius/radiometer/wm.txt'
        character(len=200), parameter :: rfi_aux_WD='aquarius/radiometer/wd.txt'

        character(len=200), parameter :: rfi_aux_TD='aquarius/radiometer/td.txt'
        character(len=200), parameter :: rfi_aux_TD_A='aquarius/radiometer/TD_asc.txt'
        character(len=200), parameter :: rfi_aux_TD_D='aquarius/radiometer/TD_desc.txt'

        character(len=200), parameter :: rfi_aux_TM='aquarius/radiometer/tm.txt'
        character(len=200), parameter :: rfi_aux_WM_c='aquarius/radiometer/wm_c.txt'
        character(len=200), parameter :: rfi_aux_WD_c='aquarius/radiometer/wd_c.txt'
        character(len=200), parameter :: rfi_aux_TD_c='aquarius/radiometer/td_c.txt'
        character(len=200), parameter :: rfi_aux_TM_c='aquarius/radiometer/tm_c.txt'


!       Prelaunch Calibration Coefficents Structure
        type prelaunch_struct 
!       ref. temp for FE losses and ND
        real(4) :: Tref_3
        real(4) :: Tref_4
        real(4) :: Tref_5

!       ref. temp for FE losses and ND
        real(4) :: Tref_CND
        real(4) :: Tref_ND

        real(4) :: Tref_ND_offset
        real(4) :: Tref_DL_offset

!       noise diodes
        real(4) :: TCND_0(2,n_rad)
        real(4) :: dTCND_dT4(2,n_rad)
        real(4) :: dTCND_dT(2,n_rad)

        real(4) :: TND_0(2,n_rad)
        real(4) :: DTND_DT(2,n_rad)
        real(4) :: TND_offset_0(2,n_rad)
        real(4) :: DTND_offset_DT(2,n_rad)
        real(4) :: TDL_offset_0(2,n_rad)
        real(4) :: dTDL_offset_DT(2,n_rad)

        real(4) :: CL1(2,n_rad)
        real(4) :: CL2A(2,n_rad)
        real(4) :: CL2B(2,n_rad)
        real(4) :: CL3(2,n_rad)
        real(4) :: CL4(2,n_rad)
        real(4) :: dL4_dT4(2,n_rad)
        real(4) :: CL5(2,n_rad)
        real(4) :: dL5_dT5(2,n_rad)
        real(4) :: CLMM(2,n_rad)

!       non-linearity correction coefficients
        real(4) :: Dnl20(4,n_rad)
        real(4) :: Dnl21(4,n_rad)
        real(4) :: Dnl22(4,n_rad)
        real(4) :: Dnl30(4,n_rad)
        real(4) :: Dnl31(4,n_rad)
        real(4) :: Dnl32(4,n_rad)
        real(4) :: Tref_nl(4,n_rad)
        end type


!       Instrument Temperature Structure

!       T1, T2A, T2B are zones 1, 2A/2B temperatures are polarization independent
!       T3, T4, T5   are zones 3,4,5 Temperatures can be different for V/H. 1=V, 2=H
!       TDL is the un-corrected Dicke load temperature [kelvin] (V/H)
!       T0 is dicke load reference temperatures [kelvin] (V/H)
!       TDL_offset is the offset correction for the Dicke load [kelvin] (V/H)
!       TCND_P is physical temperatures of CND [kelvin] (V/H)
!       TND_P is physical temperatures of ND [kelvin]  (V/H)
!       TND is noise diode injection temperature per cycle [kelvin] (V/H)
!       TND_offset is the offset correction for the ND [kelvin] (V/H) 
!       TCND is CND injection temperature per cycle [kelvin] (V/H)         
!       Tdet is detector temperatures [kelvin] (1=V, 2=P, 3=M, 4=H)
        type instrument_temp_struct 

        real(4) :: T1 (  n_rad,max_cyc)
        real(4) :: T2A(  n_rad,max_cyc)
        real(4) :: T2B(  n_rad,max_cyc)

        real(4) :: T3(2,n_rad,max_cyc)
        real(4) :: T4(2,n_rad,max_cyc)
        real(4) :: T5(2,n_rad,max_cyc)
        real(4) :: T0(2,n_rad,max_cyc)
        real(4) :: TCND_P(2,n_rad,max_cyc)
        real(4) :: TND_P( 2,n_rad,max_cyc)
        real(4) :: TCND(  2,n_rad,max_cyc)
        real(4) :: TND(   2,n_rad,max_cyc)
        real(4) :: TND_offset(   2,n_rad,max_cyc)

        real(4) :: TDL(  2,n_rad,max_cyc)
        real(4) :: TDL_offset(   2,n_rad,max_cyc)

        real(4) :: Tdet(4,n_rad,max_cyc)

        end type


!       Gain/Offset structure   
        type gain_off_struct 

        real(4) :: gvv(n_rad,max_cyc)
        real(4) :: ghh(n_rad,max_cyc)
        real(4) :: ov( n_rad,max_cyc)
        real(4) :: oh( n_rad,max_cyc)
        real(4) :: op( n_rad,max_cyc)
        real(4) :: om( n_rad,max_cyc)
        real(4) :: gpv(n_rad,max_cyc)
        real(4) :: gph(n_rad,max_cyc)
        real(4) :: gmv(n_rad,max_cyc)
        real(4) :: gmh(n_rad,max_cyc)
        real(4) :: gpp(n_rad,max_cyc)
        real(4) :: gmm(n_rad,max_cyc)
        real(4) :: gpU(n_rad,max_cyc)
        real(4) :: gmU(n_rad,max_cyc)
        real(4) :: hpU(n_rad,max_cyc)
        real(4) :: hmU(n_rad,max_cyc)
! gpU, gmU = with RFI filter, hpU, hmU = without RFI filter

        end type


!       RFI structure   
        type rfi_struct 
!       1 means this earth count has rfi
        integer(1) :: iflag_rfi(n_Sacc, n_subcyc, npol, n_rad, max_cyc)

!       total number of rfi earth counts
        integer(4) :: num_rfi(                npol, n_rad, max_cyc)

!       1 means this CND count has rfi
        integer(1) :: iflag_rfi_CND(n_subcyc, npol, n_rad, max_cyc)

!       1 means this LA is flagged for gain glitch
        integer(1) :: iflag_glitch(           npol, n_rad, max_cyc)

!       RFI aux arrays

        integer(4), dimension(181,360)          :: wm
        integer(4), dimension(181,360)          :: wd
        real(4), dimension(181,360)             :: tm
        real(4), dimension(181,360)             :: td
        real(4), dimension(181,360)             :: td_a
        real(4), dimension(181,360)             :: td_d

        integer(4), dimension(181,360)          :: wm_c
        integer(4), dimension(181,360)          :: wd_c
        real(4), dimension(181,360)             :: tm_c
        real(4), dimension(181,360)             :: td_c

        real(4), dimension(181,360,npol,n_rad)  :: stdta_rad
        real(4), dimension(npol,n_rad)          :: stdta_cnd

        integer(4)                              :: idata

        end type


        end module l2gen_aquarius_module



