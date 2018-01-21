      module static_data_module
      implicit none

      integer(4), parameter :: n_rad=3 ! number of radiometers (i.e., horns)
      integer(4), parameter :: n_zang=1441
      integer(4), parameter :: n_omega=1441
      integer(4), parameter :: nlon_lnd=2160
      integer(4), parameter :: nzang_lnd=2881
      integer(4), parameter :: nbin_w=60
      integer(4), parameter :: mbin_w=60
      integer(4), parameter :: npol_sc=4     
      integer(4), parameter :: nwin_sc=40     
      integer(4), parameter :: nbin_winx=40
      integer(4), parameter :: nbin_wavx=40
      integer(4), parameter :: mbin_winx=40
      integer(4), parameter :: mbin_wavx=40
      real(4),    parameter :: missing_val=-9999. ! indicator for missing value

      real(4),    parameter :: win0x=0.0
      real(4),    parameter :: wav0x=0.0
      real(4),    parameter :: dwinx=1.0
      real(4),    parameter :: dwavx=0.5
      real(4),    parameter :: dwin=1.0

      integer(4), parameter :: isstmax=30, nsst_arr=isstmax*100+1

      type static_data_struct
      real(4) :: apc_matrix(3,3,n_rad)
      real(4) :: apc_inverse(3,3,n_rad)
      real(4) :: sss_coef(4,451,251)
      real(8) :: time_sun(n_zang)
      real(4) :: eia_sun(n_zang,n_rad)
      real(4) :: tasun_dir_tab(n_omega,n_zang,3,n_rad)
      real(4) :: tasun_ref_tab(n_omega,n_zang,3,n_rad)
      real(4) :: tasun_bak_tab(161,26,3,n_rad)
      real(8) :: time_galaxy(n_zang)
      real(4) :: eia_galaxy(n_zang,n_rad)
      real(4) :: tagal_dir_tab(n_omega,n_zang,3,n_rad)
      real(4) :: tagal_ref_tab(n_omega,n_zang,3,n_rad,5)
      real(8) :: dtagal_ref_tab(2,n_rad,n_omega,n_zang)
      real(4) :: fpt_lnd(nlon_lnd,nzang_lnd,n_rad)
      real(4) :: frc_lnd(nlon_lnd,nzang_lnd,n_rad)
      integer(2) :: landcorr(1440,1440,12,2,n_rad)
      real(4) :: gain_ice(151,n_rad)
      integer(2) :: itau(360,180,12,91)
      real(8) :: acoef(0:2,2,3,5)
      real(8) :: bcoef(0:2,4,3,5)
      real(8) :: arr_dtbw(3,2,n_rad,nbin_w,nbin_w)
      real(8) :: arr_dtbwx(2,n_rad,nbin_winx,nbin_wavx)
      real(4) :: dsigma(2,n_rad)
      integer(1) :: iflag_dtbw(3,2,n_rad,nbin_w,nbin_w)
      integer(1) :: iflag_winx_wavx(2,n_rad,nbin_winx,nbin_wavx)
      real(8) :: wspd_max_a(0:2,2,n_rad) ! high wind speed for radiometer wind speed signal	
      real(8) :: wspd_max_b(0:2,npol_sc,n_rad) ! high wind speed for scatterometer sigma0
      real(4) :: coeffs_dI_U(4,n_rad)

      real(4) ::  sst_step, sst0, sst1
      real(8) ::  yarr(2,n_rad,nsst_arr)
      real(8) ::  T_STITCH(2,n_rad) ! computed by IDL routine. delta(T_STITCH) = 0 
      real(4) ::  W_STITCH ! average wind speed (7 m/s)

      real(4) :: estimated_error_array(7,nwin_sc)
      integer(2) :: climate_sal(360,180,12)
      integer(1) :: rfi_mask(2,180,90)
      integer(1) :: emiss_sst_sss
      end type

      end module static_data_module


