subroutine fd_dtb_roughness (itype, iwdir, irad, wspd, phir, sigv, swh, sst, sss, staticdata, ioob, dtb_rough, delta_tb_sst_wspd)
!     TM June 2013: added sst and sss to interface with fd_TM_emiss_harmonics     
!     TM September 25 2012: The 2-dimensional array is now the residual roughness depending on wind speed and sigma0_VV
!     DTB_rough = aharm(0) + DTB_rough_residual(W,sigma0_VV) + DTB_rough_residual(W,SWH)
!
!     TM August 10, 2012:
!     In V7 dtb_rough, dtb_phi and aharm are now all emissivities [*teff = 290K]. 
!     No change in this routine is necessary, but the routine fd_tbsur_sss needs to be changed accordingly 
use static_data_module

! wind induced emissivity from wind speed + scatterometer sigma0
implicit none

integer(4), intent(in)                                  ::  itype 
! 0=wspd only  1=wspd+sigma0_VV  2=wspd+sigma0_VV+swh    
        
integer(4), intent(in)                                  ::  iwdir 
! 0=no wdir correction is performed   1= wdir correction is performed

integer(4), intent(in)                                  ::  irad
real(4), intent(in)                                             ::  wspd
real(4), intent(in)                                             ::  phir
real(4), intent(in), optional                   ::  sigv        ! sigma0 VV in real units before wdir signal removed
real(4), intent(in), optional                   ::  swh         ! significant wave height [m]
real(4), intent(in), optional                   ::  sst         ! SST [C]
real(4), intent(in), optional                   ::  sss         ! salinity [psu]
        
type(static_data_struct), intent(in)    :: staticdata

integer(4), intent(inout)          ::  ioob ! OOB flag  
! 0= can use retrieval    1=underpopulated bin. do not use
real(4), dimension(2), intent(out) ::  dtb_rough  !(1=V, 2=H) 
real(4), dimension(2), intent(out) ::  delta_tb_sst_wspd  !(1=V, 2=H) 
        
real(4), dimension(0:2,2)          :: aharm  !1=V,  2=H
real(4), dimension(0:2,npol_sc)    :: bharm  !1=VV, 2=HH,  3=VH,  4=HV  
real(4), dimension(0:2,2)          :: daharm  
real(4), dimension(0:2,npol_sc)    :: dbharm   

integer(4)                         ::  ibinx1, ibiny1
integer(4)                         ::  ibinx2, ibiny2
real(4)                            ::  xpoint, ypoint
real(4)                            ::  xvec0, yvec0
real(4)                            ::  dx, dy
real(4)                            ::  x1, y1
real(4)                            ::  x2, y2
real(4)                            ::  brief_x, brief_y

real(4), dimension(2)              ::  c1, c2, z1, z2
real(4), dimension(2)              ::  xtbw_scat, xtbw_wav  !(1=V, 2=H) 

integer(4)                         ::  ioob1, ioob2

real(4)                            ::  xsigv ! sigma0 in real units after wdir signal removed
real(4), dimension(2)              ::  dtb_phi, dsigma_phi, dtb_isotropic, demiss_sst_wspd

real(4) cosd
integer, pointer :: nptr => null()

real(4), dimension(2)              ::  dew_1, dew_r

real(4), parameter :: teff=290. ! used for wind induced emissivity harmonics

if (itype>=1 .and. .not.(present(sigv))) stop ' need sigv in fd_dtb_roughness'
if (itype==2 .and. .not.(present(swh)))  stop ' need swh  in fd_dtb_roughness'

!     changed 02/03/2015 for V3.3 
!     include DTB(SST, WSPD):
!     there was a little bug that I fixed on 2/4/2015
!     I do NOT want to transform the surface emsissivity into a TB yet because that will be done 
!     in subroutine fd_tbsur_sss in the block dtb_roughness_corr(1:2,irad,icyc) = dtbc(:)
if (staticdata%W_STITCH.ne.-1) then
   call fd_TM_emiss_harmonics (irad, wspd, %VAL(0), %VAL(0), staticdata, aharm, daharm)
   dew_1(1:2) = aharm(0,1:2)/teff
   call fd_TM_emiss_harmonics (irad, staticdata%W_STITCH, %VAL(0), %VAL(0), staticdata, aharm, daharm)
   dew_r(1:2) = aharm(0,1:2)/teff

   call fd_demiss_sst_wspd (irad, sst, staticdata, dew_1, dew_r, demiss_sst_wspd)    
else
   demiss_sst_wspd = 0.0
endif

if (loc(delta_tb_sst_wspd).ne.0) then
   delta_tb_sst_wspd(1:2) = demiss_sst_wspd(1:2) * (sst+273.15) / teff
endif

if (staticdata%emiss_sst_sss .eq. 1) then
   call fd_TM_emiss_harmonics (irad, wspd, sst, sss, staticdata, aharm, daharm)
else
   call fd_TM_emiss_harmonics (irad, wspd, nptr, nptr, staticdata, aharm, daharm)
endif

call fd_TM_scat_harmonics (irad, wspd, staticdata, bharm, dbharm)

dtb_isotropic(1:2) = aharm(0,1:2)

if (iwdir==1) then
   dtb_phi(1:2)    = aharm(1,1:2)*cosd(phir) + aharm(2,1:2)*cosd(2.0*phir)         
   dsigma_phi(1:2) = bharm(1,1:2)*cosd(phir) + bharm(2,1:2)*cosd(2.0*phir)
else if (iwdir==0) then ! isotropic
   dtb_phi=0.0
   dsigma_phi=0.0
endif


if (itype==0) then  ! use wind speed only
   ioob=0
   dtb_rough(:) = aharm(0,:) + dtb_phi(:)
   return
endif
        
xsigv = sigv - dsigma_phi(1) 

ioob1  = 0
xpoint = wspd
ypoint = xsigv
xvec0  = dwin/2.0  
yvec0  = staticdata%dsigma(1,irad)/2.0
dx     = dwin  
dy     = staticdata%dsigma(1,irad)
        
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


if (staticdata%iflag_dtbw(1,1,irad,ibinx1,ibiny1)==0 .and. staticdata%iflag_dtbw(1,1,irad,ibinx2,ibiny1)==0) then
   z1 = staticdata%arr_dtbw(1, :, irad,ibinx1,ibiny1)
   z2 = staticdata%arr_dtbw(1, :, irad,ibinx2,ibiny1)
else if (staticdata%iflag_dtbw(1,1,irad,ibinx1,ibiny1)==0 .and. staticdata%iflag_dtbw(1,1,irad,ibinx2,ibiny1)==1) then
   z1 = staticdata%arr_dtbw(1, :, irad,ibinx1,ibiny1)
   z2 = z1
else if (staticdata%iflag_dtbw(1,1,irad,ibinx1,ibiny1)==1 .and. staticdata%iflag_dtbw(1,1,irad,ibinx2,ibiny1)==0) then
   z1 = staticdata%arr_dtbw(1, :, irad,ibinx2,ibiny1)
   z2 = z1
else
   ioob=1
   dtb_rough=missing_val
   return
endif

c1 = z1*(1.0-brief_x) + z2*brief_x


if (staticdata%iflag_dtbw(1,1,irad,ibinx1,ibiny2)==0 .and. staticdata%iflag_dtbw(1,1,irad,ibinx2,ibiny2)==0) then
   z1 = staticdata%arr_dtbw(1, :, irad,ibinx1,ibiny2)
   z2 = staticdata%arr_dtbw(1, :, irad,ibinx2,ibiny2)
else if (staticdata%iflag_dtbw(1,1,irad,ibinx1,ibiny2)==0 .and. staticdata%iflag_dtbw(1,1,irad,ibinx2,ibiny2)==1) then
   z1 = staticdata%arr_dtbw(1, :, irad,ibinx1,ibiny2)
   z2 = z1
else if (staticdata%iflag_dtbw(1,1,irad,ibinx1,ibiny2)==1 .and. staticdata%iflag_dtbw(1,1,irad,ibinx2,ibiny2)==0) then
   z1 = staticdata%arr_dtbw(1, :, irad,ibinx2,ibiny2)
   z2 = z1
else
   ioob1=1
   ioob=1
   dtb_rough=missing_val
   return
endif
c2 = z1*(1.0-brief_x) + z2*brief_x
xtbw_scat = c1*(1.0-brief_y) + c2*brief_y
        
if (itype ==1) then
   if (ioob == -9999) then
      ! For version 2.0 back-compatibility
      dtb_rough = dtb_phi + xtbw_scat
   else
      dtb_rough = dtb_isotropic + dtb_phi + xtbw_scat
   endif
   ioob1=0
   ioob2=0
   ioob=0
   return
endif

ioob2  = 0
xpoint = wspd
ypoint = swh
xvec0  = win0x + dwinx/2.0  
yvec0  = wav0x + dwavx/2.0
dx     = dwinx  
dy     = dwavx
        
ibinx1 = floor((xpoint-xvec0)/dx) + 1
if (ibinx1<1)           ibinx1=1
if (ibinx1>nbin_winx-1) ibinx1=nbin_winx-1      
ibinx2=ibinx1+1
x1 = xvec0 + (ibinx1-1)*dx
x2 = x1 + dx
brief_x = (xpoint-x1)/(x2-x1)

ibiny1 = floor((ypoint-yvec0)/dy) + 1
if (ibiny1<1)           ibiny1=1
if (ibiny1>nbin_wavx-1) ibiny1=nbin_wavx-1      
ibiny2=ibiny1+1
y1 = yvec0 + (ibiny1-1)*dy
y2 = y1 + dy
brief_y = (ypoint-y1)/(y2-y1)

if (staticdata%iflag_winx_wavx(1,irad,ibinx1,ibiny1)==0 .and. staticdata%iflag_winx_wavx(1,irad,ibinx2,ibiny1)==0) then
   z1 = staticdata%arr_dtbwx(:, irad,ibinx1,ibiny1)
   z2 = staticdata%arr_dtbwx(:, irad,ibinx2,ibiny1)
else if (staticdata%iflag_winx_wavx(1,irad,ibinx1,ibiny1)==0 .and. staticdata%iflag_winx_wavx(1,irad,ibinx2,ibiny1)==1) then
   z1 = staticdata%arr_dtbwx(:, irad,ibinx1,ibiny1)
   z2 = z1
else if (staticdata%iflag_winx_wavx(1,irad,ibinx1,ibiny1)==1 .and. staticdata%iflag_winx_wavx(1,irad,ibinx2,ibiny1)==0) then
   z1 = staticdata%arr_dtbwx(:, irad,ibinx2,ibiny1)
   z2 = z1
else
   ioob2=1
   ioob=1
   dtb_rough=missing_val
   return
endif
c1 = z1*(1.0-brief_x) + z2*brief_x
                
if (staticdata%iflag_winx_wavx(1,irad,ibinx1,ibiny2)==0 .and. staticdata%iflag_winx_wavx(1,irad,ibinx2,ibiny2)==0) then
   z1 = staticdata%arr_dtbwx(:, irad,ibinx1,ibiny2)
   z2 = staticdata%arr_dtbwx(:, irad,ibinx2,ibiny2)
else if (staticdata%iflag_winx_wavx(1,irad,ibinx1,ibiny2)==0 .and. staticdata%iflag_winx_wavx(1,irad,ibinx2,ibiny2)==1) then
   z1 = staticdata%arr_dtbwx(:, irad,ibinx1,ibiny2)
   z2 = z1
else if (staticdata%iflag_winx_wavx(1,irad,ibinx1,ibiny2)==1 .and. staticdata%iflag_winx_wavx(1,irad,ibinx2,ibiny2)==0) then
   z1 = staticdata%arr_dtbwx(:, irad,ibinx2,ibiny2)
   z2 = z1
else
   ioob2=1
   ioob=1
   dtb_rough=missing_val
   return
endif
c2 = z1*(1.0-brief_x) + z2*brief_x
xtbw_wav = c1*(1.0-brief_y) + c2*brief_y
        
! changed January 2015: include DTB(SST,WSPD)                
dtb_rough = dtb_isotropic + dtb_phi + xtbw_scat + xtbw_wav + demiss_sst_wspd
ioob1=0
ioob2=0
ioob=0 

    
return

end subroutine fd_dtb_roughness



! TM 01/26/2015
! changed 01/30/2015
! calculate correction demiss
! demiss (sst, wspd) = delta(sst) + beta(sst)*delta_EW(wspd)!
subroutine fd_demiss_sst_wspd (irad, sst, staticdata, dew_1, dew_r, demiss_sst_wspd)
use static_data_module
implicit none

integer(4), intent(in)                              ::  irad
real(4), intent(in)                                 ::  sst  ! Celsius
type(static_data_struct), intent(in)                ::  staticdata
real(4), intent(in), dimension(2)                   ::  dew_1, dew_r 
real(4), intent(out), dimension(2)                  ::  demiss_sst_wspd

real(4)                                             ::  xsurtep
real(4), dimension(2)                               ::  delta
integer(4)                                          ::  ipol

integer(4)                                          ::  ksst1, ksst2
real(4)                                             ::  x1, x2, brief, y1, y2

real(4), parameter :: teff=290. ! used for wind induced emissivity harmonics

xsurtep = sst+273.15

ksst1 = floor((sst - staticdata%sst0)/staticdata%sst_step) + 1
if (ksst1< 1)           ksst1 = 1
if (ksst1> nsst_arr-1)  ksst1 = nsst_arr-1
ksst2 = ksst1 + 1

do ipol=1,2
    
    x1 = staticdata%sst0 + (ksst1-1)*staticdata%sst_step
    x2 = x1   + staticdata%sst_step
    brief = (sst-x1)/staticdata%sst_step 
    y1 = staticdata%yarr(ipol,irad,ksst1)
    y2 = staticdata%yarr(ipol,irad,ksst2)
    delta(ipol) = y1*(1.0-brief) + y2*brief
    
    ! fan out as function of wind speed above T_STITCH
    if (sst <= staticdata%T_STITCH(ipol,irad)) then
        demiss_sst_wspd(ipol) = delta(ipol)
    else
        demiss_sst_wspd(ipol) = delta(ipol) * dew_1(ipol)/dew_r(ipol)
    endif

enddo ! ipol


! transform the dtb into an emissivity
demiss_sst_wspd = teff*demiss_sst_wspd/xsurtep

return
end subroutine fd_demiss_sst_wspd 

