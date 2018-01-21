#include "l12_parms.h"
#undef PI

C ======================================================================
C
C             P R O J E C T:  S P E C O P T I M - C A S E 2
C              ---------------------------------------------
C
C ----------------------------------------------------------------------
C
C Purpose:       Atmospheric Correction of Sattelite Color Imagery
C                over the Ocean in the presence of absorbing
C                and/or non-absorbing aerosol
C
C ----------------------------------------------------------------------
C 
C $Date: 2005/05/11 $
C $Authors:  R.Chomko,C.P.Kuchinke
C
C ----------------------------------------------------------------------
C Version:      4.0
C CVS-Tag:
C Computer:     DELL LINUX
C Precision:    SINGLE and DOUBLE
C Language:     FORTRAN-77
C ----------------------------------------------------------------------
C Usage:
C Arguments:
C References:
C ----------------------------------------------------------------------
C Copyright:    (c) 2000-2002, R.Chomko & H.R.Gordon, U.Miami-Physics
C Copyright:    (c) 2003-2007, C.P.Kuchinke & H.R.Gordon, U.Miami-Physics
C Distribution: Ltd
C Warranty:
C ----------------------------------------------------------------------
C
C ---------------------------- OVERVIEW --------------------------------

C Program starts with interface subroutine ATMCOR_SOA that calls main controlling subroutine 
C SPEC1GEOFIT.  After this, program is divided into 10 sections.

C SECTION 1: NIR optimization (using ZXMWD/E) to derivestimate functions nu(mi,mr) and tau(mi,mr)

C SECTION 2: Spectral optimization to derive all parameters

C SECTION 3: Data set up routines required by SPEC1GEOFIT and NIR driving routine 
C            DRV_TAUVV_ZXMWD in SECTION 1 (NIR optimization of nu,tau865 for given mi,mr 
C            using ZXMWD/E)

C SECTION 4: Data set up routines required by DRV_SPECOPT_LBFGSB (main optimization),
C            DFCNLSQ_GRAD (SD and gradiant updates of main optimization) and SPECOPT_POSTPROC
C            (post-processing of main-optimization data), all contained in SECTION 2.     

C SECTION 5: Forward WATER-model routines

C SECTION 6: Partial derivative routines

C SECTION 7: Atmosphere look-up and and water initialisation data

C SECTION 8: Radiometric correction routines.

C SECTION 9: Core optimization routines for NIR (ZXMWE, ZXMWD etc.)

C SECTION 10: Core optimization routines for main (vis) = Setulb, mainlb etc.

C **********************************************************************
C(rc) Note for SeaDAS:
C                  In C-routines pixel numbering starts from 0.
C                  It is converted in here to (pixel+1) since array
C                  numbering in Fortran begins from 1. Remember to
C                  to convert it back within soa_sma_cutils.c in the
C                  function 'get_gs97_um', otherwise everything will
C                  be shifted by 1 pixel to the left.
C                                                    
C **********************************************************************

       SUBROUTINE ATMCOR_SOA (
     &       sensorID, sensorNm, nwave_in, wave_in,
     &       scan,pixel,
     &       theta0,theta,delt_phi,
     &       lat_in,lon_in,
     &       totRad,rayRho,F0ln,Taur_in,
     &       aw_in, bbw_in,
     &       aphstar_in, adg_s_in, bbp_s_in,
     &       resAerRho, resWatRho,
     &       water_vapor,t_star,t_star0,
     &       taua,
     &       optW0,optPIG, optACDM, pcentcdm,
     &       optB0, optMR, optMI,
     &       intpNu,ierror)
c--input: 
c         scan           -current line number
c         pixel          -current pixel number
c         theta0         -solar zenith (deg)
c         theta          -sensor zenith (deg)
c         delt_phi       -sensor relative azimuth (deg)
c         lat/lon        -only used for Ches Bay Bio-optical coeff selection
c         totRad         -total observed radiance = pseudo rhot (mW/cm^2/um/Sr)
c                        *whitecap radiance is removed (see atmocor2.f)
c                        *corrected for ozone absorption (see atmocor2.f)
c         rayRho         -Rayleigh reflectance
c         F0ln           -E-S distance corrected Solar irradiance
c--output:
c         resAerRho      -resultant aerosol reflectance rhoA (= aerosol reflectance, prime)
c         resWatRho      -resultant water reflectance rhow (= water reflectance at surface, prime)
c         water_vapor    -water vapor concentration g/cm^2
c         t_star         -diffuse transmittance outbound (at theta)
c         t_star0        -diffuse transmittance inbound (at theta0)
c         taua           -retrieved aerosol optical thickness
c         optW0          -optimized single scattering albedo
c         optPIG         -optimized chlorophyll concentration (mg/l)
c         optACDM        -optimized cdm absorption at 443-nm (= absorp.gelb + absorp.det)
c         pcentcdm       -% CDM: ratio acdm(443) to atot(443); aw(443) not included
c         optB0          -optimized particulate backscattering coeff.at 443-nm
c         optMR          -optimized real part of refractive index
c         optMI          -optimized imaginary part of refractive index
c         intpNu         -optimized aerosol Junge distribution parameter
c         ierror         -error status for outside world; 0 = OK, 1 = negative or no solution

      INTEGER    NRAD,MPHI,NUM
      PARAMETER (NRAD=100,MPHI=15,NUM=75)

C ------ INPUT/PARAMETERS=SPECMATCH_GENERIC

      INTEGER    sensorID         ! sensor ID number 
      CHARACTER  sensorNm*20      ! sensor name (string)
      INTEGER    nwave_in         ! number of wavelengths
      REAL       wave_in(*)       ! list of wavelengths (nm)
      REAL       aw_in(*)         ! water absorption coeff (per wavelength)
      REAL       bbw_in(*)        ! water backscatter coeff (per wavelength)
      REAL       Taur_in(*)       ! Rayeliegh optical thickness (per wavelength)
      REAL       aphstar_in(*)    ! GSM parameter aph* (phyto specific absorp per wavelength )
      REAL       adg_s_in         ! GSM parameter S    (gelbstof absorption slope)
      REAL       bbp_s_in         ! GSM parameter eta  (particle backscatter slope)

      INTEGER    NMI,NMR,NV,
     &           NLAMBDA,NLAMBD0,
     &           NBEST,NBand,IERopt
      PARAMETER (NMI=6,NMR=2,NV=6)
      PARAMETER (NLAMBDA=8,NLAMBD0=2)
      PARAMETER (NBEST=8)

C ------ SEADAS-SPECIFIC

      INTEGER    NDIR,NMODEL,NSUN
      PARAMETER (NDIR=35,NMODEL=NV*NMR*NMI,NSUN=33)

C ------ SCAN/LINE & PIXEL VARS
c    spixel controls which pixel of a particular line is being processed

      INTEGER    pixel,scan,spixel,pix_to_screen,scan1,pixel1 

C ------ Geom/season vars

      REAL       lat_in,lon_in
      INTEGER    iseason 

C ------ INTERNAL/MISC

      REAL       air_mass
      INTEGER    ierror
      INTEGER    wavelength(NLAMBDA)
      INTEGER    i412,i443,i490,i510,i555,i670,i765,i865
      PARAMETER (i412=1,i443=2,i490=3,i510=4,i555=5,i670=6,
     &             i765=7,i865=8)
      
      REAL       mu0, index
     
C ------ RAYLEIGH PARAMETER DATA 

      REAL       tauray(NLAMBDA),rho_wat_giv(NLAMBDA)

C ------ Haze-C MODEL PARAMETER DATA

C       aerW0     Single Scat Albedo W0
C       aerC      Extinction Coeff   C
C       s11       S11 elements of Scat Muller Matrix
C       anglPhs   scat angles (not used)

      REAL       mi(NMI),mr(NMR),v(NV)

C    Look-Up Reflectance LUT
      REAL    luttau_a(NDIR,MPHI,NSUN,NLAMBDA,NMODEL),
     &        luttau_b(NDIR,MPHI,NSUN,NLAMBDA,NMODEL),
     &        luttau_c(NDIR,MPHI,NSUN,NLAMBDA,NMODEL),
     &        luttau_d(NDIR,MPHI,NSUN,NLAMBDA,NMODEL)
C    Look-Up Diffuse Transmittance LUT
      REAL    diffTrans_a(NSUN,NLAMBDA,NMODEL),
     &        diffTrans_b(NSUN,NLAMBDA,NMODEL),
     &        diffA_ri(NMR,NMI,NLAMBDA),               ! theta, for interpolated V
     &        diffB_ri(NMR,NMI,NLAMBDA),
     &        diffA0_ri(NMR,NMI,NLAMBDA),              ! theta0, for interpolated V
     &        diffB0_ri(NMR,NMI,NLAMBDA)
    
cdas 'aerW0' now becomes a matrix of W0 values for all models
      REAL    aerW0(NLAMBDA,NMI,NMR,NV), aerC(NLAMBDA,NMI,NMR,NV),
     &        cc865tab(NLAMBDA,NMI,NMR,NV),cc865(NMI,NMR,NLAMBDA),
     &        tau865(NMI,NMR),
     &        s11dummy(NUM)

C -- WATER PARAMETER DATA --------------------------------------------------

c------output
      REAL    totWatRad(NLAMBDA)
      REAL    tla(NLAMBDA)
      REAL    opt_nRhoW(NLAMBDA)
      REAL    absCoef(NLAMBDA),abs_out(NLAMBDA)

C-INPUT: GEOMETRY DATA --------------------------------------------------

      REAL    sunTAB(NSUN),theta0,theta,delt_phi

C ------- INPUT: RADIANCE DATA --------------------------------------------------

      INTEGER out_refl
      REAL    F0ln(NLAMBDA), F0_c(NLAMBDA),mu0_c,pi_c
      REAL    water_vapor
      REAL    rhot(NLAMBDA),aerwatRho(NLAMBDA),aerwatRho_c(NLAMBDA),
     &        totRad(NLAMBDA),rayRho(NLAMBDA)
      REAL    FUNCT_OXYGEN_AER

      REAL    candRho(NLAMBDA,NMI,NMR)
      REAL    rcandRho(NMI,NMR,NLAMBDA)

      INTEGER vstep,ivdummy,jidum,jrdum
      REAL    intpNu,vv_c

C ------ RETRIEVED OPTIMAL VLAUES

      REAL     optSTDDEV,optMI,optMR,optPIG,optACDM,optB0,
     &         optTAU,optW0,
     &         optTAUL(NLAMBDA),optTrans(NLAMBDA),optTrans0(NLAMBDA)

      REAL     optOC4,pcentcdm

      REAL    resAerRho(NLAMBDA),resWatRho(NLAMBDA),   ! Retrieved Reflectances
     &        resTotRho(NLAMBDA)

      REAL    t_star(NLAMBDA),t_star0(NLAMBDA)         ! transm at theta, theta0

      REAL    taua(NLAMBDA)

      REAL    fl(NLAMBDA)

C ------ SAVING VARIABLES THAT SHOULD BE FOUND ON START-UP ONLY

      LOGICAL  FIRSTB
      SAVE     FIRSTB

      REAL*4 adm_s,bbp_s,aw(NLAMBDA),bbw(NLAMBDA),aph_lin(NLAMBDA)
      COMMON /WATERPAR/ adm_s,bbp_s,aw,bbw,aph_lin

      REAL*8 d_adm_s,d_bbp_s,
     &       d_aw(NLAMBDA),d_bbw(NLAMBDA),d_aph_lin(NLAMBDA)
      COMMON /D_WATERPAR/ d_adm_s,d_bbp_s,d_aw,d_bbw,d_aph_lin

c      REAL     OC4_CHL4S
c      EXTERNAL OC4_CHL4S
C---------------------------------------------------------------------------------------------------------
C     DO NOT DELETE THESE LINES. External commons for Seadas structures.  
C     These COMMON structures will not be found anywhere anywhere else in SOA.  
C     They are commons external to this code. After msl12 finishes with a line 
C     these arrays contain all the required info for all the pixels, and 
C     everything's ready for an output (see soa_sma_utils.f). 
      INTEGER NSCANSZ
      PARAMETER (NSCANSZ=MAXPIX)
      REAL*4  chl_out(NSCANSZ),acdm_out(NSCANSZ),bbp_out(NSCANSZ),
     &        chl4s_out(NSCANSZ),pcentcdm_out(NSCANSZ)
      COMMON /GS97_OUT/ chl_out,acdm_out,bbp_out,chl4s_out,
     &                  pcentcdm_out
      REAL*4  w0_out(NSCANSZ)
      COMMON /GENERIC_OUT/ w0_out
      REAL*4  v_out(NSCANSZ)
      COMMON /SOA_OUT/ v_out
C---------------------------------------------------------------------------------------------------------
C -- COMMON
      COMMON /COMPHASE/    aerW0,aerC,cc865tab,s11dummy
      COMMON /COMRTELUT/   luttau_a,luttau_b,luttau_c,luttau_d
      COMMON /COMDIFF/     diffTrans_a,diffTrans_b

      COMMON /BAND/        NBand
      COMMON /OPTFUNCREAL/ mi,mr,cc865,tau865,
     &                     diffA_ri,diffB_ri,
     &                     diffA0_ri,diffB0_ri,
     &                     aerwatRho,candRho,
     &                     rcandRho
      COMMON /SMM_GENERIC/ wavelength
      COMMON /TAUVV/       jrdum,jidum,aerwatRho_c
      COMMON /TAUVV_MISC/  tauray,rho_wat_giv
      COMMON /V_MISC/      v,vv_c,vstep,ivdummy
      COMMON /SPEC_MISC/   F0_c,mu0_c,pi_c
      COMMON /STDDEV/      NuTmp,W0Tmp,PigTmp,TauTmp

      REAL   rhow_diff(0:4)
      COMMON /DIFF_OUT/    rhow_diff 

      COMMON /pixel_proc/  scan1,pixel1 

C -- DATA BLOCKS

      DATA FIRSTB /.TRUE./

      DATA sunTAB /
     &    0.0,  2.5,  5.0,  7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5,
     &   25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5,
     &   50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5,
     &   75.0, 77.5, 80.0                                           /

      DATA  tauray /0.3132, 0.2336, 0.1547, 0.1330,
     &              0.0947, 0.0446, 0.0256, 0.0169/

      DATA  wavelength /
     &              412, 443, 490, 510,
     &              555, 670, 765, 865 /
      
      DATA  mi /
     &              0.000, 0.001, 0.003,
     &              0.010, 0.030, 0.040 /
    
      DATA  mr /    1.33, 1.50/

      DATA  v  /    2.0, 2.5, 3.0, 3.5, 4.0, 4.5/

      DATA  NBand  / 6 /

      rad(x) = x*pi/180.0

CBAF
      if (sensorID .ne. SEAWIFS) then
          write(*,*) 'Sensor is ',sensorNm,' ID is',sensorID
          write(*,*) 'adg_s = ',adg_s_in
          write(*,*) 'bbp_s = ',bbp_s_in
          write(*,*)
          write(*,1) 'Bnd','Wave','aph*','aw','bbw','Fo','Tau_r'
          do j=1,nwave_in
              write(*,2) j,wave_in(j),aphstar_in(j),aw_in(j),bbw_in(j),
     &                   F0ln(j),Taur_in(j)
          enddo
 1        FORMAT( a3,x,a6,5(x,a13) )
 2        FORMAT( i3,x,f6.1,5(x,E13.6) )
          call exit(1)
      endif
   
C -----------------------------------------------------------------------------
C START$RUN

      spixel=pixel+1           ! for C/Fortran i/o index conversions

      ierror   = 0
      IERopt   = 0

      out_refl = 1             !for output in values of reflectance

C -------------------------------------------------------------------------
      pi    = 4.*atan(1.0)
      index = 1.334
      mu0   = cos(theta0 * pi / 180.0)

C -------------------- COMMON/SPEC_MISC/
      mu0_c = mu0
      pi_c  = pi

      DO j=1, NLAMBDA
        F0_c(j) = F0ln(j)
      ENDDO

C ------------------- season flag

c     note: ispring=1,isummer=2,ifall=3
      iseason = 2
             
c******************* PIXEL-DATA OUTPUT - INITIALISATION *****************

      pix_to_screen = 0        ! output to screen

c print-to-screen filter

      IF(pix_to_screen.eq.1) THEN  
          PRINT*,'start: Line/Pixel - ', scan, pixel+1
          PRINT*, lat_in, lon_in, theta,theta0,delt_phi
      ENDIF

      scan1=scan-998
      pixel1=pixel+1   

c****************** END PIXEL-DATA OUTPUT (INITIALISATION) ********************

         
C ---------- This should be changed when LUTs are complete for all THETA0

      IF (theta0.LT.0.0 .OR. theta0.GE.77.5 .OR.
     &    theta.GT.77.5 ) THEN
         ierror = 1
         DO j=1,NLAMBDA
           totWatRad(j)=-1.0
           t_star(j) = 1.0
           t_star0(j) = 1.0
         ENDDO
         GOTO 777
      ENDIF

C ----------- INITIALIZING ALL LOOK-UP ARRAYS ---------------------------------
      IF (FIRSTB) THEN
        PRINT*,'Loading Haze-C Model Look-Up Tables ...'
        CALL READ_MODEL_DATA(sensorNm)
        PRINT*,'Starting Atmospheric Correction'
        CALL INIT_WATER_PAR1
        FIRSTB= .FALSE.
      ENDIF
c      CALL INIT_WATER_PAR2(lat_in,lon_in,iseason)    !Chesapeake coeffs   
   
C ----------- END initializing all look-up arrays -----------------------------

c       Force delt_phi into its canonical range of [-180,+180]
        if (delt_phi.lt.-180.) then
          delt_phi = delt_phi + 360.
        else if (delt_phi.gt. 180.) then
          delt_phi = delt_phi - 360.
        end if

c ----------- apply O_2 correction to SeaWIFS band 7 (see gordon_o2.f)

      air_mass = 1./mu0 + 1./cos(rad(theta))
   
      DO j = 1,NLAMBDA
        rhot(j) = (pi*totRad(j))/(F0ln(j)*mu0)
      ENDDO
      
      DO j = 1,NLAMBDA

       IF(j.NE.i765) THEN
        aerwatRho(j) = rhot(j) - rayRho(j)
        aerwatRho_c(j) = aerwatRho(j)

        ELSE    !include O2 correction at 765
          aerwatRho(j) = (rhot(j) 
     &               - rayRho(j))*FUNCT_OXYGEN_AER(air_mass)
          aerwatRho_c(j) = aerwatRho(j)

        ENDIF

       ENDDO


C ************************************************************************  
C ******************************  SOA  ***********************************
C ************************************************************************ 

      CALL SPEC1GEOFIT(F0ln,lat_in,theta0,theta,delt_phi,
     &                 sunTAB,
     &                 water_vapor,
     &                 optMI,optMR,optPIG,optOC4,
     &                 optACDM,optB0,
     &                 optW0,optTAU,optSTDDEV,
     &                 optTAUL,optTrans,optTrans0,
     &                 resAerRho,resWatRho,resTotRho,
     &                 intpNu,IERopt,fl,absCoef)

C ************************************************************************ 
C ************************************************************************
    
C ------ error status check

      IF(IERopt.NE.0) THEN
         ierror=1
         DO j=1,NLAMBDA
           totWatRad(j)=-1.0
           t_star(j) = 1.0
           t_star0(j) = 1.0
         ENDDO
         chl_out(spixel)    = 0.0
         acdm_out(spixel)   = 0.0
         bbp_out(spixel)    = 0.0
         chl4s_out(spixel)  = 0.0
         pcentcdm_out(spixel) = 0.0 
         w0_out(spixel)     = 0.0
         v_out(spixel)      = 0.0
         GOTO 777
      ENDIF
 
c --------- prepare acdm/atot 
      pcentcdm=optACDM/absCoef(i443)
 
C --------- PREPARE OUTPUT FOR FURTHER PROCESING

      DO j=1,NLAMBDA
        taua(j)    = optTAUL(j)
        t_star(j)  = optTrans(j)
        t_star0(j) = optTrans0(j)
        abs_out(j)=  absCoef(j)
      ENDDO
c      chol     = (optPIG/1.34)**1.0173

      chl_out(spixel)    = optPIG
      acdm_out(spixel)   = optACDM
      bbp_out(spixel)    = optB0

      pcentcdm_out(spixel) = pcentcdm
     
      w0_out(spixel)     = optW0
      v_out(spixel)      = intpNu

      IF(ABS(resWatRho(5)).GT.0.000001) THEN
           chl4s_out(spixel)=optOC4
      ELSE
           chl4s_out(spixel)=0.0
      ENDIF


C -------------------- Aerosol radiances
      DO j=1,NLAMBDA
         tla(j) = resAerRho(j)*F0ln(j)*mu0/pi  !with cosine
      ENDDO

C --------- WATER-LEAVING REFLECTANCE AT SENSOR (and channel)

      DO j=1,NLAMBDA
         opt_nRhoW(j)=optTrans0(j)*resWatRho(j)*optTrans(j)
      ENDDO

C --------- WATER-LEAVING RADIANCE
C Note:
C           totWatRad(1..6) = WaterRadiance (non-aerosol bands)
C           totWatRad(7..8) = AerosolRadiance (aerosol band)

      DO j=1,NLAMBDA
crc         totWatRad(j) = resWatRho(j)*F0ln(j)*mu0/pi
         totWatRad(j) = resWatRho(j)*F0ln(j)/pi  !no need for cosine
      ENDDO


777   CONTINUE

c*********************** BEGIN PIXEL-DATA OUTPUT *******************************


c debugging output to screen

      IF(pix_to_screen.eq.1) THEN
          PRINT*,'C =',     optPIG
          PRINT*,'OC4 =',   optOC4
          PRINT*,'acdm =',  optACDM
          PRINT*,'bbp =',   optB0          
          PRINT*,'W0 =',    optW0
          PRINT*, 'trhow =', opt_nRhow(2)
      ENDIF


c********************** END PIXEL-DATA OUTPUT ********************************

      RETURN
      END

C ############################################################################
C ############################################################################

C SPEC1GEOFIT is the central contolling centre of the entire code.
C It initialises LUT's for diffuse transmittance coefficients for each fixed 
C value of (mi,mr,nu) combination at each wavelength, and for interpolated geometry 
C (theta,theta0). It then calls two main driving routines: DRV_TAUVV_ZXMWD 
C (using subroutine ZXMWD/E for NIR optimization) and DRV_SPECOPT_LBFGSB 
C (for VIS optimization). See respective routines for further details.

      SUBROUTINE SPEC1GEOFIT(F0ln,lat_in,theta0,theta,delt_phi,
     &                       sunTAB,
     &                       water_vapor,
     &                       optMI,optMR,optPIG,optOC4,
     &                       optACDM,optB0,
     &                       optW0,optTAU,optSTDDEV,
     &                       optTAUL,optTrans,optTrans0, 
     &                       resAerRho,nLw_out,resTotRho, 
     &                       optVV,IERopt,fl,absCoef)

c--input: 
c         theta0         -solar zenith (deg)
c         theta          -sensor zenith (deg)
c         delt_phi       -sensor relative azimuth (deg)
c         sunTAB         -tabulated sun angles (see DATA declaration - SUBROUTINE ATMCOR_SOA)
c--output:
c         optMI          -optimized imaginary part of refractive index
c         optMR          -optimized real part of refractive index
c         optPIG         -optimized chlorophyll concentration (mg/l)
c         optACDM        -optimized cdm absorption at 443-nm (= agelb+adet)
c         optB0          -optimized particulate backscattering coeff.at 443 nm
c         optW0          -optimized single scattering albedo
c         optTAU         -optimized aerosol optical thickness at 865-nm
c         optSTDDEV      -standard deviation of master function in main optimization (DRV_SPECOPT_LBFGSB)
c         optTAUL        -optimized aerosol optical thickness at lambda (each wavelength)
c         optTrans       -diffuse transmittance outbound (at theta)
c         optTrans0      -diffuse transmittance inbound (at theta0)
c         resAerRho      -resultant aerosol reflectance rhoA (= aerosol reflectance, prime)
c         resWatRho      -resultant water reflectance rhow (= water reflectance at surface, prime)
c         resTotRho      -resultant total reflectance rhot (= total reflectance at TOA, prime)
c         optVV          -optimized aerosol Junge distribution parameter (nu)
c         ierOPT         -logical error-flag for NIR optimization OR Newt-Raph search (for inside world)

c--other variables:
c         NMR            -#of REAL parts of REFRACTIVE INDEXES  m = mr + j*mi
c         NMI            -#of IMAG parts of REFRACTIVE INDEXES  m = mr + j*mi
c         NV             -#of Nu's
c         vstep          -Nu_grid_spacing
c         sunstep        -theta spacing in sunTAB
c         ierr           -error-flag in Newt-Raph search
c         ytmp           -3 tmp-values for QUAD-interpolation purposes
c         doLog          -logical whether to do LOG-interpolation

c Case 2 variables:
c         rhow_diff      -residual water reflectance at 865-nm
c         rhow_min       -minimum threshold for case 2 loop.
c                         *also used as minimum convergence threshold for Case 2 iterations

C START$DECLARE
      IMPLICIT NONE

      INTEGER    NMR,NMI,NV,NLAMBDA,NUM,NBAND,
     &           NMODEL,NSUN
      PARAMETER (NMI=6,NMR=2,NV=6,NSUN=33)
      PARAMETER (NMODEL=NV*NMR*NMI)
      PARAMETER (NLAMBDA=8,NBAND=6)
      PARAMETER (NUM=75)

      INTEGER    wavelength(NLAMBDA)
      INTEGER    i,ivBOT,
     &           ji,jr,jl,jv,
     &           jmodel,jMin

C ------ REAL INPUT
      REAL  theta,theta0,delt_phi,
     &      theta_c,theta0_c,delt_phi_c,
     &      sunTAB(NSUN),sunstep,
     &      aerwatRho(NLAMBDA)

      REAL    water_vapor,water_vapor_c

      REAL    aerW0(NLAMBDA,NMI,NMR,NV), aerC(NLAMBDA,NMI,NMR,NV),
     &        cc865tab(NLAMBDA,NMI,NMR,NV),
     &        s11dummy(NUM)

      REAL    diffTrans_a(NSUN,NLAMBDA,NMODEL),        ! LUT
     &        diffTrans_b(NSUN,NLAMBDA,NMODEL)
      REAL    diffA_vri(NV,NMR,NMI,NLAMBDA),           ! given ThetaView
     &        diffB_vri(NV,NMR,NMI,NLAMBDA),
     &        diffA_ri(NMR,NMI,NLAMBDA),               ! for interpolated V
     &        diffB_ri(NMR,NMI,NLAMBDA)
      REAL    diffA0_vri(NV,NMR,NMI,NLAMBDA),          ! given ThetaSun
     &        diffB0_vri(NV,NMR,NMI,NLAMBDA),
     &        diffA0_ri(NMR,NMI,NLAMBDA),              ! for interpolated V
     &        diffB0_ri(NMR,NMI,NLAMBDA)

      REAL  F0(NLAMBDA),mu0,pi
      REAL  x,phi,rad
      REAL  v(NV),mi(NMI),mr(NMR)
      REAL  w0nu(NMI,NMR)

C -- OUTPUT/OPTIMIZATION
      INTEGER  IERopt,IERR
      REAL     optSTDDEV,optMI,optMR,optPIG,optACDM,optB0,
     &         optVV,optTAU,optW0,
     &         optTAUL(NLAMBDA),optTrans(NLAMBDA),optTrans0(NLAMBDA)
      REAL     resAerRho(NLAMBDA),resWatRho(NLAMBDA),
     &         resTotRho(NLAMBDA)
      REAL     resLwN(NLAMBDA)

      REAL     OC4_CHL4S
      EXTERNAL OC4_CHL4S 
       REAL     optOC4

      REAL     fl(NLAMBDA),absCoef(NLAMBDA)

C ------ REAL FOR Nu-CALCULATIONS
      REAL  candRho(NLAMBDA,NMI,NMR), rcandRho(NMI,NMR,NLAMBDA)
      REAL  cc865(NMI,NMR,NLAMBDA)
      REAL  vv_c,vstep
      REAL  ytmp(3),tau865(NMI,NMR)

      INTEGER isDoCase2,iterCase2,NITERCASE2,totCase2Iter
      REAL    tauray(NLAMBDA),rho_wat_giv(NLAMBDA)

      INTEGER isCandData(NMI,NMR)   !flag for accepting or rejecting optimised data corresponding 
c                                    to each possible (mi,mr) combination (=jmodel)
      LOGICAL COMPCANDDATA
      EXTERNAL COMPCANDDATA

C -- COMMON
      COMMON /SPEC_GEOM/   theta0_c,theta_c,delt_phi_c
      COMMON /COMPHASE/    aerW0,aerC,cc865tab,s11dummy

      COMMON /COMDIFF/     diffTrans_a,diffTrans_b
      COMMON /COMDIFF_MISC/ diffA_vri ,diffB_vri,
     &                      diffA0_vri,diffB0_vri

      COMMON /SMM_GENERIC/ wavelength
      COMMON /OPTFUNCREAL/ mi,mr,cc865,tau865,
     &                     diffA_ri,diffB_ri,
     &                     diffA0_ri,diffB0_ri,
     &                     aerwatRho,candRho,
     &                     rcandRho
      COMMON /CAND_MISC/   isCandData,w0nu
      COMMON /SPEC_MISC/   F0,mu0,pi
      COMMON /TAUVV_MISC/  tauray,rho_wat_giv
      COMMON /CASE2_FLAG/  isDoCase2
      COMMON /V_MISC/      v,vv_c,vstep,ivBOT

      COMMON /watv/        water_vapor_c 

      REAL   rhow_diff(0:4)

      REAL rhow_min
      DATA rhow_min /0.0001/            ! Case2 minimum rhow tolerance in NIR 0.0005
      SAVE rhow_min

      COMMON /CASE2_ITER/ iterCase2,NITERCASE2

c - water pars
      REAL*4 adm_s,bbp_s,aw(NLAMBDA),bbw(NLAMBDA),aph_lin(NLAMBDA)
      COMMON /WATERPAR/ adm_s,bbp_s,aw,bbw,aph_lin

c - high accuracy aerosol and startups
      INTEGER     NN
      PARAMETER   (NN=5)

      INTEGER j
      REAL    F0ln(NLAMBDA),nLw_out(NLAMBDA)

      REAL lat,lat_in
      COMMON /latitude/lat

      rad(x)  = x*pi/180.

      phi     = rad(delt_phi)

      water_vapor_c = water_vapor

      lat=lat_in

C----
      theta0_c     = theta0
      theta_c      = theta
      delt_phi_c   = delt_phi
      IERopt       = 0           !internal error
      IERR         = 0           !     "
        
C---- controlling flags for Case2 iterations
      isDoCase2    = 1   ! 1 - assume Case-2 waters; 0 - no water-radiance at Bands-7,8
      iterCase2    = 0   ! current iteration of Case2 loop
      NITERCASE2   = 4   ! max number of iterations in Case2 loop  
      totCase2Iter = 0   ! previous (last completed) Case2 loop

C---- define step intervals in respective tables
      vstep   = 0.5      !nu step interval
      sunstep = 2.5      !geometry step interval (for theta, theta0)

      DO i=1,NITERCASE2+1
        rhow_diff(i-1) = 0.0  
      ENDDO

C ---- CALCULATE DIFFUSE TRANSMITTANCE COEFFICIENTS A(view),B(view) FOR ALL Mr,Mi,Nu (all fixed at this point) ----

c diffuse transmittance coeffs interpolated here for input theta,theta0, in preparation for
c nu-interpolation (NIR optimization) and then mi,mr interpolation (main optimization)
c Quadratic interpolation is used.

c ---------------------------DIFFTRANS-THETA (sensor angle)--------------------------------

c     select index number (of sunTAB) that is closest to (but less than) theta
c     i.e., locate minimum index number (jmin) in sunTab array that corresponds to 
c     input sensor angle theta:
      CALL BRACKET4QUAD (NSUN,sunTAB,theta,jMin)  !jmin=1::33 

      DO jv=1,NV    !1::6
      DO jr=1,NMR   !1::2
      DO ji=1,NMI   !1::6

c      initialise Haze-C model given by its parameters into a unique model index number
c      i.e. (jv,jr,ji) -> (jmodel) 
c      Each of possible 72 (mi,mr,nu) combinations are assigned a unique model number
c      Hence, we have a LUT 3D grid (6*2*6) of jmodel indices:
       CALL MODEL2IDX(jv,jr,ji,NV,NMR,NMI,jmodel)  !jmodel=1::72

       DO jl=1,NLAMBDA  !1::8

c ---------------------A-coeff
c       Initialise diffuse transmittance A-coeff for 3 consecutive jmin (theta index) steps.
c       i.e., for grid point at (theta_index,lambda,mi,mr,nu) we initialise diffTrans_a at each
c       theta index (jmin,jmin+1,jmin+2), the latter determined by real theta input.
c       jmin and jmin+1 bound the real (input) value of theta
        ytmp(1)=diffTrans_a(jMin  ,jl,jmodel)   !diffTrans_a(thetaIndex  ,lambdaIndex,modelIndex)
        ytmp(2)=diffTrans_a(jMin+1,jl,jmodel)   !     "      thetaIndex+1      "          "
        ytmp(3)=diffTrans_a(jMin+2,jl,jmodel)   !     "      thetaIndex+2      "          "

c       quadratic interpolation (for A-coeff) in 3D diffTrans_a grid based on y(x)=ax^2+bx+c
c       result is diffA_vri (in nu,mr,mi,lambda) interpolated for geometry (3 consecutive theta (index) steps)
c       methodology: jmin and jmin+1 are the sunTAB(index) boundaries of the
c       'real' input theta (see subroutine BRACKET4QUAD).
c       jmin+2 completes the coordinate trifecta for the polynomial interpolation, 
c       but its the jmin to jmin+1 curve shape that defines the interpolation.
c       once the 'real' theta point is found (for fixed  nu,mr,mi,lambda), interpolations
c       can then proceed in other dimensions (eg., nu in  DRV_TAUVV_ZXMWD)
        CALL QUADINTRP(.FALSE.,
     &                 theta,sunTAB(jMin),sunstep,
     &                 ytmp,diffA_vri(jv,jr,ji,jl)) !(allpositive,xval,start_val,2.5deg,yData,y=result)

C ---------------------B-coeff
c       initialisation of diffuse transmittance B-coeff theta grid triplet (in lambda,nu,mr,mi)
c       with theta(jmin) indices as for A-coeff above
        ytmp(1)=diffTrans_b(jMin  ,jl,jmodel)  
        ytmp(2)=diffTrans_b(jMin+1,jl,jmodel)
        ytmp(3)=diffTrans_b(jMin+2,jl,jmodel)

c       quadratic interpolation of theta for B-coeff (as above)
        CALL QUADINTRP(.FALSE.,
     &                 theta,sunTAB(jMin),sunstep,
     &                 ytmp,diffB_vri(jv,jr,ji,jl))

       ENDDO ! jl
      ENDDO  ! ji
      ENDDO  ! jr
      ENDDO  ! jv
c-------------------------------end DIFFTRANS-THETA------------------------------------------

c -----------------------------DIFFTRANS-THETA0 (SZA)----------------------------------------

c same as section above (DIFFTRANS-THETA) but for theta0 (sun angle) 

      CALL BRACKET4QUAD (NSUN,sunTAB,theta0,jMin)

      DO jv=1,NV
      DO jr=1,NMR
      DO ji=1,NMI

       CALL MODEL2IDX(jv,jr,ji,NV,NMR,NMI,jmodel)

       DO jl=1,NLAMBDA

C -----------------------A0-coeff

        ytmp(1)=diffTrans_a(jMin  ,jl,jmodel)
        ytmp(2)=diffTrans_a(jMin+1,jl,jmodel)
        ytmp(3)=diffTrans_a(jMin+2,jl,jmodel)
        CALL QUADINTRP(.FALSE.,
     &                 theta0,sunTAB(jMin),sunstep,
     &                 ytmp,diffA0_vri(jv,jr,ji,jl))

C ------------------------B0-coeff

        ytmp(1)=diffTrans_b(jMin  ,jl,jmodel)
        ytmp(2)=diffTrans_b(jMin+1,jl,jmodel)
        ytmp(3)=diffTrans_b(jMin+2,jl,jmodel)
        CALL QUADINTRP(.FALSE.,
     &                 theta0,sunTAB(jMin),sunstep,
     &                 ytmp,diffB0_vri(jv,jr,ji,jl))
       ENDDO ! jl
      ENDDO  ! ji
      ENDDO  ! jr
      ENDDO  ! jv

c-------------------------------end DIFFTRANS-THETA0----------------------------------


C ************************************************************************  
C ******************CASE2 RADIOMETRIC LOOP STARTS HERE *******************
C ************************************************************************  

      rho_wat_giv(7)=0.0   !resultant rhow at 765 (initialised)
      rho_wat_giv(8)=0.0   !         "        865       "

 777  CONTINUE

      iterCase2=iterCase2+1       ! update current iteration of Case2 loop

C -------------- INITIALIZE COMPUTED/COMPUTABLE INFO FOR ALL MODELS
      DO ji=1,NMI
      DO jr=1,NMR
        isCandData(ji,jr)=0     !flag for model (mi,mr) acceptance/rejection
      ENDDO  ! ji
      ENDDO  ! jr

C ************************************************************************  
C ***************** ESTIMATION OF FUNCTIONS NU AND TAU  ******************
C ************************************************************************ 

C COMPUTE COMPUTABLE NIR INFO FOR ALL MODELS
c NIR optimization is now performed for each possible (mi,mr) combination  = (6*2) or (6*1)
c if non-computable,progam ends and error-failure flag sent to subroutine ATMCOR_SOA
c The mi,mr 6*2 grid is 'referenced' to a nu grid in NIR optimization.  Nu and tau865 
c 'candidates' are processed and the optimal values derived by comparing the corresponding
c estimated rhoA(NIR) with rhoA(NIR) measured.  Everything is iterpolated for nu, so that 
c we have parameters that are refenced to mi,mr (exact table values) and 
c optimized nu.  This forms the basis for input into the main optimization 
c (NIR optimization: 
c for each exact mi,mr and 'optimal' nu we now have rhoA(lambda),tau865,tau(lambda)/tau(865),
c difftrans A/B coeff at theta/theta0, t0(inbound), t(outbound) and ssalbedo (w0).

      DO ji=1,NMI
      DO jr=1,NMR
        IF (.NOT.COMPCANDDATA(ji,jr)) THEN
          IERopt=1
          GO TO 999
        ENDIF
      ENDDO  ! ji
      ENDDO  ! jr

C ************************************************************************  
C ************************ SPECTRAL OPTIMIZATION *************************
C ************************************************************************ 

C OPTIMIZING STANDARD DEVIATIONS
c The main optimization derives mi,mr,C,acdm,bbp by using 'candidates' of mi,mr,C,acdm,bbp
c to first estimate MODELLED[RhoA+rhow].  All parameters in main optimization are interpolated 
c for the current mi,mr 'candidates' being processed for selection.
c LSQ Comparison with MEASURED[rhot-rhoray] quantifies the performance of the optimization and the
c final candidates are then used to compute the final values of tau(lambda),t,t0,rhoA,rhow
c and rhoA+rhow.

      CALL MAIN_OPTIM_INOUT_PROC(IERR,
     &                       optMI,optMR,optPIG,optACDM,optB0,
     &                       optVV,optTAU,optW0,optSTDDEV,
     &                       optTAUL,optTrans,optTrans0,
     &                       resAerRho,resWatRho,resTotRho,
     &                       fl,absCoef)   !(re)optimize bands 1 to 6

C ************************************************************************  
C ************************************************************************

c     IERopt=0=success of main optimization computation; 1 = failure
      IF(IERR.EQ.136) THEN
         IERopt=1
         GOTO 999
      ENDIF

      IF(isDoCase2.EQ.1.AND.iterCase2.LE.NITERCASE2) THEN
         CALL GS97_CASE2(optPIG,optACDM,optB0,resLwN)    !calculate water leaving radiance
         DO jl=7,8                  !calculate NIR water reflectance, both channels 
           resWatRho(jl)   = pi*resLwN(jl)/(F0(jl)*mu0) 
           rho_wat_giv(jl) = resWatRho(jl)
         ENDDO
     
         totCase2Iter = iterCase2-1                      ! number of Big Loop Iterations
                                                         ! in addition to first (initial) loop

c        ends case 2 loop if water reflectance at 865-nm is less than minimum threshold
         IF (resWatRho(8).LE.rhow_min) GOTO 800      ! this is Case1 water

c        take difference between re-computed rhow and previous residual 
         rhow_diff(iterCase2)=resWatRho(8)-rhow_diff(iterCase2-1)

c        if difference is less than theshold, iteration has converged.  Else, re-compute rhow
         IF (ABS(rhow_diff(iterCase2)-rhow_diff(iterCase2-1))
     &           .GT.rhow_min) THEN                                !check against threshold
            GOTO 777            ! assume Water Reflectance at Bands-7,8 and repeat the loop
         ELSE
            ! finish
         ENDIF

      ENDIF

C *************************** END CASE2 LOOP ******************************

800   CONTINUE

c store final nLw: at 'selected' accuracy
      DO j=1,NLAMBDA
        nLw_out(j)= resWatRho(j)  !lower accuracy nLw (when high_prec = 0)
      ENDDO 
      optOC4 = OC4_CHL4S(nLw_out)


999   RETURN
      END

C ############################################################################
C ############################################################################

C SECTION 1:

C This block of code uses spectral optimization (NOT epsilon) to derive nu
C and tau(865) for each combination of (mr,mi) at the NIR wavelengths.
   
      LOGICAL FUNCTION COMPCANDDATA (ji,jr)    !compute candidate date (NIR)
c main calling function for NIR optimization (DRV_TAUVV_ZXMWD); called by SOA main 
c controlling routine SPEC1GEOFIT.
c Nested in (mi,mr) loop (6*2) or (6*1)
c Purpose is to initiate NIR optimization of nu and tau [for each possible (mr,mi) 
c combination] and for all models.
c boundaries of nu are also defined for post-processing purposes.
 
c--input: 
c         ji             -loop index for NMI, #of IMAG parts of REFRACTIVE INDEXES
c         jr             -loop index for NMR, #of IMAG parts of REFRACTIVE INDEXES 
c                        *(ji,jr) = (6,2) = 12 combinations of (mi,mr)  
c--output: 
c         no explicit output; logical function = implicit output
c--other vars:
c         SlsqVTau       -standard deviation of master function in NIR optimization (DRV_TAUVV_ZXMWD)
c         ierVTau        -logical error-flag for error status (for inside world)
c                         *error status from NIR optimization (DRV_TAUVV_ZXMWD)

      INTEGER ji,jr    !loop index number for mi,mr

      INTEGER    NMR,NMI
      PARAMETER (NMI=6,NMR=2)

      INTEGER isCandData(NMI,NMR)
      REAL    vir(NMI,NMR)      ! nu derived from NIR 
      REAL    w0nu(NMI,NMR)     ! SingScatAlb corresponding to index (iMR,iMI) (given Nu) at 865nm
      COMMON /CAND_MISC/   isCandData,w0nu
      COMMON /V_EXTRA/     vir

      INTEGER ierVTau(NMI,NMR)
      REAL    SlsqErr(NMI,NMR)

      CALL DRV_TAUVV_ZXMWD 
     &         (ji,jr,w0nu,vir,SlsqErr,ierVTau)

      isCandData(ji,jr) = 1
      COMPCANDDATA      = .TRUE.

      IF(ierVTau(ji,jr).EQ.136 .OR.
     &   vir(ji,jr).LT.1.5     .OR.
     &   vir(ji,jr).GT.6.0         ) THEN        !boundary check for nu; 
cc     &   SlsqErr(ji,jr)   .GT.10.0    ) THEN
         isCandData(ji,jr) = -1                  !candidate model (mi,mr) is rejected
         COMPCANDDATA = .FALSE.
      ENDIF

      RETURN
      END

C ****************************************************************************

      SUBROUTINE DRV_TAUVV_ZXMWD(jmi,jmr,w0nu,vir,    !Drive tau/nu using ZXMWD
     &                           SlsqVTau,ierVTau )

c Main NIR driving function.
c Used in replacement of SUBROUTINE DRV_SPECOPT_ZXMWD in last Chomko code Dec 2000
c Nested in (mi,mr) loop (6*2)

c Constraints of the relationship nu(mi,mr) and tau865(mi,mr) are given by the
c HazeC LUT's (ssalbedo,extinction coeffs,rhoA~tau fourier coeffs, diff trans coeffs)

c This subroutine performs two distinctly different operations:

c 1)  Computes tau(865) and nu (at each mi,mr) using NIR values of rhoA [rhoA(765::865)].
c For this purpose a call to subroutine FCN_TAUVV_DIFFTRAN is made in the NIR optimization.
c This is the driving function for Hessian estimation.  It compares 'candidate' rhoA (derived from 
c constraints of the nu(mi,mr) and tau865(mi,mr) HazeC relationships) with real rhoA 
c [rho_total-rho_rayleigh-rho_water] to derive SD 'updates'.

c 2) Given updated/revised 'optimal' nu and tau865 (optimized by ZXMWD/E call to FCN_TAUVV_DIFFTRAN),
c performs (updated/revised) interpolation of optimal parameters (for each mi,mr) for ALL wavelengths
c and optimal nu. 

c When tau(865) and nu are not computable (from (mr,mi)) LUT interpolation 
c (of mi as a function of mr) is done.
c We refer to nu in this subroutine as 'optimal' nu (explained later)

c--input: 
c         jmi(=jji)       -loop index for NMI, #of IMAG parts of REFRACTIVE INDEXES
c         jmr(=jjr)       -loop index for NMR, #of IMAG parts of REFRACTIVE INDEXES 
c                         *(ji,jr) = (6,2) = 12 combinations of (mi,mr)
c--output:
c         w0nu           -ssalb corresponding to (iMR,iMI) (given Nu) at 865nm
c         vir            -optimal nu for given (mr,mi) combination
c         SlsqVTau       -standard deviation of master function in NIR optimization (DRV_TAUVV_ZXMWD)
c         ierVTau        -logical error-flag for error status (for inside world)
c                         *error status from NIR optimization (DRV_TAUVV_ZXMWD)
c--other variables: 
c         jmodel         -mr,mi,nu set up grid (at respective table points)
c         ivBOT          -corresponds to the location of nu in nu-mini-table
c         ivMDL          -                    // 
c         ivTOP          -                    // 
c         cc865          -Cext(lambda)/Cext(865)=tau(lambda)/tau(865) for nu
c         cc865tab       -the above for tabulated Nu-data at 865nm
c         tau865         -optimal NIR(865-nm) value of AOD for given (mr,mi) combination

      IMPLICIT NONE

      INTEGER   NN
      PARAMETER(NN=2)     !number of parameters to be optimized (nu,tau865)

      INTEGER  nsig,nsrch,IOPT,NRST,iwork(NN),IERopt  !tolerance vars ZXMWD
      REAL     x(NN),f                                !output: x(1::2)=nu,tau865; f=SD master function
      SAVE     x
      REAL     lowB(NN),upB(NN)                       !constraint (startup): lowB(1),upB(1)=nu; lowB(2),upB(2)=tau865
      REAL     work((NN+12)*NN)

      LOGICAL  FCN_TAUVV_DIFFTRAN      !surrogate for FCN in ZXMWD, FUNCT in ZXMWE
      EXTERNAL FCN_TAUVV_DIFFTRAN

      INTEGER   NV,NMI,NMR,NLAMBDA
      PARAMETER (NMI=6,NMR=2,NV=6,NLAMBDA=8)
      INTEGER   NUM
      PARAMETER (NUM=75)

      INTEGER  jmi,jmr,ierVTau(NMI,NMR)
      REAL     theta0,theta,delt_phi
      REAL     mi(NMI),mr(NMR),v(NV)
      REAL     vv_c,vstep,vir(NMI,NMR),SlsqVTau(NMI,NMR), !SlsqVTau = SD[rhoA(real@765,865)~rhoA(retieved@765,865)]
     &         cc865tab(NLAMBDA,NMI,NMR,NV)                           !tabulated tau(lambda)/tau(865)
      REAL     aerW0(NLAMBDA,NMI,NMR,NV), aerC(NLAMBDA,NMI,NMR,NV),   !aerW0 = tabulated ssalbedo (mr,mi)@lambda,mi,mr,nu
     &         s11dummy(NUM)

      REAL     aerwatRho(NLAMBDA),aerwatRho_c(NLAMBDA)  !aerwatrho = aerwatRho_c = rhot-rhor
      REAL     w0nu(NMI,NMR)                                          !nu-interpolated ssalbedo at 865-nm (mr,mi)
      REAL     candRho(NLAMBDA,NMI,NMR),rcandRho(NMI,NMR,NLAMBDA),    !candidate RhoA for (visible LBFGSB optimisation later)
     &         cc865(NMI,NMR,NLAMBDA),tau865(NMI,NMR)  !nu-intepolated tau(lambda)/tau(865);optimal tau(865)

      REAL     A(NLAMBDA),B(NLAMBDA),C(NLAMBDA),D(NLAMBDA)
      REAL     diffA_vri(NV,NMR,NMI,NLAMBDA),   ! diffuse transmittance A coef pre-interpolated for Theta
     &         diffB_vri(NV,NMR,NMI,NLAMBDA),   !            "          B                "
     &         diffA_ri(NMR,NMI,NLAMBDA),       ! diffuse transmittance A coef interpolated for Theta and nu  
     &         diffB_ri(NMR,NMI,NLAMBDA),       !            "          B                "       
     &         diffA0_vri(NV,NMR,NMI,NLAMBDA),  ! diffuse transmittance A coef pre-interpolated for Theta0
     &         diffB0_vri(NV,NMR,NMI,NLAMBDA),  !            "          B                "
     &         diffA0_ri(NMR,NMI,NLAMBDA),      ! diffuse transmittance A coef interpolated for Theta0 and nu          
     &         diffB0_ri(NMR,NMI,NLAMBDA)       !            "          B                "  

      REAL     tauA       !temp var for aer optical depth spectral interpolation [from tau(865)]

      INTEGER  jmodel,jvv,
     &         jji,jjr,ivBOT_c,ilow,jl,jv,jt

      INTEGER  ivBOT,ivMDL,ivTOP !nu mini-table indices. It beceomes superimposed 
c                                on v table (2.0,2.5,3.0,3.5,4.0,4.5) for nu interolation of all parameters
                                   
      REAL     ytmp(3),tauGrid(5),atmp(0:5),
     &         RoAtab(5,NLAMBDA,3),RoAvv(5),   !RhoA interpolated for geometry, interpolated for geometry and nu
     &         Acoef(NLAMBDA),Bcoef(NLAMBDA),  !a,b,c,d HazeC coeffs interpolated for geometry and nu
     &         Ccoef(NLAMBDA),Dcoef(NLAMBDA),
     &         vAcoef(NLAMBDA,NV),vBcoef(NLAMBDA,NV),  !a,b,c,d HazeC coeffs interpolated for geometry
     &         vCcoef(NLAMBDA,NV),vDcoef(NLAMBDA,NV)

      DATA     tauGrid / 0.05, 0.15, 0.30, 0.60, 0.80 /       !for tau lagrange polynomial fit

      LOGICAL  doLog

      COMMON /SPEC_GEOM/    theta0,theta,delt_phi
      COMMON /COMPHASE/     aerW0,aerC,cc865tab,s11dummy
      COMMON /OPTFUNCREAL/  mi,mr,cc865,tau865,
     &                      diffA_ri,diffB_ri,
     &                      diffA0_ri,diffB0_ri,
     &                      aerwatRho,candRho,
     &                      rcandRho
      COMMON /COMDIFF_MISC/ diffA_vri ,diffB_vri,
     &                      diffA0_vri,diffB0_vri
      COMMON /TAUVV/        jjr,jji,aerwatRho_c
      COMMON /V_MISC/       v,vv_c,vstep,ivBOT_c
      COMMON /V_ABCD/       vAcoef,vBcoef,vCcoef,vDcoef

      COMMON /NIR_OPTIM_CONSTRAINT/ lowB,upB

      REAL                  fi6(NMI,NMR),fi7(NMI,NMR),fi8(NMI,NMR)
      COMMON /F_INBAND/     fi6,fi7,fi8

      doLog = .TRUE.   !for log interpolation in quadratic interpolation

c this whole subroutine is inside an (mi,mr) loop, mi first, then mr (6*mr1 then 6*mr2)
c see COMPUTE COMPUTABLE NIR INFO FOR ALL MODELS in subroutine spec1geofit
c jji and jjr (= mi,mr loop) are used here for uniqueness between subroutines.
c nu and tau865 are calculated for each (mi,mr) combination.
      jji   = jmi
      jjr   = jmr

c vary these parameters------
c nsig=3;nsrch=3;NRST=0;IOPT=4 optimimum, Kuchinke May 2004

c C1.<OPTIMAL PARAMETERS>
C --------------------- GOOD QUALITY, MEDIUM SPEED   ( 12.8 pix/sec)
Crc      nsig  = 2
Crc      nsrch = 2
Crc      NRST  = 1
Crc      IOPT  = 4
Crc      IF(isNextPix) THEN
Crc         IF(ierOPT.EQ.0) IOPT=7
Crc      ENDIF
C --------------------- END: GOOD QUALITY, MEDIUM SPEED

C2.   <TOP QUALITY>
C --------------------- TOP QUALITY, HORRIBLE SPEED   ( 1.2 pix/sec)
Crc      nsig  = 2
Crc      nsrch = 10
Crc      NRST  = 0
Crc      IOPT  = 0
C --------------------- END: TOP QUALITY, HORRIBLE SPEED

C3.   <TOP SPEED>
C --------------------- TOLERABLE (few NOISY px) QUALITY,
C                                           TOP SPEED (13.6 pix/sec)
Crc      nsig  = 2
Crc      nsrch = 1
Crc      NRST  = 0
Crc      IOPT  = 4
C --------------------- END: TOLERABLE QUALITY, TOP SPEED

C4.   <GOOD QUALITY>
C --------------------- GOOD QUALITY, TOLERABLE SPEED (11.8 pix/sec)
Crc      nsig  = 2
Crc      nsrch = 3
Crc      NRST  = 2
Crc      IOPT  = 4
C --------------------- END: GOOD QUALITY, TOLERABLE SPEED

C5.     <NO USE>
C --------------------- HORRIBLE QUALITY (lots of NOISE, ugly contrast),
C                                          GOOD SPEED ( 13.3 pix/sec)
Crc      nsig  = 4
Crc      nsrch = 3
Crc      NRST  = 0
Crc      IOPT  = 4

c optimum settings
      nsig  = 3   ! # of digits of accuracy in otimization estimates (convergence criterion)
      nsrch = 3   ! # of starting points to be generated
      NRST  = 0   ! integer # of restarts
      IOPT  = 4   ! mode of Hessian matrix set up


C---------------------------------------------------------------

      lowB(1) = 1.5    !initial constraint boundaries of nu (for optimization input)
      upB(1)  = 5.0    !                    ''
      lowB(2) = 0.0    !constraint boundaries of tau865 (for optimization input)
      upB(2)  = 1.2    !                     ''

C ------------------------------- CREATE GRID (mi,mr,nu) ----------------------------------

c this loop creates the jmodel grid (mr,mi) -> (mr,mi,nu) and interpolates tau fourier coefficients
c from RhoA (original data obtained from the HazeC LUT's)
      DO jv=1,NV
        CALL MODEL2IDX(jv,jmr,jmi,NV,NMR,NMI,jmodel)   ! LUT for nu,mr,mi from each model (72 values)

C ----------------- SET UP rhoA v Tau fourier coefficients a,b,c,d -------------------------
C        parameter A,B,C,D; input grid: mi,mr,nu_table; interpolation for theta

c LUT values of the rhoA v Tau fourier coefficients a,b,c,d 
c have been computed as a quartic fit for each of 72 models (jmodel = mi,mr,v) 
c at each sun angle, lambda and phi (33*8*15*72 = 285 120 values).
c Values are taken from tabulated ones pre-computed as external LUT's
c and read in by subroutine READ_MODEL_DATA.
c rhoa (params) = a(params)tau + b(params)tau^2 + c(params)tau^3 + d(params)tau^4
c This calling function interpolates a,b,c,d fourier coefficients for geometry.
c Final values are stored in array coefficients vAcoef(1,jv)/B/C/D of dimension 
c (1,nu) where nu = 1 to 6 index value of real values

        CALL funct_a_b_c_d_fourier_interpolate (jmodel,
     &          theta0,theta,delt_phi,
     &          A,B,C,D )

        DO jl=1,NLAMBDA
            vAcoef(jl,jv) = A(jl)
            vBcoef(jl,jv) = B(jl)
            vCcoef(jl,jv) = C(jl)
            vDcoef(jl,jv) = D(jl)
        ENDDO

      ENDDO

C -------------- SOLVING A SYSTEM OF BOUNDED NON-LINEAR EQUATIONS IN v & tau865

c optimization for NIR nu,tau865 called here.  Implicit call to logical function FCN_TAUVV_DIFFTRAN (next subroutine).
 
c DRV_TAUVV_DIFFTRAN: a,b start-up constraints, N 
c lowB and upB are the constraint (startup) vectors for nu and tau (length = NN = 2)
c # of unknown parameters = NN = 2 = nu,tau865.

c Process:
c DRV_TAUVV_DIFFTRAN ---lowB,upB,N---> ZXMWD ---> ZXMWE  ---N,optimal 'candidates' x(1)=nu,x(2)=tau865---> FCN_TAUVV_DIFFTRAN
c                  <--best x(1),x(2),f--       (Hess.core)<-------f(update) feedback loop to ZXMWE--------
                     
c x (length NN) is an output array of the optimal 'candidate' values for nu,tau865
c f is the standard deviation of the master function, a comparison of rhoA(real)~rhoA(retrieved) in FCN_TAUVV_DIFFTRAN
c FCN_TAUVV_DIFFTRAN's initial call in ZXMWE provides the f startup value for Hessian estimation.
c final values of x(1)=nu,x(2)=tau865 and f are reached when ZXMWD (after repeated calls to FCN_TAUVV_DIFFTRAN)
c has finished (ZXMWE has found the minimus)

c IERopt is a logical error-flag for error status (for inside world)
c other parameters described in comments of subroutine ZXMWD

c FCN_TAUVV_DIFFTRAN:- input: N=num params,x(N)=optimal values nu,tau865;  
c                      output: f=SD[rhoA(real)~rhoA(mod)].  

c Hence for each (mi,mr) ZXMWE passes (start up) CANDIDATES nu,tau865 through FCN_TAUVV_DIFFTRAN to 
c retrieve f(=SD). This provides the initialisation (startup) for further processing in Hessian 
c estimations in ZXMWE.  Final f output from Hessian passed out through ZXMWD

      CALL ZXMWD(FCN_TAUVV_DIFFTRAN,IOPT,NN,nsig,NRST,lowB,upB,nsrch,
     &           x,f,work,iwork,IERopt)


C ---------------------- RECORDING OPTIMAL VALUES FOR NU,Tau865----------------------------------------

c From now on we refer to nu as 'optimal' nu.  Logic: it is the final optimised value originating from 
c nu startup values a(1),b(1) and processed by core Heesian routine ZXMWE until minimus is reached (values are optimal).

      vir(jmi,jmr)      = x(1)    !final output for nu (from ZXMWD)
      tau865(jmi,jmr)   = x(2)    !final output for tau865 (from ZXMWD)
      SlsqVTau(jmi,jmr) = f       !final SD output of core Hessian routine ZXMWE [rhoA(real)~rhoA(retieved)]
      ierVTau(jmi,jmr)  = IERopt  !internal error-check output


C -------------------SET UP NU MINI-TABLE FOR INTERPOLATION PURPOSES-----------------------------
C                     set up grid:ivbot,ivmin,ivtop from optimal nu

C     Locating the position of vir in the nu-table
c     select index number (of nu=v) that is closest to (but less than) vir
c     i.e., locate minimum index number (ilow) in v array that corresponds to 
c     optimal nu value vir: v[i]=1::6 ~ v = (2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
      CALL BRACKET4QUAD (NV,v,vir(jmi,jmr),ilow)

c     set up mini-table of nu indices.  Given optimised nu, ivBOT and ivMDL are index 
c     values of nu (v=1::6) that bound the optimal nu.  ivTOP completes the parabola.
c     note: from BRACKET4QUAD algorithm:
c     IF inputnu >= 4.5 THEN ilow=5 (=table nu value of 4.0, NOT 4.5)
c     IF inputnu <= 2.0 THEN ilow=1 (=table nu value of 2.0)
c     nu is out of interpolation range in both instances-could be dangerous (even for quadratic)
c     as it becomes an extrapolation!!
      IF(ilow+2.LE.NV) THEN
         ivBOT = ilow      !eg., xv=optnu=3.7,ilow=4,ivbot=4,ivmid=5,ivtop=6 (v=3.5,4,4.5)
         ivMDL = ivBOT+1   !Hence first two points (ivbot,ivmid) bound optimised nu
         ivTOP = ivBOT+2
      ELSE
         ivBOT = ilow - 1   !eg., xv=optnu=4.4,ilow=5,ivbot=4,ivmid=5,ivtop=6 (v=3.5,4,4.5)
         ivMDL = ivBOT+1    !Hence last two points (ivmid,ivtop) bound optimised nu
         ivTOP = ivBOT+2    ! if xv=4.7,ilow=4.5(from BRACKET4QUAD),ivbot=5,ivmid=6,ivtop=7
      ENDIF                
                                                                  
C -----------COMPUTING tau(lambda)/tau(865) FOR THE FOUND (VV=nu); mi,mr are given------------
C parameter tau(lambda)/tau(865); input grid: mi,mr,lambda,nu-mini-table value; interpolation for optimal nu

c     cc865tab = tabulated Cext(lambda)/Cext(865)=tau(lambda)/tau(865) for nu; 
c     cc865 = nu-interpolated Cext(lambda)/Cext(865)=tau(lambda)/tau(865) for nu; 
c     Cext = extinction coeff

      DO jl=1,NLAMBDA-1  !interpolate wavelengths 1::7

c       Initialise tabulated cc865tab [=tau(lambda)/tau(865)] for 3 consecutive nu 
c       (resultant mini-table) index steps.
c       i.e., for grid point at [lambda,mi,mr,nu(mini-table)] we initialise cc865tab
c       at each  mini-table index (ivbot,ivmid,ivtop), with the indices determined by
c       the optimal nu input. That is, optimal nu is bounded by nu mini-table values as previous  
        ytmp(1)=cc865tab(jl,jmi,jmr,ivBOT)
        ytmp(2)=cc865tab(jl,jmi,jmr,ivMDL)
        ytmp(3)=cc865tab(jl,jmi,jmr,ivTOP)

c       quadratic interpolation in cc865tab nu-mini-table based on y(x)=ax^2+bx+c
c       result is cc865 (in mi,mr,lambda) interpolated for nu (from nu mini-table)
c       Hence, cc865=tau(lambda)/tau(865) at (lambda,mi,mr) interpolated for optimal nu
        CALL QUADINTRP(doLog,
     &                 vir(jmi,jmr),v(ivBOT),vstep,
     &                 ytmp,cc865(jmi,jmr,jl))

      ENDDO

c     set tau(865)/tau(865) to unity of course
      cc865(jmi,jmr,NLAMBDA) = 1.0

C --- GIVEN (MR,MI) FINDING CANDIDATE REFLECTANCE BASED ON A,B,C,D COEFF'S FOR (NLAMBDA) & COMPUTED VV-----
C               parameter rhoA; set up at grid: mi,mr,lambda,tau_table,nu-mini-table value

c finding RhoA(lambda) at each tau-table-value, each mi,mr (implicit) and each nu-mini-table value (ivbot,ivmid,ivtop)
c Rhoa taken from quartic fit fourier coefficients pre-interpolated for geometry
c see call of funct_a_b_c_d_fourier_interpolate in DRV_TAUVV_ZXMWD
      DO jv=1,3

        jvv = ivBOT+jv-1      !ivbot,ivmid,ivtop

        DO jl=1,NLAMBDA       !lambda=412::865
          DO jt=1,5            !tau-table

c           Initialise RoAtab = Rho_A values at [tau-table,412::865,nu-mini-table] 
c           using respective fourier coeffs in larger [lambda,NV] set.   The nu-mini-table 
c           is determined by the optimal nu input. Optimal nu is bounded by nu-mini-table
c           values as given earlier.             
            RoAtab(jt,jl,jv) = vAcoef(jl,jvv)*tauGrid(jt)    +
     &                         vBcoef(jl,jvv)*tauGrid(jt)**2 +
     &                         vCcoef(jl,jvv)*tauGrid(jt)**3 +
     &                         vDcoef(jl,jvv)*tauGrid(jt)**4

          ENDDO ! jt
        ENDDO ! jl

      ENDDO  ! jv

C --- CALCULATING A,B,C,D COEFFICIENTS [quartic fit RhoA~tau(params)] FOR ALL LAMBDA ------------------
C parameter A,B,C,D; set up at grid: interpolated theta,mi,mr,lambda,interpolated_nu
C 1) using RhoA at grid point [mi,mr,tau-table,lambda,nu-mini-table value], interpolate RhoA for optimal nu 
C 2) given RhoA at grid point [mi,mr,tau-table,lambda,optimal nu], find polynomial fit (order 4) 
C    as a function of tau-table values.  This gives new A,B,C,D coeffs.

      DO jl=1,NLAMBDA       !lambda=412::865
         DO jt=1,5          !tau-table

           DO jv=1,3        !ivbot,ivmid,ivtop

c            Prepare tabulated RoAtab for quadratic interpolation in optimal nu
             ytmp(jv)=RoAtab(jt,jl,jv)
           ENDDO

c          quadratic interpolation in RoAtab nu-mini-table based on y(x)=ax^2+bx+c
c          Result is RoAvv (in tau-table), pre interpolated for geometry
c          (funct_a_b_c_d_fourier_interpolate) and for nu right now! 
c          (using nu-mini-table steps as input). We are still at the grid point (lambda,mi,mr) of course
           CALL QUADINTRP(doLog,vir(jmi,jmr),v(ivBOT),vstep,
     &                  ytmp,RoAvv(jt))
         ENDDO

c        tau-table contains 5 points.  Hence we fit Lagrange Polynomial of degree 4,
c        with x-axis=tau_grid; y-axis = corresponding RoAvv (at each tau-table value) and for
c        interpolated geometry and nu. 
c        We have effectively found a polynomial fit to
c        rhoa (params) = a(params)tau + b(params)tau^2 + c(params)tau^3 + d(params)tau^4      
c        The resultant equation is of the same order as the fourier quartic function
c        for LUT a,b,c,d coeffs in HazeC data.. 
c        Hence the Lagrange Polynomial coefficients obtained from this fit are our final 
c        interpolated coefficient values.   
         CALL LAGRINTRP(4,tauGrid,RoAvv,atmp)
         Acoef(jl)=atmp(1)
         Bcoef(jl)=atmp(2)
         Ccoef(jl)=atmp(3)
         Dcoef(jl)=atmp(4)
      ENDDO

C -------------------------- CALCULATING CANDIDATE RhoA(LAMBDA)--------------------------------
C  parameter rhoA; final grid: mi,mr,lambda,interpolated geometry(implied by A,B,C,D),
C                              interpolated optimal nu [implied by tau(lambda)/tau(865)]

      DO jl=1,NLAMBDA

c        tau(lambda) = optimal tau(865) * [tau(lambda)/tau(865)] for each (mi,mr)
c        tauA overwrites itself here (not an array)
         tauA = tau865(jmi,jmr)*cc865(jmi,jmr,jl)   !remember, cc865 has been interpolated for nu

c resultant (candidate) RhoA (at each lambda,mi,mr) pre-interpolated for geometry (in SPEC1GEOFIT) 
c and interpolated now for optimal nu
c Note: by candidate we mean candidate in main optimisation given selected (mi,mr)
         candRho(jl,jmi,jmr) = Acoef(jl) * tauA    +
     &                         Bcoef(jl) * tauA**2 +
     &                         Ccoef(jl) * tauA**3 +
     &                         Dcoef(jl) * tauA**4

c 'real' candidate RhoA.  For logic only.  Value of rhoA that is imported to visible 
c optimisation DRV_SPECOPT_LBFGSB
         rcandRho(jmi,jmr,jl)= candRho(jl,jmi,jmr)
      ENDDO

c convert SeaWiFS rhoA(670) from center band to band averaged re: Gordon (1995)
c analogy: we don't know rhoA(670) measured, so we can't convert this from band averaged to band centered.
c instead, we convert rhoA(670) modelled from band centred to band averaged.
         rcandRho(jmi,jmr,6) = candRho(6,jmi,jmr)*fi6(jmi,jmr)

C ---------------------------------------------------------------------------------------------------


C -------------- CALCULATE DIFFUSE TRANSMITTANCE COEFFICIENTS A(nL,Mi,Mr),B(...) --------------------
C                    FOR ALL NLambda,Mr,Mi (fixed) AND OPTIMAL (RETRIEVED) NU
c  A,B (for theta) and A0,B0 (for theta0) for 765-nm,Mr,Mi (fixed) and now interpolated for candidate
c  nu(mi,mr).  Mi,Mr and nu are now given of course: this subroutine is itself called in a (mi,mr) loop
c (see DRV_TAUVV_ZXMWD).

C ------------------------- A-coeff (for theta)
      DO jl=1,NLAMBDA

c       Initialise (mi,mr,lambda tabulated) diffA_vri for 3 consecutive nu 
c       (resultant mini-table) index steps.
c       i.e., for grid point at [nu(mini-table),mr,mi,412::865] we initialise diffA_vri
c       at each  nu mini-table index (ivbot,ivmid,ivtop), with the indices determined by
c       the optimal nu input. That is, optimal nu is bounded by nu mini-table values as previous
        ytmp(1)=diffA_vri(ivBOT,jmr,jmi,jl)
        ytmp(2)=diffA_vri(ivMDL,jmr,jmi,jl)
        ytmp(3)=diffA_vri(ivTOP,jmr,jmi,jl)

c       quadratic interpolation in diffA_vri nu mini-table based on y(x)=ax^2+bx+c
c       result is diffA_ri(in mr,mi,412::865) interpolated for optimal nu (from nu mini-table)
c       result was interpolated for geometry in SPEC1GEOFIT
        CALL QUADINTRP(.FALSE.,
     &                 vir(jmi,jmr),v(ivBOT),vstep,
     &                 ytmp,diffA_ri(jmr,jmi,jl))

C ----------------------------------------------------- B-coeff
c same as above but for B-coeff

        ytmp(1)=diffB_vri(ivBOT,jmr,jmi,jl)
        ytmp(2)=diffB_vri(ivMDL,jmr,jmi,jl)
        ytmp(3)=diffB_vri(ivTOP,jmr,jmi,jl)
        CALL QUADINTRP(.FALSE.,
     &                 vir(jmi,jmr),v(ivBOT),vstep,
     &                 ytmp,diffB_ri(jmr,jmi,jl))
      ENDDO ! jl

C --------------------------A0-coeff (for theta0)
c same as above (A and B) but for theta0
      DO jl=1,NLAMBDA

        ytmp(1)=diffA0_vri(ivBOT,jmr,jmi,jl)
        ytmp(2)=diffA0_vri(ivMDL,jmr,jmi,jl)
        ytmp(3)=diffA0_vri(ivTOP,jmr,jmi,jl)
        CALL QUADINTRP(.FALSE.,
     &                 vir(jmi,jmr),v(ivBOT),vstep,
     &                 ytmp,diffA0_ri(jmr,jmi,jl))

C --------------------------B0-coeff (for theta0)

        ytmp(1)=diffB0_vri(ivBOT,jmr,jmi,jl)
        ytmp(2)=diffB0_vri(ivMDL,jmr,jmi,jl)
        ytmp(3)=diffB0_vri(ivTOP,jmr,jmi,jl)
        CALL QUADINTRP(.FALSE.,
     &                 vir(jmi,jmr),v(ivBOT),vstep,
     &                 ytmp,diffB0_ri(jmr,jmi,jl))
      ENDDO ! jl

C -------------- DETERMINING W0(jmi,jmr) at 865-nm FOR JUST FOUND VV
c     Initialise (mi,mr,lambda tabulated) aerWO for 3 consecutive nu 
c     (resultant nu mini-table) index steps.
c     i.e., for grid point at [412::865, mr,mi,nu(mini-table)] we initialise aerW0
c     at each nu mini-table index (ivbot,ivmid,ivtop), with the indices determined by
c     the optimal nu input. That is, optimal nu is bounded by nu mini-table values as previous
      ytmp(1)=aerW0(NLAMBDA,jmi,jmr,ivBOT)
      ytmp(2)=aerW0(NLAMBDA,jmi,jmr,ivMDL)
      ytmp(3)=aerW0(NLAMBDA,jmi,jmr,ivTOP)

c     quadratic interpolation in aerW0 nu mini-table based on y(x)=ax^2+bx+c
c     result is w0nu(in mr,mi) and at 865-nm, interpolated for optimal nu (from nu mini-table)
      CALL QUADINTRP(.FALSE.,
     &                vir(jmi,jmr),v(ivBOT),vstep,
     &                ytmp,w0nu(jmi,jmr))

C END$RUN
      RETURN
      END

C ****************************************************************************
   
      LOGICAL FUNCTION FCN_TAUVV_DIFFTRAN (N,x,f)  !function tau/nu/diffuseTrans

c This function defines the constraints of the relationship between rhoA(765::865) and 
c nu,tau865,geometry for each mi,mr

c Nested in a big (mi,mr) loop = (6*2).
c Inputs optimised 'candidates' nu and tau865 -----> Derives RhoA(7,8).

c This is the driving function for Hessian estimation.  It compares 'candidate' rhoA (derived from 
c constraints of the nu(mi,mr) and tau865(mi,mr) HazeC relationships) with real rhoA 
c [rho_total-rho_rayleigh-rho_water] to derive SD 'updates'.
c The values of nu and tau865 are optimized using core optimization routines 
c ZXMWD and ZXMWE(nucleus for Hessian).
c As nu,tau865 are revised according to minimus in optimization, FCN_TAUVV_DIFFTRAN
c performs re-interpolation of optimal parameters (for each mi,mr) for NIR wavelengths
c and optimal nu.

c Detail:
c Rhoa(7,8) is calculated from nu-interpolated a,b,c,d coefficients (pre-interpolated for geometry
c in funct_a_b_c_d_fourier_interpolate).
c Derives diffuse transmittance coefficients (7,8) also interpolated for nu (pre-interpolated for
c geometry in SPEC1GEOFIT).
c Any rhow(7,8) (see Case 2 loop in SPEC1GEOFIT) is subtracted from rhoTotal(7,8)-rhoRayleigh(7,8) to give 
c TRUE RhoA(7,8).  This is compared with retrieved Rhoa(7,8) to give SD for 'candidate' nu

c We refer to optimised nu in this routine as 'candidate' nu.  Logic: it is updated (iteratively)
c from nu startup values a(1),b(1) by ZXMWE (including iterative calls by ZXMWE to this routine) until 
c minimus is reached (values are optimal).

c Process:
c DRV_TAUVV_DIFFTRAN ---lowB,upB,N---> ZXMWD ---> ZXMWE  ---N,optimal 'candidates' x(1)=nu,x(2)=tau865---> FCN_TAUVV_DIFFTRAN
c                  <--best x(1),x(2),f--       (Hess.core)<-------f(update) feedback loop to ZXMWE--------

c FCN_TAUVV_DIFFTRAN is passed as a parameter by ZXMWE (core optimisation routine).  
c In ZXMWE, updated 'candidates' nu and tau865 (from ZXMWD) are input into FCN_TAUVV_DIFFTRAN
c (=FCN=FUNCT calls in ZXMWD/E respectively).
c f output from FCN_TAUVV_DIFFTRAN is subsequently used in Hessian analysis and passed back to 
c ZXMWD where it is output along with final 'optimal' values of nu,tau865 once minimus is reached.

c input: N : number of input parameters;defined by N in ZXMWD which in turn is called
c            (and assigned) by NN in the ZXMWD call in SUBROUTINE DRV_TAUVV_ZXMWD.
c            current value is 2 (nu and tau865)
c        x : dimension = N; x(1)=retrieved nu, x(2)=retrieved tau865 (from ZXMWD).  In essence,
c            optnu[x(1)] and opttau [x(2)] become 'candidates' themselves for calculating
c            f = SD[rhoA(real)~rhoA(mod)] for subsequent ZXMWD analysis.            

c output: f : standard deviation of the master function, passed in the same manner as above.
c             Used as subsequent input into ZXMWD.


C Note: Regarding O2 absorption at band 765nm
C ----
C        t    c
C       R  = R  f     - correction to absorption by Rayleigh            (1)
C        r    r  r
C
C        c       t
C       R  = f  R     - correction to absorption by Aerosol             (2)
C        A    a  A
C
C        t         c                                        c
C       t   =  f  t   - correction to transmittance coeff (t = t t )    (3)
C        aer    0  aer                                            0
C
C so that, the modification to the original RTE:
C
C                given
C   (R   - f  R )      = R  / f   +  f  t t  R                          (4)
C     tot   r  r          A    a      0    0  w
C
C where t,t  are computed transmittances for given aerosol at theta,theta0
C          0
C
C Regarding <aerwatRho> variable:
C -------------------------------
C      In the main program NEWATM there was an assumption for Band-7
C                given
C   (R   - f  R )      = R  / f  , where aerwatRho(7) = R               (5)
C     tot   r  r          A    a                         A
C
C      If water reflectance is assumed as part of aerwatRho, (5) got to 
C      be modified to account for (4) through (for band-7):
C
C    R  = aerwatRho - f  f  t t  R                                      (6)
C     A                a  0    0  w



      IMPLICIT NONE

      INTEGER  N
      REAL     x(N)
      REAL     f

      INTEGER   NV,NMI,NMR,NLAMBDA
      PARAMETER (NMI=6,NMR=2,NV=6,NLAMBDA=8)
      INTEGER   NUM
      PARAMETER (NUM=75)

      REAL     theta0,theta,delt_phi
      REAL     air_mass,mu,mu0,pi
      REAL     vir(NMI,NMR)
      REAL     v(NV),vv_c,vstep
      REAL     cc78                             !calculated tau(765) for each (mi,mr)
      REAL     cc865tab(NLAMBDA,NMI,NMR,NV)     !tabulated tau(lambda)/tau(865)
      REAL     aerW0(NLAMBDA,NMI,NMR,NV), aerC(NLAMBDA,NMI,NMR,NV),
     &         s11dummy(NUM)

      REAL     aerwatRho_c(NLAMBDA)   !rhoA+rhow (same as aerwatRho)

      INTEGER  jvv,
     &         jji,jjr,ivBOT_c,ivBOT,ivMDL,ivTOP,ilow,jl,jv,jt
      REAL     ytmp(3),tauGrid(5),atmp(0:5),
     &         RoAtab(5,NLAMBDA,3),RoAvv(5)     !RhoA interpolated for geometry, interpolated for geometry and nu
      REAL     RhoA765,RhoA865,                 !calculated RhoA765,865
     &         Acoef(NLAMBDA),Bcoef(NLAMBDA),         !a,b,c,d HazeC coeffs interpolated for geometry and nu
     &         Ccoef(NLAMBDA),Dcoef(NLAMBDA),
     &         vAcoef(NLAMBDA,NV),vBcoef(NLAMBDA,NV), !a,b,c,d HazeC coeffs interpolated for geometry
     &         vCcoef(NLAMBDA,NV),vDcoef(NLAMBDA,NV)  

      REAL     fi670(NMI,NMR),fi765(NMI,NMR),fi865(NMI,NMR)   !'f' coefficients for centre band to band average conversion
       
      REAL     diffA_vri(NV,NMR,NMI,NLAMBDA), ! diffuse transmittance A coef pre-interpolated for Theta
     &         diffB_vri(NV,NMR,NMI,NLAMBDA), !                "      B                "
     &         diffA0_vri(NV,NMR,NMI,NLAMBDA),! diffuse transmittance A coef pre-interpolated for Theta0
     &         diffB0_vri(NV,NMR,NMI,NLAMBDA) !         "             B                "

c              diffuse transmittance coef (A/B) interpolated for theta and resultant nu (at mi,mr)
      REAL     diffA(NLAMBDA),diffB(NLAMBDA)

c              diffuse transmittance coef (A/B) interpolated for theta0 and resultant nu (at mi,mr)
      REAL     diffA0(NLAMBDA),diffB0(NLAMBDA)

      DATA     tauGrid / 0.05, 0.15, 0.30, 0.60, 0.80 /   !for tau lagrange polynomial fit

      INTEGER  isRELDEV
      REAL     xv,xt7,xt8,stddev7,stddev8
      REAL     xtau(NLAMBDA)
      LOGICAL  doLog

      INTEGER  iDiffTran,isDoCase2
      REAL     rhoA_true(NLAMBDA),rhoA_true_in(NLAMBDA),ttran(NLAMBDA),
     &         t(NLAMBDA),t0(NLAMBDA)
      REAL     fa,f0
      REAL     tau_ray(NLAMBDA),rho_wat_giv(NLAMBDA)
      REAL     FUNCT_OXYGEN_AER

      COMMON /SPEC_GEOM/   theta0,theta,delt_phi
      COMMON /COMPHASE/    aerW0,aerC,cc865tab,s11dummy
      COMMON /TAUVV/       jjr,jji,aerwatRho_c
      COMMON /TAUVV_MISC/  tau_ray,rho_wat_giv
      COMMON /CASE2_FLAG/  isDoCase2
      COMMON /V_MISC/      v,vv_c,vstep,ivBOT_c
      COMMON /V_ABCD/      vAcoef,vBcoef,vCcoef,vDcoef
      COMMON /V_EXTRA/     vir
      COMMON /COMDIFF_MISC/ diffA_vri ,diffB_vri,
     &                      diffA0_vri,diffB0_vri
      COMMON /F_INBAND/     fi670,fi765,fi865

 
C END$DECLARE

      doLog      = .TRUE.    !for log interpolation in quadratic interpolation

C SD master function control  
      isRELDEV   = 0     ! if isRELDEV=0 or 2 do absolute DIFF.FCN, otherwise (=1), STD.DEV.FCN (rel diff)
C------

      iDiffTran  = 2     ! 0 - no water-leaving radiance
                         ! 1 - simple Rayleigh,
                         ! 2 - exact diff-tran from LUT's

c switch off for testing.................
      IF(isDoCase2.EQ.0) iDiffTran=0   !assumes no Case 2 loop in main controlling 
c                                       routine SPEC1GEOFIT


C --------------------------GET AEROSOL INFORMATION------------------------------------------

      xv  = x(1)               !  xv = candidate nu input (from ZXMWD)...it is the CANDIDATE value here!!!!
      xt8 = x(2)               ! xt8 = candidate tau865 input (from ZXMWD).............".................

C -------------------SET UP NU MINI-TABLE FOR INTERPOLATION PURPOSES-----------------------------
C                     set up grid:ivbot,ivmin,ivtop from candidate nu

C     Locating the position of xv in the nu-table
c     select index number (of nu=v) that is closest to (but less than) xv
c     i.e., locate minimum index number (ilow) in v array that corresponds to 
c     candidate nu value xv: v[i]=1::6 ~ v = (2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
      CALL BRACKET4QUAD (NV,v,xv,ilow)

c     set up mini-table of nu indices.  Given nu candidate, ivBOT and ivMDL are index 
c     values of nu (v=1::6) that bound the candidate nu.  ivTOP completes the parabola.
c     note: from BRACKET4QUAD algorithm:
c     IF inputnu >= 4.5 THEN ilow=5 (=table nu value of 4.0, NOT 4.5)
c     IF inputnu <= 2.0 THEN ilow=1 (=table nu value of 2.0)
c     nu is out of interpolation range in both instances-could be dangerous (even for quadratic)
c     as it becomes an extrapolation!!
      IF(ilow+2.LE.NV) THEN  !eg., xv=optnu=3.7,ilow=4,ivbot=4,ivmid=5,ivtop=6 (v=3.5,4,4.5)
         ivBOT = ilow        !Hence first two points (ivbot,ivmid) bound candidate nu
         ivMDL = ivBOT+1
         ivTOP = ivBOT+2
      ELSE                   !eg., xv=optnu=4.4,ilow=5,ivbot=4,ivmid=5,ivtop=6 (v=3.5,4,4.5)
         ivBOT = ilow - 1    !Hence last two points (ivmid,ivtop) bound candidate nu
         ivMDL = ivBOT+1     ! if xv=4.7,ilow=5(from BRACKET4QUAD),ivbot=4,ivmid=5,ivtop=6
         ivTOP = ivBOT+2     ! Hence nu value is out of interpolation range!!!!!
      ENDIF

C -----------COMPUTING tau(765)/tau(865) FOR THE FOUND (VV=nu); mi,mr are given------------
C parameter tau(765)/tau(865); input grid: mi,mr,lambda,nu-mini-table value; interpolation for candidate nu

c     cc865tab = tabulated Cext(lambda)/Cext(865)=tau(lambda)/tau(865) for VV; 
c     Cext = extinction coeff

c     Initialise tabulated cc865tab [=tau(765)/tau(865)] for 3 consecutive nu 
c     (resultant nu-mini-table) index steps.
c     i.e., for grid point at [765,mi,mr,nu(mini-table)] we initialise cc865tab
c     at each  mini-table index (ivbot,ivmid,ivtop), with the indices determined by
c     the candidate nu input. That is, candidate nu is bounded by nu-mini-table values as previous
      ytmp(1)=cc865tab(7,jji,jjr,ivBOT)
      ytmp(2)=cc865tab(7,jji,jjr,ivMDL)
      ytmp(3)=cc865tab(7,jji,jjr,ivTOP)

c     quadratic interpolation in cc865tab nu-mini-table based on y(x)=ax^2+bx+c
c     result is cc78 (in 765,mr,mi) interpolated for candidate nu (from nu-mini-table)
c     Hence, cc78=tau(765)/tau(865) at (lambda,mr,mi) interpolated for nu (at lambda=765 of course)
      CALL QUADINTRP(doLog,
     &               xv,v(ivBOT),vstep,
     &               ytmp,cc78)

c     tau(765) =  tau(865) * [tau(765)/tau(865)] for each (mi,mr)
      xt7 = xt8 * cc78

C ------------------------------------------------------------------------------------------


C ---------- GIVEN (MR,MI) FINDING CANDIDATE REFLECTANCE(765,865) BASED ON A,B,C,D COEFF'S -----
C                                  FOR (NLAMBDA) & COMPUTED (CANDIDATE) VV
C                parameter rhoA(7,8); set up at grid: mi,mr,765::865,tau-table,nu-mini-table value

c finding RhoA(7,8) at each tau-table value, each (mi,mr) and each nu-mini-table value (ivbot,ivmid,ivtop)
C Rhoa taken from quartic fit fourier coefficients pre-interpolated for geometry
c see call of funct_a_b_c_d_fourier_interpolate in DRV_TAUVV_ZXMWD

      DO jv=1,3               
        jvv = ivBOT+jv-1      !ivbot,ivmid,ivtop

        DO jl=7,8             !lambda=765,865
          DO jt=1,5           !tau grid

c           Initialise RoAtab = tabulated Rho_A values at [tau-table,765::865,nu-mini-table] 
c           using respective fourier coeffs in larger [lambda,NV] set.   The nu-mini-table
c           is determined by the candidate nu input. Candidate nu is bounded by nu-mini-table
c           values as given earlier              
            RoAtab(jt,jl,jv) = vAcoef(jl,jvv)*tauGrid(jt)    +
     &                         vBcoef(jl,jvv)*tauGrid(jt)**2 +
     &                         vCcoef(jl,jvv)*tauGrid(jt)**3 +
     &                         vDcoef(jl,jvv)*tauGrid(jt)**4

          ENDDO ! jt
        ENDDO ! jl

      ENDDO  ! jv

C --- CALCULATING A,B,C,D COEFFICIENTS [quartic fit RhoA~tau(params)] FOR 765,865 ------------------
C parameter A,B,C,D; set up at grid: interpolated theta,mi,mr,765::865,interpolated candidate nu
C 1) using RhoA at grid point [mi,mr,tau_grid,765::865,nu-mini-table value], interpolate RhoA for candidate nu 
C 2) given RhoA at grid point [mi,mr,tau_grid,765::865,candidate nu], find polynomial fit (order 4) 
C    as a function of tau-table values.  This gives new A,B,C,D coeffs.

      DO jl=7,8               !lambda=765,865       
         DO jt=1,5            !tau-table

           DO jv=1,3          !ivbot,ivmid,ivtop

c            Prepare tabulated RoAtab for quadratic interpolation in candidate nu
             ytmp(jv)=RoAtab(jt,jl,jv)

           ENDDO

c          quadratic interpolation in RoAtab nu-mini-table based on y(x)=ax^2+bx+c
c          Result is RoAvv (in tau-table), pre interpolated for geometry
c          (funct_a_b_c_d_fourier_interpolate) and for nu right now! 
c          (using nu mini-table steps as input). We are still at the grid point (765,mi,mr) of course
           CALL QUADINTRP(doLog,xv,v(ivBOT),vstep,ytmp,RoAvv(jt))

         ENDDO

c        tau-table contains 5 points.  Hence we fit Lagrange Polynomial of degree 4,
c        with x-axis=tau-table values; y-axis = corresponding RoAvv (at each tau-table value) and for
c        interpolated geometry and candidate nu. 
c        We have effectively found a polynomial fit to
c        rhoa (params) = a(params)tau + b(params)tau^2 + c(params)tau^3 + d(params)tau^4      
c        The resultant equation is of the same order as the fourier quartic function
c        for LUT a,b,c,d coeffs in HazeC data.. 
c        Hence the Lagrange Polynomial coefficients obtained from this fit are our final 
c        interpolated coefficient values.   
         CALL LAGRINTRP(4,tauGrid,RoAvv,atmp)
         Acoef(jl)=atmp(1)
         Bcoef(jl)=atmp(2)
         Ccoef(jl)=atmp(3)
         Dcoef(jl)=atmp(4)
      ENDDO

C -------------------------- CALCULATING CANDIDATES RhoA(765::865)--------------------------------
C  parameter rhoA; final grid: mi,mr,interpolated geometry(implied by A,B,C,D),
C                              interpolated candidate nu [implied by tau(lambda)/tau(865)]

      RhoA765 = Acoef(7) * xt7    +
     &          Bcoef(7) * xt7**2 +
     &          Ccoef(7) * xt7**3 +
     &          Dcoef(7) * xt7**4

      RhoA865 = Acoef(8) * xt8    +
     &          Bcoef(8) * xt8**2 +
     &          Ccoef(8) * xt8**3 +
     &          Dcoef(8) * xt8**4

C ------------------------------------------------------------------------------------------


C ---------------------COMPUTING DIFFUSE TRANSMITTANCE COEFFICIENTS-------------------------
C           FOR ALL Mr,Mi (fixed), LAMBDA = 765-nm AND CANDIDATE (RETRIEVED) NU
c  A,B (for theta) and A0,B0 (for theta0) for 765-nm,Mr,Mi (fixed) and now interpolated for candidate
c  nu(mi,mr).  Mi,Mr and nu are now given of course: this subroutine is itself called in a (mi,mr) loop
c (see DRV_TAUVV_ZXMWD).

      IF (iDiffTran.NE.0) THEN  ! there is water-leaving radiance at Bands-7,8
        pi       = 4.0*atan(1.0)
        mu = COS(theta * pi / 180.0)
        mu0 = COS(theta0 * pi / 180.0)
        air_mass = 1./mu0 + 1./mu

        IF (iDiffTran.EQ.2) THEN ! compute aerosol transmittances from LUT's

C ------------------------- A-coeff (for theta)
          DO jl=7,8

c           Initialise (mr,mi,765 tabulated) diffA_vri for 3 consecutive nu 
c           (resultant nu-mini-table) index steps.
c           i.e., for grid point at [nu(mini-table),mr,mi,765::865] we initialise diffA_vri
c           at each  nu-mini-table index (ivbot,ivmid,ivtop), with the indices determined by
c           the candidate nu input. That is, candidate nu is bounded by nu-mini-table values as previous
            ytmp(1)=diffA_vri(ivBOT,jjr,jji,jl)
            ytmp(2)=diffA_vri(ivMDL,jjr,jji,jl)
            ytmp(3)=diffA_vri(ivTOP,jjr,jji,jl)

c           quadratic interpolation in diffA_vri nu-mini-table based on y(x)=ax^2+bx+c
c           result is diffA(in 765::865,mr,mi) interpolated for candidate nu (from nu-mini-table)
c           result was interpolated for geometry in SPEC1GEOFIT
            CALL QUADINTRP(.FALSE.,
     &                   vir(jji,jjr),v(ivBOT),vstep,
     &                   ytmp,diffA(jl))

C --------------------------B-coeff (for theta)
c same as above but for B-coeff

            ytmp(1)=diffB_vri(ivBOT,jjr,jji,jl)
            ytmp(2)=diffB_vri(ivMDL,jjr,jji,jl)
            ytmp(3)=diffB_vri(ivTOP,jjr,jji,jl)
            CALL QUADINTRP(.FALSE.,
     &                   vir(jji,jjr),v(ivBOT),vstep,
     &                   ytmp,diffB(jl))
          ENDDO ! jl=7,8


C --------------------------A0-coeff (for theta0)
c same as above (A and B) but for theta0
          DO jl=7,8


            ytmp(1)=diffA0_vri(ivBOT,jjr,jji,jl)
            ytmp(2)=diffA0_vri(ivMDL,jjr,jji,jl)
            ytmp(3)=diffA0_vri(ivTOP,jjr,jji,jl)
            CALL QUADINTRP(.FALSE.,
     &                   vir(jji,jjr),v(ivBOT),vstep,
     &                   ytmp,diffA0(jl))

C --------------------------B0-coeff (for theta0)

            ytmp(1)=diffB0_vri(ivBOT,jjr,jji,jl)
            ytmp(2)=diffB0_vri(ivMDL,jjr,jji,jl)
            ytmp(3)=diffB0_vri(ivTOP,jjr,jji,jl)
            CALL QUADINTRP(.FALSE.,
     &                   vir(jji,jjr),v(ivBOT),vstep,
     &                   ytmp,diffB0(jl))
          ENDDO ! jl=7,8

C ----------------------------------------------

          xtau(7) = xt7     !new allocation for tau(765) [at each (mi,mr)]
          xtau(8) = xt8     !new allocation for tau(865) [      "        ]

          DO jl=7,8
            t(jl)  = diffA(jl) *EXP(diffB(jl) * xtau(jl))   !diffuse transmittance at lambda 7,8 and theta (outbound)
            t0(jl) = diffA0(jl)*EXP(diffB0(jl)* xtau(jl))   !diffuse transmittance at lambda 7,8 and theta0 (inbound)
          ENDDO ! jl=7,8

          ttran(7) = t(7)*t0(7)   !t*t0 at 765
          ttran(8) = t(8)*t0(8)   !t*t0 at 865

        ELSE IF (iDiffTran.EQ.1) THEN ! simple Rayleigh Aerosol Transmittance

          ttran(7) = EXP(-tau_ray(7)*air_mass/2.0)   !pseudo t*t0 at 765
          ttran(8) = EXP(-tau_ray(8)*air_mass/2.0)   !pseudo t*t0 at 865

        ELSE ! assume Rayleigh diff-tran (iDiffTran.EQ.0.....not executed anyway)

          ttran(7) = EXP(-tau_ray(7)*air_mass/2.0)   !pseudo t*t0 at 765
          ttran(8) = EXP(-tau_ray(8)*air_mass/2.0)   !pseudo t*t0 at 865
          
        ENDIF

C --------------------------FINAL COMPARISON WITH TRUE DATA -----------------------------
c if iDiffTran = 1 or 2:
        fa = FUNCT_OXYGEN_AER(air_mass)   !converts aerosol at 765 with oxygen absorption to 
c                                          aerosol at 765 without oxygen absorption 
c        print*, theta, mu
c        print*, theta0, mu0
c        print*, air_mass
c        print*, faC

        f0 = 0.88-0.2*(air_mass-2.0)

c       compute rhoA            =[rhoTotal(real)-rhoRayleigh]-rhow
c       rho_wat_giv             = resultant rhow (at surface).  See Case 2 loop in SPEC1GEOFIT
c       fa*f0*ttran(7)*rho_wat_giv = rhow residual
c       Case 2 rhow at 7 and 8 removed here
c       band-averaged values
        rhoA_true_in(7) = aerwatRho_c(7) - fa*f0*ttran(7)*rho_wat_giv(7)  !fa*f0*ttran(7)*rho_wat_giv(7)
        rhoA_true_in(8) = aerwatRho_c(8) - ttran(8)*rho_wat_giv(8)

        CALL INBAND(rhoA_true_in(7),rhoA_true_in(8),
     &                    fi670(jji,jjr),fi765(jji,jjr),fi865(jji,jjr))

c       band-centred values
        rhoA_true(7) = rhoA_true_in(7)/fi765(jji,jjr)
        rhoA_true(8) = rhoA_true_in(8)/fi865(jji,jjr)

      ELSE   !if iDiffTran = 0

        CALL INBAND(aerwatRho_c(7),aerwatRho_c(8),
     &                    fi670(jji,jjr),fi765(jji,jjr),fi865(jji,jjr))

c       band-centred values
        rhoA_true(7) = aerwatRho_c(7)/fi765(jji,jjr)   !aerwatRho_c = rhoTotal-rhoRayleigh (main input)
        rhoA_true(8) = aerwatRho_c(8)/fi865(jji,jjr)   !              (see ATMCOR_SOA)        

      ENDIF   ! for whichever iDiffTran is selected

C ------------------------------------------ MASTER FUNCTION
c Remember: resultant RhoA765,865(at each mi,mr) is pre-interpolated for geometry and 
c interpolated now for candidate nu

      IF (isRELDEV.EQ.1) THEN
        stddev7 = RhoA765/rhoA_true(7) - 1.0
        stddev8 = RhoA865/rhoA_true(8) - 1.0
        f       = 100.0 * SQRT(stddev7**2 + stddev8**2)   !standard deviation function (relative differences)
      ELSE IF (isRELDEV.EQ.2) THEN
        stddev7 = RhoA765-rhoA_true(7)
        stddev8 = RhoA865-rhoA_true(8)
        f       = 100.0 * (ABS(stddev7) + ABS(stddev8))   !absolute differences
      ELSE
        stddev7 = RhoA765-rhoA_true(7)
        stddev8 = RhoA865-rhoA_true(8)
        f       = 10000.0 * SQRT(stddev7**2 + stddev8**2) !classical LSQ (absolute differences)
      ENDIF

      FCN_TAUVV_DIFFTRAN = .TRUE.

      RETURN
      END


C ############################################################################
C ############################################################################

C SECTION 2:

c This section of code performs the main spectral optimization of mi,mr,C,acdm,bbp.
c DRV_SPECOPT_LBFGSB = Drive SOA Large Broyden-Fletcher-Golfarb-Shanno Bounded.

c This subroutine contains the main optimization pre-processing (using channels 1 to 6) 
c and the call to controlling loop of main optimization (DRV_SPECOPT_LBFGSB).  
c It allows for all testing of contstraints and boundaries and can be removed later.

      SUBROUTINE MAIN_OPTIM_INOUT_PROC(IERopt,
     &                             optMI,optMR,optPIG,optACDM,optB0,
     &                             optVV,optTAU,optW0,optSTDDEV,
     &                             optTAUL,optTrans,optTrans0,
     &                             resAerRho,resWatRho,resTotRho,
     &                             fl,absCoef) 

C mi converted to mi^0.25 before optimisation input

      IMPLICIT NONE

C -- INPUT/PARAMETERS=SPECMATCH_GENERIC

      INTEGER  NMI,NMR,NV,
     &         NLAMBDA,NLAMBD0,
     &         NBEST,NBand

      INTEGER  i670,i765,i865
      PARAMETER (i670=6,i765=7,i865=8)

C -----------------------------------------------------------------------------
C File: SPECMATCH_GENERIC.PAR
C -----------------------------------------------------------------------------
      PARAMETER (NMI=6,NMR=2,NV=6)
      PARAMETER (NLAMBDA=8,NLAMBD0=2)
      PARAMETER (NBEST=8)
C -----------------------------------------------------------------------------

      REAL     fl(NLAMBDA)

      INTEGER  isCandData(NMI,NMR)
      REAL     aerwatRho(NLAMBDA),
     &         candRho(NLAMBDA,NMI,NMR)
      REAL     rcandRho(NMI,NMR,NLAMBDA),
     &         diffA_ri(NMR,NMI,NLAMBDA),diffB_ri(NMR,NMI,NLAMBDA),
     &         diffA0_ri(NMR,NMI,NLAMBDA),diffB0_ri(NMR,NMI,NLAMBDA)
      REAL     mi(NMI),mr(NMR)
      REAL     w0nu(NMI,NMR),tau865(NMI,NMR),cc865(NMI,NMR,NLAMBDA)
      REAL     vir(NMI,NMR)

      REAL     F0(NLAMBDA),mu0,pi

      INTEGER  IERopt
      REAL     optSTDDEV,optMI,optMR,optPIG,optACDM,optB0,optTAU,optW0,
     &         optVV,
     &         optTAUL(NLAMBDA),optTrans(NLAMBDA),optTrans0(NLAMBDA)
      REAL     resAerRho(NLAMBDA),resWatRho(NLAMBDA),
     &         resTotRho(NLAMBDA)

      REAL     absCoef(NLAMBDA),rhoW(NLAMBDA)

c      REAL     nLw_resid(NLAMBDA)   ! final nLw at each lambda (if required)

      INTEGER  ji,jr,jl

      INTRINSIC REAL

C -- EXTERNAL
      LOGICAL  DFCNLSQ_GRAD  !surrogate 'FCN' in main optizisation
      EXTERNAL DFCNLSQ_GRAD

      COMMON /BAND/        NBand
      COMMON /OPTFUNCREAL/ mi,mr,cc865,tau865,    !single COMMONS imported from DRV_TAUVV_ZXMWD
     &                     diffA_ri,diffB_ri,     !mi,mr imported from ATMCOR_SOA
     &                     diffA0_ri,diffB0_ri,
     &                     aerwatRho,candRho,
     &                     rcandRho            
      COMMON /CAND_MISC/   isCandData,w0nu
      COMMON /V_EXTRA/     vir
      COMMON /SPEC_MISC/   F0,mu0,pi
c      COMMON /STG2_MISC/   nLw_resid   !use if residual nLw at each lambda is required

      REAL*8   dmi_4(NMI),dmr(NMR),                                       !mi**0.25, mr
     &         dcandRho(NMI,NMR,NLAMBDA),                                 !RhoA
     &         ddiffA_ri(NMR,NMI,NLAMBDA),ddiffB_ri(NMR,NMI,NLAMBDA),     !diff trans A-coeff/B-coeff for theta
     &         ddiffA0_ri(NMR,NMI,NLAMBDA),ddiffB0_ri(NMR,NMI,NLAMBDA),   !diff trans A-coeff/B-coeff for theta0
     &         dcc865(NMI,NMR,NLAMBDA),dtau865(NMI,NMR),                  !tau(lambda)/tau(865),tau865
     &         dw0nu(NMI,NMR),dvir(NMI,NMR),                              !ssalbedo, nu
c    partial derivatives of above variables as a function of (mi,mr)
     &         pdCandRho(3,NMI,NMR,NLAMBDA),
     &         pdA(3,NMR,NMI,NLAMBDA),pdB(3,NMR,NMI,NLAMBDA),
     &         pdA0(3,NMR,NMI,NLAMBDA),pdB0(3,NMR,NMI,NLAMBDA),
     &         pdcc865(3,NMI,NMR,NLAMBDA),pdTau865(3,NMI,NMR),
     &         pdw0nu(3,NMI,NMR),pdvir(3,NMI,NMR)


C ------ NIR OPTIMIZATION - BOUNDARY CONSTRAINTS FOR NU & TAU865
      INTEGER   NNV
      PARAMETER (NNV=2)

      INTEGER     NN,n
      PARAMETER   (NN=5)
      REAL*8      xx(NN),l(NN),u(NN)
  
      COMMON /boundaries/ l,u  
      COMMON /dimension/ n

      n       = 5   !dimensions of problem (number of parameters)

C ------------------ SPECOPT-PREPROCESSING: PARTIAL DERIVATIVES OF TABULATED DATA -----------------------

c     declare mi,mr as doubles; mi converted to better conditioned form (**0.25)
      DO ji=1,NMI
        dmi_4(ji) = (DBLE(mi(ji)))**0.25D0
      ENDDO
      DO jr=1,NMR
        dmr(jr) = DBLE(mr(jr))
      ENDDO    
 
c-------- convert imported (COMMON) vars to double (from NIR optimization for each mr,mi combination) 
      DO ji=1,NMI
      DO jr=1,NMR
        DO jl=1,NLAMBDA  !lambda dependent vars imported from NIR optimization (DRV_TAUVV_ZXMWD)        
          ddiffA_ri(jr,ji,jl) =DBLE(diffA_ri(jr,ji,jl))    !    "     diffuse transmittance A coef interpolated for Theta and nu  
          ddiffB_ri(jr,ji,jl) =DBLE(diffB_ri(jr,ji,jl))    !    "     diffuse transmittance B coef interpolated for Theta and nu  
          ddiffA0_ri(jr,ji,jl)=DBLE(diffA0_ri(jr,ji,jl))   !    "     diffuse transmittance A coef interpolated for Theta0 and nu  
          ddiffB0_ri(jr,ji,jl)=DBLE(diffB0_ri(jr,ji,jl))   !    "     diffuse transmittance B coef interpolated for Theta0 and nu  
          dcc865(ji,jr,jl)    =DBLE(cc865(ji,jr,jl))       !    "     tau(lambda)/tau(865)
        ENDDO !non-spectral vars from DRV_TAUVV_ZXMWD

        !candidate rhoA(412::865) from NIR calculations and Haze C model interpolation
        DO jl=1,NLAMBDA    !bands 1 to 8
          dcandRho(ji,jr,jl)  =DBLE(rcandRho(ji,jr,jl))
        ENDDO
c       these are passed in from NIR commons (tau865,w0nu,vir)
        dtau865(ji,jr)      =DBLE(tau865(ji,jr))         !candidate tau at 865
        dw0nu(ji,jr)        =DBLE(w0nu(ji,jr))           !    "     ssalbedo at 865
        dvir(ji,jr)         =DBLE(vir(ji,jr))            !    "     nu       
      ENDDO
      ENDDO

c derive partial derivatives at each wavelength [z(x,y)---->zx,zy,zxy] 
c of all candidates as a function of mi,mr.
c [dim1,dim2,value(dim1),value(dim2),value(dim1*dim2)---->value(3*dim1*dim2)]
c note: first indices of multi dimensional arrays declared in call statements.  The loop iteration
c of all indices (mi,mr and respective vars) is performed in the called function (DRGPD3P).  
c Hence, an (mi,mr) loop is 'implicitly' nested inside the wavelength loop here.
      DO jl=1,NLAMBDA
c       do i,j=1 to mi(i),mr(j)
        CALL DRGPD3P(NMI,NMR,dmi_4,dmr,dcandRho(1,1,jl),     
     &               pdCandRho(1,1,1,jl))                                
        CALL DRGPD3P(NMI,NMR,dmi_4,dmr,dcc865(1,1,jl),
     &               pdcc865(1,1,1,jl))
        CALL DRGPD3P(NMR,NMI,dmr,dmi_4,ddiffA_ri(1,1,jl),
     &               pdA(1,1,1,jl))
        CALL DRGPD3P(NMR,NMI,dmr,dmi_4,ddiffB_ri(1,1,jl),
     &               pdB(1,1,1,jl))
        CALL DRGPD3P(NMR,NMI,dmr,dmi_4,ddiffA0_ri(1,1,jl),
     &               pdA0(1,1,1,jl))
        CALL DRGPD3P(NMR,NMI,dmr,dmi_4,ddiffB0_ri(1,1,jl),
     &               pdB0(1,1,1,jl))
      ENDDO
        CALL DRGPD3P(NMI,NMR,dmi_4,dmr,dtau865,pdTau865)
        CALL DRGPD3P(NMI,NMR,dmi_4,dmr,dw0nu,  pdW0nu)
        CALL DRGPD3P(NMI,NMR,dmi_4,dmr,dvir,   pdVir)         
c       enddo mi,mr

c bounds of optimization variables; l-lower, u-upper bound; note: mi = mi**0.25

      l(1)  = 0.D0    
      u(1)  = DBLE(mi(NMI))**0.25    ! w.r.t. mi4 (i.e., mi**0.25)

      l(2)  = DBLE(mr(1))            ! mr
      u(2)  = DBLE(mr(NMR))        
        
      l(3)  = 0.02D0                   ! Chlorophyll concentration
      u(3)  = 60.D0                   !
  
      l(4)  = 0.0008D0                ! Acdm(443) 
      u(4)  = 1.1D0                  !
    
      l(5)  = 0.0003D0                 ! Bbp(443)
      u(5)  = 0.1D0   
  
      CALL DRV_SPECOPT_LBFGSB(NBand,aerwatRho,
     &                          dmi_4,dmr,
     &                          dcandRho,pdCandRho,
     &                          dtau865,pdTau865,
     &                          dvir,pdvir,            !nu for internal iterations only
     &                          dw0nu,pdW0nu,          !w0 for internal iterations only
     &                          dcc865,pdcc865,
     &                          ddiffA_ri, pdA ,ddiffB_ri ,pdB ,
     &                          ddiffA0_ri,pdA0,ddiffB0_ri,pdB0,
     &                          F0,mu0,pi,
     &                          xx,
     &                          optMI,optMR,optPIG,optACDM,
     &                          optB0,ierOPT,optSTDDEV,fl) 

      IF (IERopt.EQ.136) GO TO 911   ! non-computable FCN (=DFCNLSQ_GRAD), stop process
      
      CALL SPECOPT_POSTPROC (N,xx,aerwatRho,
     &                         dmi_4,dmr,
     &                         dcandRho,pdCandRho,
     &                         dtau865,pdTau865,dcc865,pdcc865,
     &                         dw0nu,pdW0nu,dvir,pdvir,
     &                         ddiffA_ri ,pdA ,ddiffB_ri ,pdB ,
     &                         ddiffA0_ri,pdA0,ddiffB0_ri,pdB0,
     &                                               F0,mu0,pi,
     &                     optW0,optVV,
     &                     optTrans0,optTrans,
     &                     optTAU,optTAUL,
     &                     resAerRho,resWatRho,resTotRho,
     &                     absCoef  )
           

c not output here but inserted for checking mechanism if required
      DO jl=1,NLAMBDA
c              nLw_resid(jl) =
c     &         ((F0(jl)*mu0)/pi)*
c     &           (aerwatRho(jl)/optTrans(jl)/optTrans0(jl) -
c     &            resAerRho(jl)/optTrans(jl)/optTrans0(jl)   )  ! avoid small values

        rhoW(jl) = resWatRho(jl)*optTrans(jl)*optTrans0(jl) !at sensor
      ENDDO

c derivation of 'residual' normalised water leaving radiance nLw AT SURFACE (division by t and t0).
c since [(rhot-rhoray)/t/t0] ~ real (rhoA+rhow)/t/t0,
c [(rhot-rhoray)/t/t0](real)  - rhoA(retieved)/t/t0 =  rhow
c rhow*f0*cos(theta)/pi= nLw at surface

C END$RUN

 911  CONTINUE

      RETURN
      END

C*********************************************************************************

c DRV_SPECOPT_LBFGSB = Drive SOA Large Broyden-Fletcher-Golfarb-Shanno Bounded.
c Controlling routine for optimization of SeaWiFS bands 1 to 6.
c Startups: 
c for high ACDM waters and Chl 0.1 to 1.0 (mgm-3), startups require 0.1 Chl resolution.
c Tolerances:
c gradient and LSQ tolerances are set
c lsq tolerance contollled by factr, a double precision variable.
c       On entry factr >= 0 is specified by the user.  The iteration
c         will stop when
c
c         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c         where epsmch is the machine precision, which is automatically
c         generated by the code. Typical values for factr: 1.d+12 for
c         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
c         high accuracy.
c       On exit factr is unchanged.
c Stopping filters:
c designed to halt optimization at best local minimum. Purpose: unwise to run optimisation 
c until fact (1.d+7) is reached.  This results in unrealistic fit (rhoA+rhow)~(rhot-rhor).
c fact = 1.d+12 unwise.  Local mimimum not always reached = inconsistency across pixels.

c call to master function and gradient calcs (FUNCTION DFCNLSQ_GRAD)

      SUBROUTINE DRV_SPECOPT_LBFGSB(NBand,aerwatRho,
     &                          dmi_4,dmr,
     &                          dcandRho,pdCandRho,
     &                          dtau865,pdTau865,
     &                          dvir,pdvir,            !nu for internal iterations only 
     &                          dw0nu,pdW0nu,          !w0 for internal iterations only
     &                          dcc865,pdcc865,
     &                          ddiffA_ri ,pdA ,ddiffB_ri ,pdB ,
     &                          ddiffA0_ri,pdA0,ddiffB0_ri,pdB0,
     &                          F0,mu0,pi,
c output:
     &                          xx,
     &                          optMI,optMR,optPIG,optACDM,
     &                          optB0,ierOPT,optSTDDEV,fl)


c--input: none.  All input is from COMMONS
c        
c--output:
c         ierOPT         -integer error-flag for Newt-Raph search (for inside world)
c                         *error status is converted to logical form (1,0) in SPEC1GEOFIT (by IERR declaration)
c         optMI          -optimized imaginary part of refractive index
c         optMR          -optimized real part of refractive index
c         optPIG         -optimized chlorophyll concentration (mg/l)
c         optACDM        -optimized cdm absorption at 443-nm (= agelb+adet)
c         optB0          -optimized particulate backscattering coeff.at 443 nm
c         optVV          -optimized aerosol Junge distribution parameter (nu)
c         optTAU         -optimized aerosol optical thickness at 865-nm
c         optW0          -optimized single scattering albedo
c         optSTDDEV      -standard deviation of master function in main optimization (DRV_SPECOPT_LBFGSB)
c         optTAUL        -optimized aerosol optical thickness at lambda (each wavelength)
c         optTrans       -diffuse transmittance outbound (at theta)
c         optTrans0      -diffuse transmittance inbound (at theta0)
c         resAerRho      -resultant aerosol reflectance rhoA (= aerosol reflectance, prime)
c         resWatRho      -resultant water reflectance rhow (= water reflectance at surface, prime)
c         resTotRho      -resultant total reflectance rhot (= total reflectance at TOA, prime)

      IMPLICIT NONE

C -- INPUT/PARAMETERS=SPECMATCH_GENERIC

      INTEGER  NMI,NMR,NV,
     &         NLAMBDA,NLAMBD0,
     &         NBand

C -----------------------------------------------------------------------------
C File: SPECMATCH_GENERIC.PAR
C -----------------------------------------------------------------------------
      PARAMETER (NMI=6,NMR=2,NV=6)
      PARAMETER (NLAMBDA=8,NLAMBD0=2)
C -----------------------------------------------------------------------------
     
      INTEGER  isCandData(NMI,NMR)
      REAL     aerwatRho(NLAMBDA)
     
      REAL     w0nu(NMI,NMR)
      REAL     vir(NMI,NMR)

      REAL     F0(NLAMBDA),mu0,pi

      INTEGER  IERopt
      REAL     optSTDDEV,
     &         optMI,optMR,optPIG,optACDM,optB0

      INTEGER   NN,NM
      PARAMETER(NN=5,NM=5)

      INTEGER  j

      INTRINSIC REAL
 
C ---------------------------------------------------------------------------
C parameter startups
C ---------------------------------------------------------------------------
      INTEGER          f_start,f_fin,count
      REAL             fBest,f_avl,f_temp
     
c     flags for best starting-point selection (of mi,mr,C,acdm,bbp input)
      INTEGER          do12best,bestiter 
      INTEGER          idxf,jidx,N6BEST
      PARAMETER        (N6BEST=6)           !The #of starting points provided in the
                                              !arrays <cb>,<acdmb>,<bbpm>;             
      REAL*8           mi4b(N6BEST),mrb(N6BEST),
     &                 cb(N6BEST),acdmb(N6BEST),bbpb(N6BEST)
    
C ---------------------------------------------------------------------------

C -- EXTERNAL
      LOGICAL  DFCNLSQ_GRAD  !surrogate 'FCN' in main optizisation
      EXTERNAL DFCNLSQ_GRAD

      COMMON /CAND_MISC/   isCandData,w0nu
      COMMON /V_EXTRA/     vir
c      COMMON /SPEC_MISC/   F0,mu0,pi
c      COMMON /STG2_MISC/   nLw_resid   !use if residual nLw at each lambda is required

      REAL*8   dmi_4(NMI),dmr(NMR),                                          !mi**0.25, mr
     &         dcandRho(NMI,NMR,NLAMBDA),                                    !RhoA
     &         ddiffA_ri(NMR,NMI,NLAMBDA),ddiffB_ri(NMR,NMI,NLAMBDA),        !diff trans A-coeff/B-coeff for theta
     &         ddiffA0_ri(NMR,NMI,NLAMBDA),ddiffB0_ri(NMR,NMI,NLAMBDA),      !diff trans A-coeff/B-coeff for theta0
     &         dcc865(NMI,NMR,NLAMBDA),dtau865(NMI,NMR),                     !tau(lambda)/tau(865),tau865
c    partial derivatives of above variables as a function of (mi,mr)
     &         pdCandRho(3,NMI,NMR,NLAMBDA),
     &         pdA(3,NMR,NMI,NLAMBDA),pdB(3,NMR,NMI,NLAMBDA),
     &         pdA0(3,NMR,NMI,NLAMBDA),pdB0(3,NMR,NMI,NLAMBDA),
     &         pdcc865(3,NMI,NMR,NLAMBDA),pdTau865(3,NMI,NMR),
     &
     &         dvir(NMI,NMR)           ,pdVir(3,NMI,NMR),
     &         dw0nu(NMI,NMR)          ,pdW0nu(3,NMI,NMR)

      REAL     tau865_iter,v_iter,w0_iter
      REAL     trans(NLAMBDA),trans0(NLAMBDA),tau(NLAMBDA)
      REAL     resAerRho(NLAMBDA),resWatRho(NLAMBDA),
     &         resTotRho(NLAMBDA)

C    STARTING CONDITIONS Case 2: Kuchinke

      DATA mi4b / 0.01D0,0.01D0,0.01D0,
     &            0.01D0,0.01D0,0.01D0/,

     &      mrb / 1.34D0,1.34D0,1.34D0,
     &            1.34D0,1.34D0,1.34D0/,
   
     &       cb / 0.1D0,0.1D0,
     &            1.0D0,1.0D0,1.0D0,1.0D0/, 

     &    acdmb / 0.001D0,0.01D0,  
     &            0.001D0,0.01D0,0.1D0,0.5D0/,  
  
     &     bbpb / 0.01D0,0.01D0,
     &            0.001D0,0.005D0,0.01D0,0.03D0/ 


      CHARACTER*60     task, csave
      LOGICAL          lsave(4)
      INTEGER          n, m, iprint,
     +                 nbd(NN), iwa(3*NN), isave(44)
      REAL*8           ff, factr, pgtol,
     +                 xx(NN),l(NN), u(NN), gg(NN), dsave(29), 
     +                 wa(2*NM*NN+4*NN+12*NM*NM+12*NM)

      REAL*8           f_lambda(NLAMBDA)
      REAL             fl(NLAMBDA)

      INTEGER      uf           ! history file unit
      LOGICAL      lHistory     ! do output for history?
      CHARACTER*20 fHistName    ! history file name

      COMMON /boundaries/ l,u        
      COMMON /dimension/ n

      REAL lat
      COMMON /latitude/ lat  

c     note: epsmch = 2.220446049x10-16 for current machine

C ------------------------------ SPECOPT-INITIALIZATION; Large BFGS bounded input ---------------------------------------

      lHistory  = .false.  !.true.        !no output of optimization parameters = FALSE
      fHistName = 'lbfgs.out'    !file name of 'optimization output' data
    
      IF (lHistory) CALL ZXMWH_INI (uf,fHistName)

c LBFGSB computation settings (subroutine setulb)
c      n       = 5   !dimensions of problem; now declared in SUBROUTINE MAIN_OPTIM_INOUT_PROC

c maximum number of variable metric corrections used to define the limited memory matrix
      m       = 5 

c no iteration-summary output is generated for iprint < 0
      iprint  = -1  !no iteration-summary output is generated for iprint < 0

c error status of from DFCNLSQ_GRAD call (ff,gg computation)
      IERopt  = 0

c tolerances on corresponding steps/iterations in optimization (see main comments SPEC1GEOFIT)     
c factr tells iteration to stop as a function of machine precision (auto generated by code)
c Typical values for factr: 1.d+12 = low accuracy; 1.d+7 = moderate accuracy; 1.d+1 for very high accuracy.     
      factr   = 1.0d+7   !DO NOT CHANGE

c The iteration will stop when max{|proj g_i | i = 1, ..., n} <= pgtol
c where pg_i is the ith component of the projected gradient. 
      pgtol   = 1.0d-5    !DO NOT CHANGE 

c BFGSB input nbd is an integer array of dimension n. On startup, nbd = type of bounds imposed on the variables:
c nbd(i)=0 if x(i) is unbounded,
c        1 if x(i) has only a lower bound,
c        2 if x(i) has both lower and upper bounds, and
c        3 if x(i) has only an upper bound.
      DO j=1,NN
        nbd(j)=2  ! x(i) has both lower and upper bounds
      ENDDO

c range of flsq iterations for 'best' startup selection
      f_start = 5
      f_fin   = 15

c To use best starting points (=.TRUE.) or to allow the manual selection by the
c user (=.FALSE.). If .TRUE., see the variables below.
      do12best       =  1
      bestiter       =  0
   
c index of current starting point (1::N6Best) for each mi,mr,C,acdm,bbp 
      jidx      = 1   !index of counter 'j'

c index of manually selected startup values        (if Do12best is INITIALLY set to FALSE)
c index of final (best SD producing) startup values (if Do12best is INITIALLY set to TRUE)
c (initialisation created by kuchinke: error in original Chomko code)
      idxf      = 1   !index 'final' OR index corresponding to 'best' or 'selected' f (=SD)

c dummies in startup selection;   stores 'best' isave(30)
c comparisons (with other candidate startups)
      fBest    = 1000.

C ----------------------- SOLVING A SYSTEM OF BOUNDED NON-LINEAR EQUATIONS ---------------------------------
C                            WITH LSQ-OBJECTIVE (commence optimization)
   
C *********LOOP 100 START: SELECTION OF OPTIMUM STARTUP-VALUES FOR MI,MR,C,ACDM,BBP************

c If Do12Best = TRUE: This loop inputs the first set of startup values into L-BFGS-B and calculates SD.  
c It them moves to the next startup values and calculates revised SD.  If revised SD is lower than previous, 
c then original startup values are replaced by new ones.  Loop continues until all optimization SD's from 
c all startup values have been compared. 
c If Do12Best = FALSE: The manually selected index (idxf setting) is used to select one set of start up values

 100  CONTINUE

        IF (do12best.eq.1) THEN    !input first group of startup (initialisation=approx solution) values
          xx(1) = mi4b(jidx)    
          xx(2) = mrb(jidx)     
          xx(3) = cb(jidx)     
          xx(4) = acdmb(jidx)     
          xx(5) = bbpb(jidx)   
        ELSE                
          xx(1) = mi4b(idxf) !final start up values to be used (best SD has been selected)
          xx(2) = mrb(idxf)  !if Do12best=false (initially) then this is the manually selected startup.
          xx(3) = cb(idxf)
          xx(4) = acdmb(idxf)
          xx(5) = bbpb(idxf)
        ENDIF
          
      task='START'
       
c------ call the L-BFGS-B code and  SOLVE A SYSTEM OF BOUNDED NON-LINEAR EQUATIONS IN mi,mr,C,acdm and bbp.

c Implicit call to logical function DFCNLSQ_GRAD (next subroutine).
 
c DRV_TAUVV_DIFFTRAN: a,b start-up constraints, N 
c a and b are the constraint (startup) vectors for nu and tau (length = NN = 2)
c # of unknown parameters = NN = 2 = nu,tau865.

c Process:
c DRV_SPECOPT_LBFGSB ---xx,l,u,---> Setulb ---> mainlb (core BFGS) ---N=5,optimal 'candidates' xx---> DFCNLSQ_GRAD 
c                    <--final x(1)..x(5),f----                    <------ff,gg(update) feedback loop to mainlb----

c setulb: partitions the working arrays wa and iwa, and then uses the limited memory BFGS method to solve the 
c bound constrained optimization problem by calling mainlb. (setulb='set upper-lower boundary'?)

c mainlb: this subroutine solves bound constrained optimization problems by using the compact formula of the 
c limited memory BFGS updates (mainlb='main library'?)
c DFCNLSQ_GRAD (='double function least-squares/gradient' or 'derivative function least-squares/gradient'? )

c xx (length N) is an input array of the initial 'candidate' values for mi,mr,C,acdm,bbp

c ff is the standard deviation of the master function, a comparison of 
c [rhot(real)-rhoray]~[rhoA(retrieved)+rhow(retrieved)] in DFCNLSQ_GRAD
c Note: On first entry ff is unspecified.  On final exit f is the value of the function at xx. 

c gg is also computed by DFCNLSQ_GRAD. It is the gradient with respect to each of the 5 input parameters
c C,acdm,bbp,mi,mr.  On first entry g is unspecified. On final exit g is the value 
c of the gradient at x.

c output variable 'task' from setulb provides call for ff,gg retrieval (call to DFCNLSQ_GRAD)
c final values of xx(1::N) are reached when mainlb computes minimus

C CALL TO L-BFGS-B

c###########################################################################

 111  continue

      call setulb(n,m,xx,l,u,nbd,ff,gg,factr,pgtol,wa,iwa,task,
     +            iprint,csave,lsave,isave,dsave)


c      input: n,m,xx,l,u,nbd,ff,gg,factr,pgtol
c      output: wa,iwa,task,iprint,csave,lsave,isave,dsave

CCCCCCCCCCC explanation of optimization outputs (see setulb and mainlb code)

c     wa is a double precision working array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax.
c
c     iwa is an integer working array of length 3nmax.
c
c     task is a working string of characters of length 60 indicating
c       the current job when entering and quitting this subroutine.

c     csave is a working string of characters of length 60.
c
c     lsave is a logical working array of dimension 4.
c       On exit with 'task' = NEW_X, the following information is 
c                                                             available:
c         If lsave(1) = .true.  then  the initial X has been replaced by
c                                     its projection in the feasible set;
c         If lsave(2) = .true.  then  the problem is constrained;
c         If lsave(3) = .true.  then  each variable has upper and lower
c                                     bounds;
c
c     isave is an integer working array of dimension 44.
c       On exit with 'task' = NEW_X, the following information is 
c                                                             available:
c         isave(22) = the total number of intervals explored in the 
c                         search of Cauchy points;
c         isave(26) = the total number of skipped BFGS updates before 
c                         the current iteration;
c         isave(30) = the number of current iteration;
c         isave(31) = the total number of BFGS updates prior the current
c                         iteration;
c         isave(33) = the number of intervals explored in the search of
c                         Cauchy point in the current iteration;
c         isave(34) = the total number of function and gradient 
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         if isave(37) = 0  then the subspace argmin is within the box;
c         if isave(37) = 1  then the subspace argmin is beyond the box;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints in the current
c                         iteration;
c         n + 1 - isave(40) = the number of variables leaving the set of
c                           active constraints in the current iteration;
c         isave(41) = the number of variables entering the set of active
c                         constraints in the current iteration.
c
c     dsave is a double precision working array of dimension 29.
c       On exit with 'task' = NEW_X, the following information is
c                                                             available:
c         dsave(1) = current 'theta' in the BFGS matrix;
c         dsave(2) = f(x) in the previous iteration;
c         dsave(3) = factr*epsmch;
c         dsave(4) = 2-norm of the line search direction vector;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(7) = the accumulated time spent on searching for
c                                                         Cauchy points;
c         dsave(8) = the accumulated time spent on
c                                                 subspace minimization;
c         dsave(9) = the accumulated time spent on line search;
c         dsave(11) = the slope of the line search function at
c                                  the current point of line search;
c         dsave(12) = the maximum relative step length imposed in
c                                                           line search;
c         dsave(13) = the infinity norm of the projected gradient;
c         dsave(14) = the relative step length in the line search;
c         dsave(15) = the slope of the line search function at
c                                 the starting point of the line search;
c         dsave(16) = the square of the 2-norm of the line search
c                                                      direction vector.
CCCCCCCCCCC

c       OPEN(unit=63,file='task.txt',status='UNKNOWN',
c     &        form='FORMATTED',access='APPEND')
c         write(63,*) 'task =', task
c       CLOSE (63)

C SECTION 1: START ITERATION
c     the minimization routine has returned to request the
c     function ff and gradient gg values at the current xx. 

      
      if (task(1:2) .eq. 'FG') then

c        compute function value ff for the sample problem.
        IF (.NOT.DFCNLSQ_GRAD (N,xx,NBand,
     &                         aerwatRho,
     &                         dmi_4,dmr,
     &                         dcandRho,pdCandRho,
     &                         dtau865,pdTau865,
     &                         dvir,pdvir,            !input pdnu for internal iterations only 
     &                         dw0nu,pdW0nu,          !input pdw0 for internal iterations only 
     &                         dcc865,pdcc865,
     &                         ddiffA_ri,pdA,ddiffB_ri,pdB,
     &                         ddiffA0_ri,pdA0,ddiffB0_ri,pdB0,
     &                         F0,mu0,pi,
     &                         tau865_iter,v_iter,w0_iter,    !outputs for internal iterations only 
     &                         trans0,trans,                  !outputs for internal iterations only   
     &                         tau,                           !outputs for internal iterations only  
     &                         resAerRho,resWatRho,resTotRho, !outputs for internal iterations only  
     &                         ff, gg , f_lambda)) THEN

          IERopt = 136
          task   = 'STOP'
          !can't compute, go back to the minimization routine and try again
          GO TO 111
        ENDIF

c       finished retrieving current ff,gg updates, re-enter minimization routine.
        GO TO 111
      ENDIF


C SECTION 2: TEST FOR NEW ITERATION     
c     Convergence criteria met for current x-values.  Try new direction,x-point?

      IF (task(1:5) .eq. 'NEW_X')  THEN

c          All optim.parameters can be output here.
c          Set to .FALSE. to have less garbage on the disk.
        IF (lHistory)
     &     CALL DZXMWH_WRITE (lat,uf,NN,isave(30),ff,xx,gg,
     &                        f_lambda,tau865_iter,v_iter,w0_iter,
     &                        aerwatrho(1),aerwatrho(2),aerwatrho(3),
     &                        aerwatrho(4),aerwatrho(5),aerwatrho(6),
     &                        jidx, idxf)

c###########

        IF (do12best .eq. 1) THEN            
          IF (isave(30) .eq. f_start) then
            f_avl  = real(ff)  
            count=2
          ENDIF  
          IF (isave(30) .gt. f_start .and. isave(30) .le. f_fin) THEN    !range 5 to 14
            f_temp  = real(ff)         
            f_avl  = (real(count-1)*f_avl/real(count) 
     &                                     + f_temp*1./real(count))                     
            count = count+1
          ENDIF
          IF (isave(30) .ge. f_fin+1) goto 167                 !exit at 15th iteration
        ENDIF   !do12best
 
        GO TO 111     !continue iterations
 
      ENDIF                                          !NEW_X

c########

      IF ((do12best.eq.1) .and. lHistory) go to 167    !continue 'write' if startup  not forced to stop early

      IF (lHistory) CALL ZXMWH_END (uf)

 167   continue

c###########################################################################        
c     If 'use custom starting points' and current starting point index <= total num starting sets

      IF ((do12best.eq.1) .AND. jidx .LE. N6BEST) THEN

cif: final SD from current starting set < large dummy (first loop) OR stored SD < final SD from 
cprevious loop (previous start ups), then
c note:    f_avl = average flsq for iterations f_start to f_fin
              
          IF (f_avl .LT. fBest) THEN        !f_avl or bestcount
            !update with 'revised' SD from current starting set     
            fbest = f_avl               
            !store index of 'best' startups
            idxf  = jidx               !overwrites initial idxf with best 'f' startup index
          ENDIF

        ! if current loop has reached last startups, stop checking
        IF (jidx .EQ. N6BEST) THEN
          do12best = 0
          bestiter = 1
        ! else move to next set of startups
        ELSE
          jidx = jidx+1
        ENDIF
        GOTO 100
      ENDIF

c###########################################################################

C ***********LOOP 100 END**************

c mi converted back to original form [ (mi^0.25)^4 ]
      optMI     = REAL(xx(1))**4   ! w.r.t. mi4
      optMR     = REAL(xx(2))
      optPIG    = REAL(xx(3))
      optACDM   = REAL(xx(4))
      optB0     = REAL(xx(5))
      optSTDDEV = REAL(ff)

      DO j = 1,NBAND 
        fl(j)=abs(real(f_lambda(j)))    !rhot-rhor v rhoA+t*rhow absolute comparison
      ENDDO

      RETURN
      END

C *****************************************************************************

C Master function routine for main optimization.

C Note: Any change to Master Function will require modifications to exact derivatives
C in respective called routines here

C This subroutine takes the startup (first run - optimization) or current (subsequent runs - optimization)
C values of mi,mr,C,acdm,bbp as input.
C It then inteprolates the following for mi and mr: tau865, rhoA(lambda),rhow(lambda),tau(lambda)/tau(865) and
C diffuse transmittance  A/B coeff at theta/theta0.
C It then computes tau(lambda), diffuse transmittance (outbound;inbound) at each lambda and MODELLED(rhoA+rhow)
C All are therefore interpolated for mi,mr (and nu from NIR optimization of course)
C note: diffuse transmittance coeffs have also been interpolated for theta/theta0 in SPEC1GEOFIT.
C The LSQ function for LBFGSB input is finally calculated, as is the gradient of each mi,mr,C,acdm,bbp
C The gradient is defined as D_function_LSQ/D_parameter.  The respective partial derivatives
C D_[rhoA+rhow] / D_(mi,mr,C,acdm,bbp) and D_LwN / D_(C,acdm,bbp) required for gradients are 
C calculated in two (called) derivative subroutines.

      LOGICAL FUNCTION DFCNLSQ_GRAD (N,x,NBand, !Double function LSQ and Grad
     &                               aerwatRho,    
     &                               dmi_4,dmr,
     &                               dcandRho,pdCandRho,
     &                               dtau865,pdTau865,
     &
     &                               dvir,pdvir,            !input pdnu for internal iterations only 
     &                               dw0nu,pdW0nu,          !input pdw0 for internal iterations only 
     &
     &                               dcc865,pdcc865,
     &                               ddiffA_ri ,pdA ,ddiffB_ri ,pdB ,
     &                               ddiffA0_ri,pdA0,ddiffB0_ri,pdB0,
     &                               F0,mu0,pi,
c output
     &                               tau865_iter,v_iter,w0_iter,      !outputs for internal iterations only 
     &                               trans0,trans,                    !outputs for internal iterations only   
     &                               tau,                             !outputs for internal iterations only  
     &                               resAerRho,resWatRho,resTotRho,   !outputs for internal iterations only  
     &                               f, g, f_lambda)      !f_lambda:  for internal iterations only  

c input: 
c           N - number of input parameters; current value is 5 (C,acdm,bbp,mi,mr)
c           x - (length N) is an input array of the initial 'candidate' values for mi,mr,C,acdm,bbp
c       Nband - number of bands
c   aerwatRho - real total reflectance - rayleigh reflectance [rhot-rhoray]
c       dmi_4 - mi^4 (double); 6 values
c         dmr - mr (double) 2 values
c    dcandRho - candidate RhoA (at each lambda,mi,mr) interpolated for geometry (SPEC1GEOFIT) 
c               and interpolated for candidate nu (optimal nu given mi,mr in DRV_TAUVV_ZXMWD)
c   pdCandRho - partial derivative of dcandRho(lambda) [=rhoA(lambda)] as a function of mi,mr
c     dtau865 - candidate tau865 at each mi,mr = optimal tau865 GIVEN mi,mr
c    pdTau865 - partial derivative of tau865 as a function of mi,mr
c      dcc865 - candidate tau(lambda)/tau(865) =  optimal tau(lambda)/tau865 GIVEN mi,mr
c               and interpolated for candidate nu (= interpolated for optimal nu GIVEN mi,mr)
c     pdcc865 - partial derivative of tau(lambda)/tau865 as a function of mi,mr, interpolated for 'candidate' nu

c   ddiffA_ri - candidate diffuse transmittance A coeff interpolated for Theta and candidate nu
c               (= optimal nu GIVEN mi,mr)
c         pdA - partial derivative of ddiffA_ri as a function of mi,mr, interpolated for 'candidate' nu
c   ddiffB_ri - as above for B coeff
c         pdB - as above for B coeff
c  ddiffA0_ri - as above for A coeff and Theta0
c        pdA0 - as above for A coeff and Theta0
c  ddiffB0_ri - as above for B coeff and Theta0
c        pdB0 - as above for B coeff and Theta0

c   F0,mu0,pi - passed explicitly here.  Input to DRV_SPECOPT_LBFGSB by implicit COMMON /SPEC_MISC/

c output:
c           f - (or ff) Standard deviation of the master function; a comparison of 
c               [rhot(real)-rhoray]~[rhoA(retrieved)+rhow(retrieved)] in DFCNLSQ_GRAD
c                Note: On first entry ff is unspecified.  On final exit f is the value of the function at xx. 
c           g - (or gg) Computed by DFCNLSQ_GRAD. On first entry g is unspecified. On final exit g is the value 
c               of the gradient at x.

c Note: z prefix denotes interpolated values of respective variables (in mi,mr subspace)

      IMPLICIT NONE

      INTEGER  N
      REAL*8   x(N)
      REAL*8   f,g(N)

      INTEGER  NMI,NMR,NLAMBDA,NBand
      PARAMETER (NMI=6,NMR=2,NLAMBDA=8)
      INTEGER  iFcnLSQ                     !SD form

      INTEGER    NMAX
      PARAMETER (NMAX=5)

      REAL     F0(NLAMBDA),mu0,pi

      REAL     aerwatRho(NLAMBDA)
      REAL     tau865_iter
      REAL     v_iter,w0_iter                 !outputs for internal iterations only 
      REAL     trans(NLAMBDA),trans0(NLAMBDA),tau(NLAMBDA)
      REAL     resAerRho(NLAMBDA),resWatRho(NLAMBDA),
     &         resTotRho(NLAMBDA)

      REAL*8   dmi_4(NMI),dmr(NMR)
      REAL*8   dcandRho(NMI,NMR,NLAMBDA),pdCandRho(3,NMR,NMI,NLAMBDA),
     &         ddiffA_ri(NMR,NMI,NLAMBDA) ,pdA(3,NMR,NMI,NLAMBDA) ,
     &         ddiffB_ri(NMR,NMI,NLAMBDA) ,pdB(3,NMR,NMI,NLAMBDA) ,
     &         ddiffA0_ri(NMR,NMI,NLAMBDA),pdA0(3,NMR,NMI,NLAMBDA),
     &         ddiffB0_ri(NMR,NMI,NLAMBDA),pdB0(3,NMR,NMI,NLAMBDA),
     &         dcc865(NMI,NMR,NLAMBDA) ,pdcc865(3,NMI,NMR,NLAMBDA),
     &         dtau865(NMI,NMR)        ,pdTau865(3,NMI,NMR),
     &         dvir(NMI,NMR)           ,pdvir(3,NMI,NMR),
     &         dw0nu(NMI,NMR)          ,pdW0nu(3,NMI,NMR)

      INTEGER  IMI1(1),IMI2(1),IMR1(1),IMR2(1)
      REAL*8   dtrans(NLAMBDA),dtrans0(NLAMBDA)   !t(outbound);t(inbound), interpolated for candidate mi,mr

      REAL     bandwt(NLAMBDA)

      REAL*8   f_lambda(NLAMBDA)

      INTEGER  jl
      REAL*8   xmi,ymr(1),C,acdm,bbp,xmi_4(1),
     &         zw0_iter(1)   , pdZw0(2),     !outputs for internal iterations only
     &         zv_iter(1)    , pdZv(2),      !outputs for internal iterations only
     &         zTau865_iter(1) , pdZtau865(2), 
     &         zTau(NLAMBDA) , 
     &         zRho(NLAMBDA) , pdZrho(2,NLAMBDA),
     &         zA(NLAMBDA)   , pdZa(2,NLAMBDA),
     &         zB(NLAMBDA)   , pdZb(2,NLAMBDA),
     &         zA0(NLAMBDA)  , pdZa0(2,NLAMBDA),
     &         zB0(NLAMBDA)  , pdZb0(2,NLAMBDA),
     &         zC865(NLAMBDA), pdZc865(2,NLAMBDA),
     &                         pdRmi(NLAMBDA),
     &                         pdRmr(NLAMBDA),
     &         zwRho(NLAMBDA),
     &                         pdRc(NLAMBDA),
     &                         pdRacdm(NLAMBDA),
     &                         pdRbbp(NLAMBDA) 

      REAL*8   LwN(NLAMBDA)

      INTEGER  ixMi,ixMr,ixC,ixAcdm,ixBbp
      PARAMETER (ixMi=1,ixMr=2,ixC=3,ixAcdm=4,ixBbp=5)
      REAL*8   bw2,rhoDiff(NLAMBDA),
     &         modelRho(NLAMBDA),pfmu0(NLAMBDA),pp,pp2,rconst,
     &         sconst,sconststd,stddev

c selection of master function band weight
      DATA bandwt /1.0, 1.0, 1.0, 1.0,
     &             1.0, 1.0,
     &             1.0, 1.0 /                       ! equal weights
     
      DATA     sconststd /100.D0/
      DATA     sconst /1.D+4/      !for ifcnlsq = 0

c      DATA     sconst /1.D+5/     !for ifcnlsq = 2
cc      DATA     sconst /1000.D0/

c SD master function control     
c      iFcnLSQ = 0    ! SQRT(WEIGHTED DIFFERENCE**2) Slsq
c      iFcnLSQ = 2   ! WEIGHTED DIFFERENCE**2 [sum_(dev_i**2)],  Slsq, used for Case-1 waters
      iFcnLSQ = 1   ! REL.STD.DEV Slsq [sqrt((1/N-1)*sum_(dev_i/giv_i-1))], std.dev (other tests)

c initial 'candidate' values re-initialised
      xmi  = x(ixMi)     !x(1)=mi
      ymr(1)  = x(ixMr)     !x(2)=mr
      C    = x(ixC)      !x(3)=C
      acdm = x(ixAcdm)   !x(4)=acdm
      bbp  = x(ixBbp)    !x(5)=bbp

c recall that mi is actually mi^4; reset variable name for clarity
      xmi_4(1) = xmi 
      
c conversion function radiance --> reflection (for LwN->Rho)
      DO jl=1,NBand
        pfmu0(jl) = DBLE(pi/(F0(jl)*mu0))
      ENDDO

C ------- LOCATE MI,MR GRID SQUARE (IN LARGER MI*MR GRID) THAT 'CANDIDATES' MI,MR LIE WITHIN -----
c index of (IMI,IMR)
c we have a mi,mr grid of resolution 6*2 (x-y plane)
c output given as interval coordinate of square (x-y)
c 
c DRGLCTN  IN: 6,2,mi_table_values(6),mr_table_values(2),#output points=1, candidate_mi(1),candidate_mr(1)
c         OUT: x_position INTERVAL NUMBER(1) of candidate mi, y_position INTERVAL NUMBER(1) of candidate mr
c The resultant x-y position of 'candidates' is given as an INTERVAL NUMBER (integer) in the mi,mr grid for 
c subsequent use in interpolation procedure in the z direction (call to subroutine DRGPLNL).  i.e, IMI,IMR describe
c the interval (bounded by adjacent mi,mr table values) that the candidate mi,mr lies within.
c x-plane = mi

      CALL DRGLCTN(NMI,NMR,dmi_4,dmr  ,1,xmi_4,ymr  ,IMI1,IMR2)
c x-plane = mr
      CALL DRGLCTN(NMR,NMI,dmr  ,dmi_4,1,ymr  ,xmi_4,IMR1,IMI2)

C ------- MASTER FUNCTION; GIVEN (MI,MR),(C,acdm,bbp) LOOKING FOR REFLECTANCE (RHOA+RHOW) --------
c         AT ALL WAVELENGTHES AND CALCULATING SLSQ-FUNCTION (FOR ALL WAVELENGTHES<=(NBand)

c     initialise SD and absolute percentage difference in master function
      stddev = 0.D0 

c     LwN - given initial/startup values of  C,acdm,bbp on first entrance into optimization LBFGSB (subroutine setulb)
c           or
c         - given 'current candidates' C,acdm,bbp in subsequent searches
      CALL DGS97 (C,acdm,bbp,LwN)


C -------INTERPOLATION OF TAU865 AS A FUNCTION OF 'CANDIDATE' MI,MR (xmi,ymr)

c The input 'candidates' mi,mr lie in one of the 5*1=5 smaller rectangular cells of the 
c large 6*2 grid of mi,mr table values (see declarations in subroutine ATMCOR_SOA).

c Example: "see call to DRGPLNL with comment XXXX at start"
c A 'candidate' tau865 is given at each mi,mr grid point(6*2) - optimized value derived in DRV_TAUVV_ZXMWD.
c A polynomial is fit through the candidate mi,mr point within the smaller rectangular cell.
c If candidate mi or mr (or both) lie outside the 6*2 table, the rectangular cell is considered
c to approach infinity in either the mi or mr directions (or both).
c tau865 is then interpolated to 'candidates' mi,mr by evaluating the polynomial for rectangular-
c grid bivariate interpolation and surface fitting.
c DRGPLNL IN: 6,2,mi_table_values(6),mr_table_values(2),
c             candidate (optimal) tau865(6,2),partial derivative of tau865(3,6,2) as a function of mi,mr
c                (note: for tau865(3,6,2), 3 = partial derivatives at zx,zy,zxy; 6,2=mi,mr of course),
c             #output points=1,candidate_mi(1),candidate_mr(1),
c             x_position INTERVAL NUMBER(1) of candidate mi, y_position INTERVAL NUMBER(1) of candidate mr
c        OUT: interpolated tau865(1),interpolated pdtau865(2,1)
c                [note: interpolated pdtau865 has two values, the interpolated partial derivatives of tau865
c                       at both candidate_mi(1) and candidate_mr(1)]

c #################################################

cccc  this call to DRGPLNL adds output W0 and nu at each iteration of the optimization
cccc  not required for processing       
                       
c                                                                    ----- w0_865(xmi,ymr)
      CALL DRGPLNL(NMI,NMR,dmi_4,dmr,dw0nu,pdW0nu,1,xmi_4,ymr, 
     &             IMI1,IMR2,zw0_iter,pdZw0)

      w0_iter     = REAL(zW0_iter(1))        !for iterative output only

C                                                                    ----- v(xmi,ymr)
      CALL DRGPLNL(NMI,NMR,dmi_4,dmr,dvir,pdVir,1,xmi_4,ymr,
     &             IMI1,IMR2,zv_iter,pdZv)

      v_iter     = REAL(zV_iter(1))         !for iterative output only

c ################################################

c x-axis (plane) as MI
c                                                                    ----- tau865(xmi,ymr)
      CALL DRGPLNL(NMI,NMR,dmi_4,dmr,  dtau865,pdTau865,1,xmi_4,ymr, 
     &             IMI1,IMR2,zTau865_iter,pdZtau865)

      tau865_iter = REAL(zTau865_iter(1)) !for iterative output only


C ----INTERPOLATION OF WAVELENGTH-DEPENDENT PARAMS; CALCULATION OF CANDIDATE WATER REFLECTANCE;
C     CALCULATION OF SD AND LSQ MASTER FUNCTION : measured(RHOT+RHORAY)~modelled(RHOA+RHOW);
c     ALL AS A FUNCTION OF 'CANDIDATE' MI,MR (xmi,ymr)
                        
c comments: same as for tau865 (above) but for respective parameters; all interpolated for 'candidate' mi,mr
                                                             
      DO jl=1,NBand


c all following have x-axis (plane) as MI
c       candidate rhoA [= candRho(xmi,ymr,jl)]; partial derivative rhoA
        CALL DRGPLNL(NMI  ,NMR, dmi_4, dmr ,
     &               dcandRho(1,1,jl),pdCandRho(1,1,1,jl),1,
     &               xmi_4,ymr, IMI1 , IMR2,zRho(jl),pdZrho(1,jl))

c       candidate watRho(jl); including reflection conversion, Case 1 mode
        zwRho(jl) = pfmu0(jl)*LwN(jl)    ! nLw >> wRho

c       candidate tau(lambda)/tau(865) [= cc865(xmi,ymr,jl)]; partial derivative tau(lambda)/tau(865)
        CALL DRGPLNL(NMI,NMR  ,dmi_4,dmr,
     &               dcc865(1,1,jl), pdcc865(1,1,1,jl),1,
     &               xmi_4,ymr,IMI1,IMR2, zC865(jl),pdZc865(1,jl))

c all following have x-axis (plane) as MR 
c       diffuse transmittance A coeff at theta [zA(xmi,ymr,jl)] and partial derivative       
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffA_ri(1,1,jl), pdA(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zA(jl),pdZa(1,jl))

c       diffuse transmittance B coeff at theta [zB(xmi,ymr,jl)] and partial derivative
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffB_ri(1,1,jl), pdB(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zB(jl),pdZb(1,jl))

c       diffuse transmittance A coeff at theta0 [zA0(xmi,ymr,jl)] and partial derivative
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffA0_ri(1,1,jl), pdA0(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zA0(jl),pdZa0(1,jl))

c       diffuse transmittance A coeff at theta0 [zA0(xmi,ymr,jl)] and partial derivative                                             
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffB0_ri(1,1,jl), pdB0(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zB0(jl),pdZb0(1,jl))


c       optical depth at each wavelength, interpolated for 'candidate' mi,mr
        zTau(jl)    = zTau865_iter(1)*zC865(jl)

c       diffuse transmittance (outbound;inbound) at each wavelength, interpolated for 'candidate' mi,mr
c       t*(theta)~Acoeff*exp[-B(theta)*tau]
        dtrans(jl)  = zA(jl) *DEXP(zB(jl)*zTau(jl))
        dtrans0(jl) = zA0(jl)*DEXP(zB0(jl)*zTau(jl))

c ####################################################################
c for iterative output only

        tau(jl)     = REAL(zTau(jl))
        trans(jl)   = REAL(dtrans(jl))
        trans0(jl)  = REAL(dtrans0(jl))

        resAerRho(jl) = REAL(zRho(jl))
        resWatRho(jl) = REAL(zwRho(jl))
        resTotRho(jl) = resAerRho(jl) +
     &           REAL(dtrans(jl))* resWatRho(jl)*REAL(dtrans0(jl))

c ####################################################################

c       MODELLED(rhoA+rhow) = rhoA+ t*t0*rhow, interpolated for 'candidate' mi,mr
c       rhoA(6) is band averaged using gordon(1995) 
        modelRho(jl)    = zRho(jl) + dtrans0(jl)*zwRho(jl)*dtrans(jl)

c       reflection difference = MODELLED(rhoA+rhow) - MEASURED(rhot-rhoray)
        rhoDiff(jl) = modelRho(jl)-DBLE(aerwatRho(jl))

c calculation of differences; spectral weighting incorporated (bandwt)
c stddev 'accumulates' with each wavelength, dependent on the modelRho~real comparison
c (see statement above) at each wavelength, the spectral weighting of the wavelength and the type of the 
c 'difference' function (ifcnlsq=0,1 or 2)
        IF (iFcnLSQ.EQ.1) THEN
           pp  = DBLE(bandwt(jl))*
     &                  (modelRho(jl)/DBLE(aerwatRho(jl))-1.D0)
           pp2 = pp*pp         
        ELSE IF (iFcnLSQ.EQ.0.OR.iFcnLSQ.EQ.2) THEN
           pp  = DBLE(bandwt(jl))*(modelRho(jl)-DBLE(aerwatRho(jl)))
           pp2 = pp*pp
        ENDIF
        stddev = stddev + pp2
      
        f_lambda(jl) = 
     &    ((modelRho(jl)/DBLE(aerwatRho(jl))-1.D0))*sconststd  !*100

      ENDDO  ! jl

c LSQ master function; final calculation after all wavelengths considered
      IF (iFcnLSQ.EQ.1) THEN
        f       = sconststd * DSQRT(stddev/DBLE(NBand-1))      !*100
      ELSE IF (iFcnLSQ.EQ.0) THEN                              !*10,000 
        f       = sconst * DSQRT(stddev)
      ELSE IF (iFcnLSQ.EQ.2) THEN                              !*100,000
        f       = sconst * stddev
      ENDIF
        
C ------------------- CALCULATE GRADIENTS FOR EACH MI,MR,C,ACDM,BBP ------------------------
C                            gradiant(param) = D_Flsq/D_param
C -------------------
C *** JACOBIAN ***
C -------------------
c ------------------------------------------------------------------------------------------
c D_[rhoA+rhow] / D_(mi,mr,C,acdm,bbp)

      CALL DPDR_MIMR   ( NBand, zwRho, pdZrho,
     &                   zTau865_iter(1),pdZtau865,zTau,
     &                   zC865,pdZc865,
     &                   dtrans0,zA0,zB0,pdZa0,pdZb0,
     &                   dtrans ,zA ,zB ,pdZa ,pdZb ,
     &                   pdRmi,pdRmr                 )
c output:
c  pdRmi = interpolated partial derivative rhoA+rhow (as a function of candidate (optimal) mi)
c  pdRmr = interpolated partial derivative rhoA+rhow (as a function of candidate (optimal) mr)   

c ------------------------------------------------------------------------------------------
c D_LwN / D_(C,acdm,bbp) (note: this is not reflectance yet)

      CALL DPDLWN_GS97 ( NBand,
     &                   C,acdm,bbp,
     &                   pdRc,pdRacdm,pdRbbp )
c output:
c  pdLc    = partial derivative LwN as a function of C
c  pdLacdm =                  "                      acdm
c  pdLbbp  =                  "                      bbp

c      mi4_34 = 0.25D0/xmi**0.75D0    ! (reversal ^4) w.r.t. mi if required

c conversion function: surface radiance --> TOA reflection (for LwN->Rho)
      DO jl=1,NBand
        rconst    = pfmu0(jl)*dtrans0(jl)*dtrans(jl)  !surface radiance --> TOA reflectance
c        pdRmi  (jl) = mi4_34*pdRmi(jl)       ! w.r.t. mi if required
        pdRmi  (jl) = pdRmi(jl)                 ! w.r.t mi4
        pdRc   (jl) = rconst*pdRc(jl)         !TOA p.d. rhow as a function of C
        pdRacdm(jl) = rconst*pdRacdm(jl)      !TOA p.d. rhow as a function of acdm
        pdRbbp (jl) = rconst*pdRbbp(jl)       !TOA p.d. rhow as a function of bbp
      ENDDO
c ------------------------------------------------------------------------------------------


C g(p) = D_Flsq/D_p = SUM(D_Flsq(jl)/D_p) ***

      g(ixMi)  =0.D0             !g(1)=0.0 for mi
      g(ixMr)  =0.D0             !g(2)=0.0 for mr
      g(ixC)   =0.D0             !g(3)=0.0 for C
      g(ixAcdm)=0.D0             !g(4)=0.0 for acdm
      g(ixBbp) =0.D0             !g(5)=0.0 for bbp

c note: rhoDiff = MODELLED(rhoA+rhow) - MEASURED(rhot-rhoray)
c of course,  g(param) 'accumulates' with each wavelength, dependent on the partial derivative (param)
c at each wavelength, the spectral weighting of the wavelength and the type of the lsq 
c function (ifcnlsq=0,1 or 2)
      IF (iFcnLSQ.EQ.0) THEN

        DO jl=1,NBand
          bw2       = DBLE(bandwt(jl))*DBLE(bandwt(jl)) !(spectral weighting)^2
          rconst    = sconst*sconst*bw2*rhoDiff(jl)/f   ![MOD(rhoA+rhow)-MEAS(rhot-rhoray)/LSQ]
          g(ixMi)   = g(ixMi)   + rconst*pdRmi(jl)     
          g(ixMr)   = g(ixMr)   + rconst*pdRmr(jl)      !g(mr)=g(mr)+[MOD(rhoA+rhow)-MEAS(rhot-rhoray)/LSQ]
          g(ixC)    = g(ixC)    + rconst*pdRc(jl)       !                       *p.d.rhoA+rhow (as a function of candidate (optimal) mr)
          g(ixAcdm) = g(ixAcdm) + rconst*pdRacdm(jl)
          g(ixBbp)  = g(ixBbp)  + rconst*pdRbbp(jl)
        ENDDO

      ELSE IF (iFcnLSQ.EQ.1) THEN

        DO jl=1,NBand
          bw2       = DBLE(bandwt(jl))*DBLE(bandwt(jl))            !(spectral weighting)^2
          rconst    = sconststd*bw2*DSQRT(1.D0/DBLE(NBand-1))*     !sqrt(1/5)*100*(weight)^2* [MOD(rhoA+rhow)-MEAS(rhot-rhoray)]/ MEAS(rhot-rhoray)]^2
     &                rhoDiff(jl)/DBLE(aerwatRho(jl)*aerwatRho(jl))
          g(ixMi)   = g(ixMi)   + rconst*pdRmi(jl)
          g(ixMr)   = g(ixMr)   + rconst*pdRmr(jl)
          g(ixC)    = g(ixC)    + rconst*pdRc(jl)
          g(ixAcdm) = g(ixAcdm) + rconst*pdRacdm(jl)
          g(ixBbp)  = g(ixBbp)  + rconst*pdRbbp(jl)
        ENDDO

      ELSE IF (iFcnLSQ.EQ.2) THEN

        DO jl=1,NBand
          bw2       = DBLE(bandwt(jl))*DBLE(bandwt(jl)) !(spectral weighting)^2
          rconst    = 2.D0*sconst*bw2*rhoDiff(jl)       !2*[MOD(rhoA+rhow)-MEAS(rhot-rhoray)] 
          g(ixMi)   = g(ixMi)   + rconst*pdRmi(jl)
          g(ixMr)   = g(ixMr)   + rconst*pdRmr(jl)
          g(ixC)    = g(ixC)    + rconst*pdRc(jl)
          g(ixAcdm) = g(ixAcdm) + rconst*pdRacdm(jl)
          g(ixBbp)  = g(ixBbp)  + rconst*pdRbbp(jl)
        ENDDO

      ENDIF

      DFCNLSQ_GRAD = .TRUE.                !function is computable at this point
      RETURN

      END
     
C *****************************************************************************

C Post-processing of optimised data for final output in SUBROUTINE ATMCOR_SOA.F
C A repetition of the mi,mr  grid setup and interpolations undertaken in FUNCTION DFCNLSQ_GRAD,
C but using the final optimized values of mi,mr from the optimization and the extra inclusion of 
C final mi.mr-interpolated W0,nu and tau865.
C Also includes final calculations of t*,tau(lambda),rhoW,rhoA and (rhoA+rhow)

      SUBROUTINE SPECOPT_POSTPROC (N,x,aerwatRho,
     &                             dmi_4,dmr,
     &                             dcandRho,pdCandRho,
     &                             dtau865,pdTau865,dcc865,pdcc865,
     &                             dw0nu,pdW0nu,dvir,pdVir,
     &                             ddiffA_ri ,pdA ,ddiffB_ri ,pdB ,
     &                             ddiffA0_ri,pdA0,ddiffB0_ri,pdB0,
     &                                                   F0,mu0,pi,
     &                             w0,v,
     &                             trans0,trans,
     &                             tau865,tau,
     &                             resAerRho,resWatRho,resTotRho,
     &                             absCoef  )

C BEGIN$DECLARE

      IMPLICIT NONE

      INTEGER  N
      REAL*8   x(N)

      INTEGER  NMI,NMR,NLAMBDA
      PARAMETER (NMI=6,NMR=2,NLAMBDA=8)

      INTEGER    NMAX
      PARAMETER (NMAX=5)

      REAL     aerwatRho(NLAMBDA),F0(NLAMBDA),mu0,pi
      REAL*8   dmi_4(NMI),dmr(NMR)
      REAL*8   dcandRho(NMI,NMR,NLAMBDA),pdCandRho(3,NMR,NMI,NLAMBDA),
     &         ddiffA_ri(NMR,NMI,NLAMBDA) ,pdA(3,NMR,NMI,NLAMBDA) ,
     &         ddiffB_ri(NMR,NMI,NLAMBDA) ,pdB(3,NMR,NMI,NLAMBDA) ,
     &         ddiffA0_ri(NMR,NMI,NLAMBDA),pdA0(3,NMR,NMI,NLAMBDA),
     &         ddiffB0_ri(NMR,NMI,NLAMBDA),pdB0(3,NMR,NMI,NLAMBDA),
     &         dcc865(NMI,NMR,NLAMBDA) ,pdcc865(3,NMI,NMR,NLAMBDA),
     &         dtau865(NMI,NMR)        ,pdTau865(3,NMI,NMR),
     &         dvir(NMI,NMR)           ,pdVir(3,NMI,NMR),
     &         dw0nu(NMI,NMR)          ,pdW0nu(3,NMI,NMR)

      REAL     w0,v,
     &         trans(NLAMBDA),trans0(NLAMBDA),tau865,tau(NLAMBDA)
      REAL*8   dtrans(NLAMBDA),dtrans0(NLAMBDA)
      REAL     resAerRho(NLAMBDA),resWatRho(NLAMBDA),
     &         resTotRho(NLAMBDA)

      REAL*8   absCoefNw(NLAMBDA)
      REAL     absCoef(NLAMBDA)

      INTEGER  IMI1(1),IMI2(1),IMR1(1),IMR2(1)
      INTEGER  jl

      REAL*8   xmi,ymr(1),C,acdm,bbp,xmi_4(1),
     &         zRho(NLAMBDA) , pdZrho(2,NLAMBDA),
     &         zA(NLAMBDA)   , pdZa(2,NLAMBDA),
     &         zB(NLAMBDA)   , pdZb(2,NLAMBDA),
     &         zA0(NLAMBDA)  , pdZa0(2,NLAMBDA),
     &         zB0(NLAMBDA)  , pdZb0(2,NLAMBDA),
     &         zTau(NLAMBDA) , 
     &         zTau865(1)       , pdZtau865(2),
     &         zW0(1)           , pdZw0(2),
     &         zV(1)            , pdZv(2),
     &         zC865(NLAMBDA), pdZc865(2,NLAMBDA),
     &         zwRho(NLAMBDA)

      REAL*8   LwN(NLAMBDA),pfmu0(NLAMBDA)

      INTEGER    ixMi  ,ixMr  ,ixC  ,ixAcdm  ,ixBbp
      PARAMETER (ixMi=1,ixMr=2,ixC=3,ixAcdm=4,ixBbp=5)

C END$DECLARE

      xmi  = x(ixMi)
      ymr(1)  = x(ixMr)
      C    = x(ixC)
      acdm = x(ixAcdm)
      bbp  = x(ixBbp)
C
crc      xmi_4 = xmi**0.25D0       ! w.r.t. mi
      xmi_4(1) = xmi                  ! w.r.t. mi4

      DO jl=1,NLAMBDA
        pfmu0(jl) = DBLE(pi/(F0(jl)*mu0))
      ENDDO
C                                                               ----- idx (IMI,IMR)
      CALL DRGLCTN(NMI,NMR,dmi_4,dmr  ,1,xmi_4,ymr  ,IMI1,IMR2)
      CALL DRGLCTN(NMR,NMI,dmr  ,dmi_4,1,ymr  ,xmi_4,IMR1,IMI2)
C                                                               ----- LwN
      CALL DGS97_CASE2 (C,acdm,bbp,LwN,absCoefNw)
C                                                               ----- tau865(xmi,ymr)
      CALL DRGPLNL(NMI,NMR,dmi_4,dmr,dtau865,pdTau865,1,xmi_4,ymr,
     &             IMI1,IMR2, zTau865,pdZtau865)
C                                                               ----- w0_865(xmi,ymr)
      CALL DRGPLNL(NMI,NMR,dmi_4,dmr,dw0nu,pdW0nu,1,xmi_4,ymr,
     &             IMI1,IMR2, zW0,pdZw0)
C                                                               ----- v(xmi,ymr)
      CALL DRGPLNL(NMI,NMR,dmi_4,dmr,dvir,pdVir,1,xmi_4,ymr,
     &             IMI1,IMR2,zV,pdZv)

      tau865 = REAL(zTau865(1))
      w0     = REAL(zW0(1))
      v      = REAL(zV(1))
      
C *** WAVELENGTH-DEPENDENT STUFF ***

      DO jl=1,NLAMBDA
C                                                               ----- candRho(xmi,ymr,jl)
        CALL DRGPLNL(NMI  ,NMR, dmi_4, dmr ,
     &               dcandRho(1,1,jl),pdCandRho(1,1,1,jl),1,
     &               xmi_4,ymr, IMI1 , IMR2,zRho(jl),pdZrho(1,jl))
C                                                               ----- watRho(jl)
        zwRho(jl) = pfmu0(jl)*LwN(jl)    ! nLw >> wRho

C absorption coefficient non-water parameters
        absCoef(jl) = REAL(absCoefNw(jl))
C                                                               ----- cc865(xmi,ymr,jl)
        CALL DRGPLNL(NMI,NMR  ,dmi_4,dmr,
     &               dcc865(1,1,jl), pdcc865(1,1,1,jl),1,
     &               xmi_4,ymr,IMI1,IMR2, zC865(jl),pdZc865(1,jl))
C                                                               ----- zA(xmi,ymr,jl)
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffA_ri(1,1,jl), pdA(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zA(jl),pdZa(1,jl))
C                                                               ----- zB(xmi,ymr,jl)
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffB_ri(1,1,jl), pdB(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zB(jl),pdZb(1,jl))
C                                                               ----- zA0(xmi,ymr,jl)
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffA0_ri(1,1,jl), pdA0(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zA0(jl),pdZa0(1,jl))
C                                                               ----- zB0(xmi,ymr,jl)
        CALL DRGPLNL(NMR,NMI  ,dmr,dmi_4,
     &               ddiffB0_ri(1,1,jl), pdB0(1,1,1,jl),1,
     &               ymr,xmi_4,IMR1,IMI2, zB0(jl),pdZb0(1,jl))


        zTau(jl)    = zTau865(1)*zC865(jl)
        dtrans(jl)  = zA(jl) *DEXP(zB(jl)*zTau(jl))
        dtrans0(jl) = zA0(jl)*DEXP(zB0(jl)*zTau(jl))
        tau(jl)     = REAL(zTau(jl))
        trans(jl)   = REAL(dtrans(jl))
        trans0(jl)  = REAL(dtrans0(jl))

        resAerRho(jl) = REAL(zRho(jl))
c        resWatRho(jl) = REAL(dtrans(jl)*zwRho(jl)*dtrans0(jl))
        resWatRho(jl) = REAL(zwRho(jl))
        resTotRho(jl) = resAerRho(jl) +
     &           REAL(dtrans(jl))* resWatRho(jl)*REAL(dtrans0(jl))

      ENDDO  ! jl

      RETURN
      END

C ############################################################################
C ############################################################################

C SECTION 3:
                
C Data set up and interpolation routines required by :
C SPEC1GEOFIT (main controlling routine) and
C DRV_TAUVV_ZXMWD driving routine (NIR optimization of nu, tau, mi, mr using ZXMWD)

C ****************************************************************************

c     Locate the minimum index-element in an array such that
c     XTAB(i) :: X :: XTAB(i+2) (or outside, if out of the range).
c     To be used in conjunction with the QUADRATIC Interpolation Routines.
c     Note: Xtab >= 3

      SUBROUTINE BRACKET4QUAD (NX,Xtab,X,jMin)

c input:
c        NX,Xtab,X
c output:
c        jmin

      IMPLICIT NONE
      INTEGER  NX,jMin
      REAL     X,Xtab(NX)

      jMin=1
      IF(X.GE.Xtab(NX)) THEN
         jMin=NX-1               !if X > max boundary
      ELSE IF(X.LE.Xtab(1)) THEN
         jMin=1                  !if X < min boundary
      ELSE
        jMin=0
        DO WHILE (X.GT.Xtab(jMin+1))     !main loop
           jMin  = jMin+1
        ENDDO
      ENDIF

      RETURN
      END

C ****************************************************************************

c converts Haze-C model given by its parameters into a unique model index number
c i.e. (jv,jr,ji) -> (jmodel).  Inverse to IDX2MODEL subroutine

      SUBROUTINE MODEL2IDX(jv,jr,ji,NV,NMR,NMI,jmodel)

c     jmodel = ji + NMI*(jr-1) + NMR*NMI*(jv-1)
c     where the first changing index is 'ji' -> then 'jr' -> then 'jv'

c input:
c        jv,jr,ji    -loop indices for nu,mr,mi
c        NV,NMR,NMI  -number of total elements of nu,mr,mi 
c output:
c        jmodel      -unique model index number

      IMPLICIT NONE
      INTEGER*4 jv,jr,ji,NV,NMR,NMI,jmodel

      jmodel = ji + NMI*(jr-1) + NMR*NMI*(jv-1)

      RETURN
      END

C ****************************************************************************

c converts Haze-C model INDEX given by (jmodel) into parameter
c indices i.e. (jmodel)->(jv,jr,ji). Inverse to MODEL2IDX subroutine

      SUBROUTINE IDX2MODEL(jmodel,NV,NMR,NMI,jv,jr,ji)

c     (jv-1)= IFIX(jmodel/NMR/NMI)
c     (jr-1)= IFIX((jmodel-(jv-1)*NMR*NMI)/NMI)
c     (ji)   =jmodel-(NMI*(jr-1) + NMR*NMI*(jv-1))

c input:
c        jmodel         -unique model index number
c        NV,NMR,NMI     -number of total elements of nu,mr,mi 
c output:
c        jv,jr,ji       -loop indices for nu,mr,mi

      IMPLICIT NONE
      INTEGER*4 jv,jr,ji,NV,NMR,NMI,jmodel
      INTEGER*4 jj

      jj = MOD(jmodel,NMR*NMI)
      IF(jj.EQ.0) jj=NMR*NMI
      jv = 1 + (jmodel-jj)/(NMR*NMI)   ! fixing (jv)
      ji = MOD(jj,NMI)                 ! fixing (ji)
      IF(ji.EQ.0) ji=NMI
      jr = 1 + (jj-ji)/NMI             ! fixing (jr)

      RETURN
      END

C ****************************************************************************

      SUBROUTINE QUADINTRP(doLog,x,x0,h,ydat,y)

c input: dolog.x.x0,h,ydat
c output: y

C Purpose:      given THREE data-pairs (x0,y0),(x1,y1),(x2,y2) find coeff's
C               of approximating Quadratic Polynomial (Parabola) that
C               passes through given points (the best usage with 
C               LOGDEV_DRIVER, i.e. like LGDRV_QUADINTRP).

C Arguments:
C               h        - grid-step
C               a        - a-coeff for parabola
C               b        - b-coeff
C               c        - c-coeff
C               x0       - supplied starting x_value (x0,x0+h,x0+2*h)
C               ydat     - supplied y_data(0:2)
C               ylog     - temporal for log-interpolation regime
C               doLog    - logical whether to do log-interpolation
C               lNegFlag - check, if a negative y(i) was encounted => shift
C               yneg     - the least negative value
C               ymove    - all y(i)'s are shifted (ymove) into positive values
C               x        - x-value of interest
C               y        - y=f(x) (y = exp(flog(x))-yshift in LOG-interpolation)

C Algorithm:    There are two regimes: REGULAR and LOG fit. If (doLog) is set,
C               checking if all y(i)'s are positive to assure log(y) is OK.
C               If not, then making a shift of all y's into positives. Then,
C               in BOTH cases:
C               Trying to fit parabola through given points (this procedure
C               is equivalent to Lagrange Interpolation but more speedy since
C               points are taken to be equally spaced such that direct
C               analytical formulas are found via system of equations y1,y2,y3
C                      y2 -->                               o
C        
C                      y1 -->                   o
C                      y0 -->       o
C                             ----- * --------- * --------- * -------
C                                   x0       x1=x0+h     x2=x0+2*h
C
C               Parabola:    y(x) = a*x^2 + b*x + c ,  where
C
C        y2-2*y1+y0        y1-y0
C    a = ---------- ;  b = ----- - 2*a*x0 - a*h ;  c = y0 - a*x0^2 - b*x0 ;
C          2*h^2             h
C
C Comments:     the accuracy as f(x) = f(x0) + f'(x0)*(x-x0) + f"(x0)*(x-x0)^2/2
C               Taking into account that f'(x),f"(x) are not necessarily known,
C               its been shown that errors in calculating deriv's along with the
C               series expansion provide inadequate accuracy and speed.

      IMPLICIT NONE
      REAL     rMaxExp
      PARAMETER (rMaxExp=37.0)
      INTEGER  i
      REAL     a,b,c,abc
      REAL     x,x0,h,ydat(0:2),ylog(0:2),y
      REAL     yneg,ymove
      LOGICAL  doLog,lNegFlag,doLogLocal

      doLogLocal = doLog

10    CONTINUE                     ! GOTO here if Exp crashes

      IF(.NOT.doLogLocal) THEN
C ----------- FITTING THE PARABOLA AND FINDING INTERPOLATED POINT
        a = (ydat(2) - 2.0*ydat(1) + ydat(0))/(2.0*h**2)
        b = (ydat(1) - ydat(0))/h - 2.0*a*x0 - a*h
        c = ydat(0) - a*x0**2 - b*x0
        y = a*x**2 + b*x + c
      ELSE
C ----------- SEARCHING FOR THE MOST NEGATIVE VALUE
        lNegFlag = .FALSE.
        yneg = 0.0
        DO i=0,2
          IF(ydat(i).LE.0.0) THEN
             lNegFlag = .TRUE.
             IF(ydat(i).LT.yneg) yneg = ydat(i)
          ENDIF
        ENDDO
C ----------- MOVING ALL Y's A DISTANCE (1.0-YNEG) TO ASSURE CORRECT LOG'S
        IF(lNegFlag) THEN
           ymove = 1.0 - yneg
        ELSE
           ymove = 0.0
        ENDIF
        DO i=0,2
          ylog(i) = ALOG(ymove+ydat(i))
        ENDDO
C ----------- FITTING PARABOLA TO POINTS
        a = (ylog(2) - 2.0*ylog(1) + ylog(0))/(2.0*h**2)
        b = (ylog(1) - ylog(0))/h - 2.0*a*x0 - a*h
        c = ylog(0) - a*x0**2 - b*x0
C ----------- CALCULATING THE INTERPOLATED POINT WITH A SHIFT
        abc = a*x**2 + b*x + c
        IF (ABS(abc).LE.rMaxExp) THEN
          y = EXP(a*x**2 + b*x + c) - ymove
        ELSE
crc          PRINT*, '*SpecOptim-Warning* QUADINTRP: ExpReal > rMaxExp'
crc          PRINT*, '     *Switching to Linear Interpolation*        '
          doLogLocal = .FALSE.
          GOTO 10
        ENDIF
      ENDIF

      RETURN
      END

C *****************************************************************************

      SUBROUTINE LAGRINTRP(N,xdat,ydat,a)

c input: N,xdat,ydat
c output: a

C Purpose:      given data-pairs (x0,y0),(x1,y1),...,(xND,yND) find coeff's
C               of approximating Lagrange Polynomial of ND's order using
C               NEWTON DIVIDED DIFFERENCES.
C
C Arguments:    N   - the order of approximating Lagrange Polynomial
C               a   - array of coeffs LagrPoly(x)=A(0)+A(1)*x^1+...+A(N)*x^N
C               c   - tmp_array of coeffs LagrPoly(x)=C(0)+C(1)*x^1+...+C(N)*x^N
C               dd  - array of NEWTON DIVIDED DIFFERENCES, such that
C                       the storage is as large as
C                         NDD = (N+1)*D0+N*D1+...+1*DN = (N+1)*(N+2)/2, where
C                       D0 - the #of zero-order divided diff's,
C                       D1 - the #of first-order divided diff's,...,
C                       DN - the #of N's-order divided diff's
C               ND   - (internal) current m's-order calculation
C               xdat - supplied x_data(0:N)
C               ydat - supplied y_data(0:N)
C               NLO  - (internal) the Number_Location_Old(zero_element)
C                                 of already calculated (m-1)'s-DD-order 
C               NLN  - (internal) the Number_Location_New in DD() assuming
C                                 current is the (m)'s-order calculation of DD

C Algorithm:    Newton Divided Differences:
C
C    [x0]=DD(0,0)
C                   [x0,x1]=DD(0,1)
C    [x1]=DD(1,0)                    [x0,x1,x2]=DD(0,2)
C                   [x1,x2]=DD(1,1)                     [x0,x1,x2,x3]=DD(0,3)
C    [x2]=DD(2,0)                    [x1,x2,x3]=DD(1,2)
C                   [x2,x3]=DD(2,1)
C    [x3]=DD(3,0)
C
C               where I introduced for convinience in algorithmic formulation
C               the above notation DD(i,m)= i's DividedDifference of m's order.
C               Next,
C                      [x0,x1]    = (y0-y1)/(x0-x1)
C                      [x0,x1,x2] = ([x0,x1]-[x1,x2])/(x0-x2)
C                                ...
C    Then the Lagrange Polynomial is:
C
C               f(x) =        [x0]          +
C                     pi0    *[x0,x1]       +
C                     pi1    *[x0,x1,x2]    +
C                         ...
C                     pi(N-1)*[x0,x1,...,xN]
C               where
C                     piK = (x-x0)(x-x1)(x-x2)...(x-xK) and is calculated
C               using the subroutine PILAGRANGE.
C
C      Note:    Since DD is stored in a single array the location of a specific
C               DD(i,m) can be found through (after fooling around with numbers
C               for a while):
C
C                 LOCATION_DD(i,m)_in_DD = m*(N+1)-m*(m-1)/2 + i
C
C -----------------------------------------------------------------------------
C Comments:     the #of data_points should be (N+1), enumeration starts from 0 !
C -----------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER  N,NMAX
      PARAMETER (NMAX=4)                     ! for this particular application
      INTEGER  ND,NLO,NLN,i,id,ia
      REAL     a(0:N),c(0:NMAX),dd(0:(NMAX+1)*(NMAX+2)/2 - 1)
      REAL     xdat(0:N),ydat(0:N)

c      IF(N.LT.0) PRINT*,'LAGRINTRP: supplied N < 0'   ! as a consequence above
c      IF(N.EQ.0) THEN
c        dd(0) = ydat(0)
c        RETURN
c      ENDIF

C ------------ LOOKING FOR NEWTON DIVIDED DIFFERENCES
      DO i=0,N
        dd(i) = ydat(i)
      ENDDO
      ND  = N
      NLO = 0
      DO id=1,N
        NLN = id*(N+1)-id*(id-1)/2
        ND  = ND-1
        DO i=0,ND
          dd(NLN+i)=(dd(NLO+i)-dd(NLO+1+i))/(xdat(i)-xdat(i+id))
        ENDDO
        NLO=NLN
      ENDDO
C ------------ INITIALIZING ARRAYS A(i) and C(i)
      DO i=0,N
        a(i) = 0.0
        c(i) = 0.0
      ENDDO
C ------------ CALCULATING THE RESULTING POLYNOMIAL COEFFICIENTS
      a(0) = dd(0)
      IF(N.GE.1) THEN
        c(0) = -xdat(0)
        c(1) = 1.
        a(0) = a(0)+c(0)*dd(N+1)
        a(1) =      c(1)*dd(N+1)
C ------------ CALCULATING PIn-functions and adding the result to desired coeffs
        DO i=2,N
c          CALL PILAGRANGE(i,xdat(i-1),c)            ! it sucks, really
          CALL PILAGRANGE(i,N,xdat(i-1),c)
          DO ia=0,i
            a(ia) = a(ia) + c(ia)*dd(i*(N+1)-i*(i-1)/2)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END

C *****************************************************************************

      SUBROUTINE PILAGRANGE(ND,NSTRG,x0,a)

c input: ND,NSTRG
c output: x0,a

C Purpose:      an auxiliary subroutine to calculate PIn(n+1) of degree (n+1)
C               given coefficients A[n],A[n-1],A[n-2],...,A[n-n]=A[0] of
C               the polynomial
C                     PIn(n) = A[n]*X^n + A[n-1]*X^(n-1) +...+ A[0]*1
C               of degree (n) and a constant value (X0) such that
C                     PIn(n+1) = (X-X0)*PIn(n)
C -----------------------------------------------------------------------------
C Comments:     the INPUT A[array]-values are changed at OUTPUT
C -----------------------------------------------------------------------------
C Arguments:    ND       -   the degree of PI(n+1)
C               NSTRG    -   required storage (to avoid conflicts with the call)
C               x0       -   the value of X0 in (X-X0)
C               a        -   array of coeff's in old Polynomial a[n]*X^n+...
C               atmp     -   temporal storage for coefficients

      IMPLICIT NONE
      INTEGER  ND,NSTRG,NSTRGMAX,NDMAX,i
      PARAMETER (NDMAX=7)
      PARAMETER (NSTRGMAX=4)                    ! for portability
      REAL     x0
      REAL     a(0:NSTRG),atmp(0:NSTRGMAX)

      IF(ND.GT.NDMAX)
     & PRINT*,'The degree of required Poly in PILAGRANGE exceeds 7'
      DO i=0,ND-1
        atmp(i) = -x0*a(i)
      ENDDO
      DO i=ND,1,-1
        a(i) = a(i-1)
      ENDDO
      a(0) = 0.0
      DO i=0,ND-1
        a(i) = a(i) + atmp(i)
      ENDDO

      RETURN
      END

C *****************************************************************************

        SUBROUTINE FUNCT_A_B_C_D_FOURIER_INTERPOLATE
     &                               (iphase, sun, theta, phi, 
     &                                acost,bcost,ccost,dcost)

c input:  
c        iphase        : equivalent to jmodel = LUT index value for 
c                        each of 72 models (mi,mr,nu)
c        sun,theta,phi : sun/sensor geometry; sun=theta0=SZA
c output: 
c        acost,bcost,ccost,dcost : coefficients, array size = 6 = nu step values

C Interpolate a,b,c,d fourier coefficients (for geometry) from tabulated ones.

c LUT values of the rhoA v Tau fourier coefficients a,b,c,d 
c have been computed as a quartic fit for each of 72 models (mi,mr,v) 
c at each sun angle, lambda and phi (33*8*15*72 = 285 120 values).
c They are pre-calculated into external LUT's and read into
c arrays (mini LUT's) luttau_a/b/c/d by subroutine READ_MODEL_DATA.
c They are passed as COMMON VARS to this subroutine:
c rhoa (params) = a(params)tau + b(params)tau^2 + c(params)tau^3 + d(params)tau^4

c Gordon notes:
c It uses 2-D linear interpolation for sun and viewing (polar) angles which are not 
c covered in the tables.  If theta < smallest angle in the table, then a, b and c 
c are computed for the smallest angle and averaged over the three azimuth angles
c phi = 0., 90. and 180.  If the viewing direction is alligned to viewing solar 
c photons specularly reflected from the sea surface, for these phi=0.  The 'phs_name'
c must have type as '_mh70' for the Maritime with RH = 70%.
c Note: the azimuth of the solar beam is taken to be zero.

c It assumes that the data are tabulated for theta0 = 0(2.5)80 deg.

c Input angles are in DEGREES

        INTEGER ndir,ndir_out,mphi,num,mxmodels,msun,nlambda
        PARAMETER (ndir=35, ndir_out=26, mphi=15, num=75)
        INTEGER nv,nmr,nmi,nmodel
        PARAMETER (mxmodels=1, msun=33, nlambda=8)
        PARAMETER (nv=6,nmr=2,nmi=6,nmodel=nv*nmr*nmi)

        REAL thetav(ndir)

        REAL c_coef_acost(ndir_out,mphi,2),c_coef_bcost(ndir_out,mphi,2)
        REAL c_coef_ccost(ndir_out,mphi,2),c_coef_dcost(ndir_out,mphi,2)

        REAL acost(nlambda), bcost(nlambda)
        REAL ccost(nlambda), dcost(nlambda)

        REAL mu0
        REAL as00, as01, as10, as11
        REAL ai00, ai01, ai10, ai11
        REAL ac00, ac01, ac10, ac11
        REAL ad00, ad01, ad10, ad11
        REAL factor, a_value, b_value, c_value, d_value
        REAL rad, sun, theta, rphi, x, sun_max
        REAL hview, hsun, p, q, pi, phi
        INTEGER iisun(2)
        INTEGER m, jv1, jv2, irem, ndang, isun, morder
        INTEGER k, j, nsun, i1, j1, jjsun, iphase

        REAL luttau_a(ndir,mphi,msun,nlambda,nmodel),
     &       luttau_b(ndir,mphi,msun,nlambda,nmodel),
     &       luttau_c(ndir,mphi,msun,nlambda,nmodel),
     &       luttau_d(ndir,mphi,msun,nlambda,nmodel)

        COMMON /comrtelut/ luttau_a,luttau_b,luttau_c,luttau_d

        DATA thetav / 1.021250, 2.344183, 3.674450, 5.005859,
     &                7.667694, 10.325176, 12.976083, 15.618372,
     &               18.250134, 20.869442, 23.474306, 26.062740,
     &               28.632664, 31.181959, 33.708408, 36.209732,
     &               38.683578, 41.127476, 43.538868, 45.915089,
     &               48.253349, 50.550777, 52.804333, 55.010883,
     &               57.167152, 59.269760, 61.315170, 63.299725,
     &               65.219673, 67.071136, 68.850151, 70.552650,
     &               72.174538, 73.711670, 75.159866 /

        rad(x) = x*pi/180.

        pi = 4.*atan(1.0)
        ndang = 25
        sun_max  = 2.5 * float(msun-1)

        if(sun .le. sun_max) then

        isun = 10*ifix(sun)
        rphi = rad(phi)

        if(sun .ne. sun_max) then
                irem = mod(isun,ndang)
                iisun(1) = isun - irem
                iisun(2) = iisun(1) + ndang
        else
                iisun(1) = 10*(sun_max-ndang)
                iisun(2) = 10*sun_max
        endif

c convert it the correct index for the tables.

        nsun = 1 + iisun(1)/ndang

        do 100 j = 1, nlambda            ! wavelength
          do 10  nsun = 1, 2               ! sun angle

c convert sun angle index to the correct index for the tables.
            jjsun = 1 + iisun(nsun)/ndang

            mu0 = cos(rad(float(iisun(nsun))/10.))
           if(mu0 .gt. 0.9999) then
                morder=1
           else
                morder=15               !maximum # you can choose is 15
           endif

           do i1 = 1, morder
             do j1 = 1, ndir_out
               c_coef_acost(j1,i1,nsun) =
     &                  luttau_a(j1,i1,jjsun,j,iphase)
               c_coef_bcost(j1,i1,nsun) = 
     &                  luttau_b(j1,i1,jjsun,j,iphase)
               c_coef_ccost(j1,i1,nsun) = 
     &                  luttau_c(j1,i1,jjsun,j,iphase)
               c_coef_dcost(j1,i1,nsun) = 
     &                  luttau_d(j1,i1,jjsun,j,iphase)
             enddo
           enddo

10         continue

c find the viewing angles for interpolation
  
           do 20 j1 = 1, ndir_out
             jv2 = j1
20           if(theta .lt. thetav(j1)) goto 21
21           jv1 = jv2 - 1

c theta is between thetav(jv1) and thetav(jv2)
        
c set up the variables required for interpolation

             if(jv1 .ne. 0) then
               hsun = float(iisun(2) - iisun(1))
               hview = thetav(jv2) - thetav(jv1)
               p = ( 10.*sun - float(iisun(1)) )/hsun
               q = ( theta - thetav(jv1) )/hview

c compute the acost, bcost and ccost

               acost(j) = 0.
               bcost(j) = 0.
               ccost(j) = 0.
               dcost(j) = 0.
        
               do m = 1, morder

                 if(m .eq. 1) then
                        factor = 1.0
                 else
                        factor = 2.0
                 endif

          if(cos(rad(float(iisun(1))/10.)) .gt. .9999     
     &                  .and. m .gt. 1) then   !if mu0=1. only one fourier
                                               !order is involved
                   as00 = 0.0
                   as10 = c_coef_acost(jv1,m,2)
                   as01 = 0.0
                   as11 = c_coef_acost(jv2,m,2)
        
                   ai00 = 0.0
                   ai10 = c_coef_bcost(jv1,m,2)
                   ai01 = 0.0
                   ai11 = c_coef_bcost(jv2,m,2)

                   ac00 = 0.0
                   ac10 = c_coef_ccost(jv1,m,2)
                   ac01 = 0.0
                   ac11 = c_coef_ccost(jv2,m,2)

                   ad00 = 0.0
                   ad10 = c_coef_dcost(jv1,m,2)
                   ad01 = 0.0
                   ad11 = c_coef_dcost(jv2,m,2)

                 else

                   as00 = c_coef_acost(jv1,m,1)   !see abromowitz and stegun
                   as10 = c_coef_acost(jv1,m,2)   !handbook of mathematical
                   as01 = c_coef_acost(jv2,m,1)   !functions. page 882
                   as11 = c_coef_acost(jv2,m,2)   !formula 25.2.66

                   ai00 = c_coef_bcost(jv1,m,1)
                   ai10 = c_coef_bcost(jv1,m,2)
                   ai01 = c_coef_bcost(jv2,m,1)
                   ai11 = c_coef_bcost(jv2,m,2)

                   ac00 = c_coef_ccost(jv1,m,1)
                   ac10 = c_coef_ccost(jv1,m,2)
                   ac01 = c_coef_ccost(jv2,m,1)
                   ac11 = c_coef_ccost(jv2,m,2)

                   ad00 = c_coef_dcost(jv1,m,1)
                   ad10 = c_coef_dcost(jv1,m,2)
                   ad01 = c_coef_dcost(jv2,m,1)
                   ad11 = c_coef_dcost(jv2,m,2)

                 endif

                 a_value = (1.-p)*(1.-q)*as00 + p*q*as11
     &                       + p*(1.-q)*as10 + q*(1.-p)*as01
                 b_value = (1.-p)*(1.-q)*ai00 + p*q*ai11
     &                       + p*(1.-q)*ai10 + q*(1.-p)*ai01
                 c_value = (1.-p)*(1.-q)*ac00 + p*q*ac11
     &                       + p*(1.-q)*ac10 + q*(1.-p)*ac01
                 d_value = (1.-p)*(1.-q)*ad00 + p*q*ad11
     &                       + p*(1.-q)*ad10 + q*(1.-p)*ad01

                 acost(j) = acost(j) +
     &                      factor*a_value*cos(rphi*float(m-1))
                 bcost(j) = bcost(j) + 
     &                      factor*b_value*cos(rphi*float(m-1))
                 ccost(j) = ccost(j) + 
     &                      factor*c_value*cos(rphi*float(m-1))
                 dcost(j) = dcost(j) + 
     &                      factor*d_value*cos(rphi*float(m-1))

               enddo

             else

               hsun = float( iisun(2) - iisun(1) )
               p = ( 10.*sun - float(iisun(1)) )/hsun

               acost(j) = 0.
               bcost(j) = 0.
               ccost(j) = 0.
               dcost(j) = 0.
        
               do k = 1, 3
               rphi = rad (float(k-1)*90.)
                 do m = 1, morder

                   if(m .eq. 1) then
                     factor = 1.0
                   else
                     factor = 2.0
                   endif

                   if(cos(rad(float(iisun(1))/10.)) .gt. .9999
     &                  .and. m .gt. 1) then

                     as00 = 0.0
                     as10 = c_coef_acost(1,m,2)

                     ai00 = 0.0
                     ai10 = c_coef_bcost(1,m,2)

                     ac00 = 0.0
                     ac10 = c_coef_ccost(1,m,2)

                     ad00 = 0.0
                     ad10 = c_coef_dcost(1,m,2)

                   else

                     as00 = c_coef_acost(1,m,1)
                     as10 = c_coef_acost(1,m,2)

                     ai00 = c_coef_bcost(1,m,1)
                     ai10 = c_coef_bcost(1,m,2)

                     ac00 = c_coef_ccost(1,m,1)
                     ac10 = c_coef_ccost(1,m,2)

                     ad00 = c_coef_dcost(1,m,1)
                     ad10 = c_coef_dcost(1,m,2)

                  endif

                  a_value = (1.-p)*as00 + p*as10
                  b_value = (1.-p)*ai00 + p*ai10
                  c_value = (1.-p)*ac00 + p*ac10
                  d_value = (1.-p)*ad00 + p*ad10

                  acost(j) = acost(j) +
     &                       factor*a_value*cos(rphi*float(m-1))
                  bcost(j) = bcost(j) + 
     &                       factor*b_value*cos(rphi*float(m-1))
                  ccost(j) = ccost(j) + 
     &                       factor*c_value*cos(rphi*float(m-1))
                  dcost(j) = dcost(j) + 
     &                       factor*d_value*cos(rphi*float(m-1))

                enddo   
              enddo

              acost(j) = acost(j)/3.
              bcost(j) = bcost(j)/3.          ! average value
              ccost(j) = ccost(j)/3.          ! average value
              dcost(j) = dcost(j)/3.          ! average value

            endif

100         continue

          else
        
          print*, 'the sun angle must be less or equal to', sun_max

        endif

        return
        end

C ############################################################################
C ############################################################################

C SECTION 4:

C Data set up and interpolation routines required by main optimization :
C DRV_SPECOPT_LBFGSB driving routine (main optimization)
C FUNCTION DFCNLSQ_GRAD (f and g updates for main optimization)
C SPECOPT_POSTPROC (postprocessing of main optimization output)

C *********************************************************************

C Required for interpolation of parameters as a function of optimization
C'candidaates' mi,mr.  Used in subroutine DRV_SPECOPT_LBFGSB.
C Large 6*2 mi,mr grid square has been reduced to a smaller rectangle whose
C vertices are the nearest mi,mr table values that bound the 'candidates' mi,mr from the
C optimization iteration (see SUBROUTINE DRGLCTN below).
C This subroutine determines a polynomial in x and y for a rectangle
C of the input grid in the x-y plane and calculates the z value for
C the desired points by evaluating the polynomial for rectangular-
C grid bivariate interpolation and surface fitting.             

      SUBROUTINE DRGPLNL(NXD,NYD,XD,YD,ZD,PDD,NIP,XI,YI,INXI,INYI,
     &                   ZI,PDZI)
*
* Polynomials for rectangular-grid bivariate interpolation and
* surface fitting
* (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)
*
* Added output of partial derivatives D_(x,y) at XI,YI in PDZI(2,NIP)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08

* The input arguments are
*   NXD  = number of the input-grid data points in the x
*          coordinate (must be 2 or greater),
*   NYD  = number of the input-grid data points in the y
*          coordinate (must be 2 or greater),
*   XD   = array of dimension NXD containing the x coordinates
*          of the input-grid data points (must be in a
*          monotonic increasing order),
*   YD   = array of dimension NYD containing the y coordinates
*          of the input-grid data points (must be in a
*          monotonic increasing order),
*   ZD   = two-dimensional array of dimension NXD*NYD
*          containing the z(x,y) values at the input-grid data
*          points,
*   PDD  = three-dimensional array of dimension 3*NXD*NYD
*          containing the estimated zx, zy, and zxy values
*          at the input-grid data points,
*   NIP  = number of the output points at which interpolation
*          is to be performed,
*   XI   = array of dimension NIP containing the x coordinates
*          of the output points,
*   YI   = array of dimension NIP containing the y coordinates
*          of the output points,
*   INXI = integer array of dimension NIP containing the
*          interval numbers of the input grid intervals in the
*          x direction where the x coordinates of the output
*          points lie,
*   INYI = integer array of dimension NIP containing the
*          interval numbers of the input grid intervals in the
*          y direction where the y coordinates of the output
*          points lie.
*
* The output argument is
*   ZI   = array of dimension NIP, where the interpolated z
*          values at the output points are to be stored.
*  PDZI  = array of dimension (2,NIP) where the interpolated 
*          partial derivatives D_(x,y) at XI,YI are to be stored.
*
* Specification statements
*     .. Scalar Arguments ..
      INTEGER          NIP,NXD,NYD
*     ..
*     .. Array Arguments ..
      REAL*8           PDD(3,NXD,NYD),XD(NXD),XI(NIP),YD(NYD),YI(NIP),
     +                 ZD(NXD,NYD),ZI(NIP),PDZI(2,NIP)
      INTEGER          INXI(NIP),INYI(NIP)
*     ..
*     .. Local Scalars ..
      REAL*8           A,B,C,D,DX,DXSQ,DY,DYSQ,P00,P01,P02,P03,P10,P11,
     +                 P12,P13,P20,P21,P22,P23,P30,P31,P32,P33,
     +                 Q0,Q1,Q2,Q3,Q0Y,Q1Y,Q2Y,Q3Y,
     +                 U,V,X0,XII,Y0,YII,Z00,Z01,Z0DX,Z0DY,Z10,Z11,
     +                 Z1DX,Z1DY,ZDXDY,ZII,ZIIX,ZIIY,
     +                 ZX00,ZX01,ZX0DY,ZX10,ZX11,
     +                 ZX1DY,ZXY00,ZXY01,ZXY10,ZXY11,ZY00,ZY01,ZY0DX,
     +                 ZY10,ZY11,ZY1DX
      INTEGER          IIP,IXD0,IXD1,IXDI,IXDIPV,IYD0,IYD1,IYDI,IYDIPV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MAX
*     ..
* Calculation
* Outermost DO-loop with respect to the output point
      DO 10 IIP = 1,NIP
          XII = XI(IIP)
          YII = YI(IIP)
          IF (IIP.EQ.1) THEN
              IXDIPV = -1
              IYDIPV = -1
          ELSE
              IXDIPV = IXDI
              IYDIPV = IYDI
          END IF
          IXDI = INXI(IIP)
          IYDI = INYI(IIP)
* Retrieves the z and partial derivative values at the origin of
* the coordinate for the rectangle.
          IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
              IXD0 = MAX(1,IXDI)
              IYD0 = MAX(1,IYDI)
              X0 = XD(IXD0)
              Y0 = YD(IYD0)
              Z00 = ZD(IXD0,IYD0)
              ZX00 = PDD(1,IXD0,IYD0)
              ZY00 = PDD(2,IXD0,IYD0)
              ZXY00 = PDD(3,IXD0,IYD0)
          END IF
* Case 1.  When the rectangle is inside the data area in both the
* x and y directions.
          IF ((IXDI.GT.0.AND.IXDI.LT.NXD) .AND.
     +        (IYDI.GT.0.AND.IYDI.LT.NYD)) THEN
* Retrieves the z and partial derivative values at the other three
* vertexes of the rectangle.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  IXD1 = IXD0 + 1
                  DX = XD(IXD1) - X0
                  DXSQ = DX*DX
                  IYD1 = IYD0 + 1
                  DY = YD(IYD1) - Y0
                  DYSQ = DY*DY
                  Z10 = ZD(IXD1,IYD0)
                  Z01 = ZD(IXD0,IYD1)
                  Z11 = ZD(IXD1,IYD1)
                  ZX10 = PDD(1,IXD1,IYD0)
                  ZX01 = PDD(1,IXD0,IYD1)
                  ZX11 = PDD(1,IXD1,IYD1)
                  ZY10 = PDD(2,IXD1,IYD0)
                  ZY01 = PDD(2,IXD0,IYD1)
                  ZY11 = PDD(2,IXD1,IYD1)
                  ZXY10 = PDD(3,IXD1,IYD0)
                  ZXY01 = PDD(3,IXD0,IYD1)
                  ZXY11 = PDD(3,IXD1,IYD1)
* Calculates the polynomial coefficients.
                  Z0DX = (Z10-Z00)/DX
                  Z1DX = (Z11-Z01)/DX
                  Z0DY = (Z01-Z00)/DY
                  Z1DY = (Z11-Z10)/DY
                  ZX0DY = (ZX01-ZX00)/DY
                  ZX1DY = (ZX11-ZX10)/DY
                  ZY0DX = (ZY10-ZY00)/DX
                  ZY1DX = (ZY11-ZY01)/DX
                  ZDXDY = (Z1DY-Z0DY)/DX
                  A = ZDXDY - ZX0DY - ZY0DX + ZXY00
                  B = ZX1DY - ZX0DY - ZXY10 + ZXY00
                  C = ZY1DX - ZY0DX - ZXY01 + ZXY00
                  D = ZXY11 - ZXY10 - ZXY01 + ZXY00
                  P00 = Z00
                  P01 = ZY00
                  P02 = (2.0* (Z0DY-ZY00)+Z0DY-ZY01)/DY
                  P03 = (-2.0*Z0DY+ZY01+ZY00)/DYSQ
                  P10 = ZX00
                  P11 = ZXY00
                  P12 = (2.0* (ZX0DY-ZXY00)+ZX0DY-ZXY01)/DY
                  P13 = (-2.0*ZX0DY+ZXY01+ZXY00)/DYSQ
                  P20 = (2.0* (Z0DX-ZX00)+Z0DX-ZX10)/DX
                  P21 = (2.0* (ZY0DX-ZXY00)+ZY0DX-ZXY10)/DX
                  P22 = (3.0* (3.0*A-B-C)+D)/ (DX*DY)
                  P23 = (-6.0*A+2.0*B+3.0*C-D)/ (DX*DYSQ)
                  P30 = (-2.0*Z0DX+ZX10+ZX00)/DXSQ
                  P31 = (-2.0*ZY0DX+ZXY10+ZXY00)/DXSQ
                  P32 = (-6.0*A+3.0*B+2.0*C-D)/ (DXSQ*DY)
                  P33 = (2.0* (2.0*A-B-C)+D)/ (DXSQ*DYSQ)
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V* (P01+V* (P02+V*P03))
              Q1 = P10 + V* (P11+V* (P12+V*P13))
              Q2 = P20 + V* (P21+V* (P22+V*P23))
              Q3 = P30 + V* (P31+V* (P32+V*P33))
              ZII = Q0 + U* (Q1+U* (Q2+U*Q3))
C Evaluates partial derivatives at ZI of the polynomial.
              Q0Y = P01 + V*(2.0*P02+3.0*P03*V)
              Q1Y = P11 + V*(2.0*P12+3.0*P13*V)
              Q2Y = P21 + V*(2.0*P22+3.0*P23*V)
              Q3Y = P31 + V*(2.0*P32+3.0*P33*V)
              ZIIX = Q1 + U*(2.0*Q2+3.0*Q3*U)
              ZIIY = Q0Y + U* (Q1Y+U*(Q2Y+U*Q3Y))
* End of Case 1
* Case 2.  When the rectangle is inside the data area in the x
* direction but outside in the y direction.
          ELSE IF ((IXDI.GT.0.AND.IXDI.LT.NXD) .AND.
     +             (IYDI.LE.0.OR.IYDI.GE.NYD)) THEN
* Retrieves the z and partial derivative values at the other
* vertex of the semi-infinite rectangle.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  IXD1 = IXD0 + 1
                  DX = XD(IXD1) - X0
                  DXSQ = DX*DX
                  Z10 = ZD(IXD1,IYD0)
                  ZX10 = PDD(1,IXD1,IYD0)
                  ZY10 = PDD(2,IXD1,IYD0)
                  ZXY10 = PDD(3,IXD1,IYD0)
* Calculates the polynomial coefficients.
                  Z0DX = (Z10-Z00)/DX
                  ZY0DX = (ZY10-ZY00)/DX
                  P00 = Z00
                  P01 = ZY00
                  P10 = ZX00
                  P11 = ZXY00
                  P20 = (2.0* (Z0DX-ZX00)+Z0DX-ZX10)/DX
                  P21 = (2.0* (ZY0DX-ZXY00)+ZY0DX-ZXY10)/DX
                  P30 = (-2.0*Z0DX+ZX10+ZX00)/DXSQ
                  P31 = (-2.0*ZY0DX+ZXY10+ZXY00)/DXSQ
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V*P01
              Q1 = P10 + V*P11
              Q2 = P20 + V*P21
              Q3 = P30 + V*P31
              ZII = Q0 + U* (Q1+U* (Q2+U*Q3))
C Evaluates partial derivatives at ZI of the polynomial.
              Q0Y = P01
              Q1Y = P11
              Q2Y = P21
              Q3Y = P31
              ZIIX = Q1 + U*(2.0*Q2+3.0*Q3*U)
              ZIIY = Q0Y + U* (Q1Y+U*(Q2Y+U*Q3Y))
* End of Case 2
* Case 3.  When the rectangle is outside the data area in the x
* direction but inside in the y direction.
          ELSE IF ((IXDI.LE.0.OR.IXDI.GE.NXD) .AND.
     +             (IYDI.GT.0.AND.IYDI.LT.NYD)) THEN
* Retrieves the z and partial derivative values at the other
* vertex of the semi-infinite rectangle.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  IYD1 = IYD0 + 1
                  DY = YD(IYD1) - Y0
                  DYSQ = DY*DY
                  Z01 = ZD(IXD0,IYD1)
                  ZX01 = PDD(1,IXD0,IYD1)
                  ZY01 = PDD(2,IXD0,IYD1)
                  ZXY01 = PDD(3,IXD0,IYD1)
* Calculates the polynomial coefficients.
                  Z0DY = (Z01-Z00)/DY
                  ZX0DY = (ZX01-ZX00)/DY
                  P00 = Z00
                  P01 = ZY00
                  P02 = (2.0* (Z0DY-ZY00)+Z0DY-ZY01)/DY
                  P03 = (-2.0*Z0DY+ZY01+ZY00)/DYSQ
                  P10 = ZX00
                  P11 = ZXY00
                  P12 = (2.0* (ZX0DY-ZXY00)+ZX0DY-ZXY01)/DY
                  P13 = (-2.0*ZX0DY+ZXY01+ZXY00)/DYSQ
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V* (P01+V* (P02+V*P03))
              Q1 = P10 + V* (P11+V* (P12+V*P13))
              ZII = Q0 + U*Q1
C Evaluates partial derivatives at ZI of the polynomial.
              Q0Y = P01 + V*(2.0*P02+3.0*P03*V)
              Q1Y = P11 + V*(2.0*P12+3.0*P13*V)
              ZIIX = Q1
              ZIIY = Q0Y + U*Q1Y
* End of Case 3
* Case 4.  When the rectangle is outside the data area in both the
* x and y direction.
          ELSE IF ((IXDI.LE.0.OR.IXDI.GE.NXD) .AND.
     +             (IYDI.LE.0.OR.IYDI.GE.NYD)) THEN
* Calculates the polynomial coefficients.
              IF (IXDI.NE.IXDIPV .OR. IYDI.NE.IYDIPV) THEN
                  P00 = Z00
                  P01 = ZY00
                  P10 = ZX00
                  P11 = ZXY00
              END IF
* Evaluates the polynomial.
              U = XII - X0
              V = YII - Y0
              Q0 = P00 + V*P01
              Q1 = P10 + V*P11
              ZII = Q0 + U*Q1
C Evaluates partial derivatives at ZI of the polynomial.
              Q0Y = P01
              Q1Y = P11
              ZIIX = Q1
              ZIIY = Q0Y + U*Q1Y
          END IF
* End of Case 4
          ZI(IIP)     = ZII
          PDZI(1,IIP) = ZIIX
          PDZI(2,IIP) = ZIIY
   10 CONTINUE
      RETURN
      END

C *****************************************************************************

C mi,mr GRID SET UP AND INTERVAL SEARCH ROUTINE required by DRV_SPECOPT_LBFGSB
C This routine reduces the large 6*2 mi,mr grid square to a smaller rectangle whose
C vertices are the nearest mi,mr table values that bound the 'candidates' mi,mr from the
C optimization iteration.  All in 2D space of course (x-y)

      SUBROUTINE DRGLCTN(NXD,NYD,XD,YD,NIP,XI,YI, INXI,INYI)  !'desired rectangular grid location'
                                                              !or 'drive grid location'

* Location of the desired points in a rectangular grid
* (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08
*
* The grid lines can be unevenly spaced.
*
* The input arguments are
*   NXD  = number of the input-grid data points in the x
*          coordinate (must be 2 or greater),
*   NYD  = number of the input-grid data points in the y
*          coordinate (must be 2 or greater),
*   XD   = array of dimension NXD containing the x coordinates
*          of the input-grid data points (must be in a
*          monotonic increasing order),
*   YD   = array of dimension NYD containing the y coordinates
*          of the input-grid data points (must be in a
*          monotonic increasing order),
*   NIP  = number of the output points to be located (must be
*          1 or greater),
*   XI   = array of dimension NIP containing the x coordinates
*          of the output points to be located,
*   YI   = array of dimension NIP containing the y coordinates
*          of the output points to be located.
*
* The output arguments are
*   INXI = integer array of dimension NIP where the interval
*          numbers of the XI array elements are to be stored,
*   INYI = integer array of dimension NIP where the interval
*          numbers of the YI array elements are to be stored.
* The interval numbers are between 0 and NXD and between 0 and NYD,
* respectively.
*
*
* Specification statements
*     .. Scalar Arguments ..
      INTEGER          NIP,NXD,NYD
*     ..
*     .. Array Arguments ..
      REAL*8           XD(NXD),XI(NIP),YD(NYD),YI(NIP)
      INTEGER          INXI(NIP),INYI(NIP)
*     ..
*     .. Local Scalars ..
      REAL*8           XII,YII
      INTEGER          IIP,IMD,IMN,IMX,IXD,IYD,NINTX,NINTY
*     ..
* DO-loop with respect to IIP, which is the point number of the
* output point
      DO 30 IIP = 1,NIP
          XII = XI(IIP)
          YII = YI(IIP)
* Checks if the x coordinate of the IIPth output point, XII, is
* in a new interval.  (NINTX is the new-interval flag.)
          IF (IIP.EQ.1) THEN
              NINTX = 1
          ELSE
              NINTX = 0
              IF (IXD.EQ.0) THEN
                  IF (XII.GT.XD(1)) NINTX = 1
              ELSE IF (IXD.LT.NXD) THEN
                  IF ((XII.LT.XD(IXD)) .OR.
     +                (XII.GT.XD(IXD+1))) NINTX = 1
              ELSE
                  IF (XII.LT.XD(NXD)) NINTX = 1
              END IF
          END IF
* Locates the output point by binary search if XII is in a new
* interval.  Determines IXD for which XII lies between XD(IXD)
* and XD(IXD+1).
          IF (NINTX.EQ.1) THEN
              IF (XII.LE.XD(1)) THEN
                  IXD = 0
              ELSE IF (XII.LT.XD(NXD)) THEN
                  IMN = 1
                  IMX = NXD
                  IMD = (IMN+IMX)/2
   10             IF (XII.GE.XD(IMD)) THEN
                      IMN = IMD
                  ELSE
                      IMX = IMD
                  END IF
                  IMD = (IMN+IMX)/2
                  IF (IMD.GT.IMN) GO TO 10
                  IXD = IMD
              ELSE
                  IXD = NXD
              END IF
          END IF
          INXI(IIP) = IXD
* Checks if the y coordinate of the IIPth output point, YII, is
* in a new interval.  (NINTY is the new-interval flag.)
          IF (IIP.EQ.1) THEN
              NINTY = 1
          ELSE
              NINTY = 0
              IF (IYD.EQ.0) THEN
                  IF (YII.GT.YD(1)) NINTY = 1
              ELSE IF (IYD.LT.NYD) THEN
                  IF ((YII.LT.YD(IYD)) .OR.
     +                (YII.GT.YD(IYD+1))) NINTY = 1
              ELSE
                  IF (YII.LT.YD(NYD)) NINTY = 1
              END IF
          END IF
* Locates the output point by binary search if YII is in a new
* interval.  Determines IYD for which YII lies between YD(IYD)
* and YD(IYD+1).
          IF (NINTY.EQ.1) THEN
              IF (YII.LE.YD(1)) THEN
                  IYD = 0
              ELSE IF (YII.LT.YD(NYD)) THEN
                  IMN = 1
                  IMX = NYD
                  IMD = (IMN+IMX)/2
   20             IF (YII.GE.YD(IMD)) THEN
                      IMN = IMD
                  ELSE
                      IMX = IMD
                  END IF
                  IMD = (IMN+IMX)/2
                  IF (IMD.GT.IMN) GO TO 20
                  IYD = IMD
              ELSE
                  IYD = NYD
              END IF
          END IF
          INYI(IIP) = IYD
   30 CONTINUE

      RETURN
      END


C ############################################################################
C ############################################################################
C ############################################################################
C ############################################################################

C SECTION 5:

C Partial derivative routines.
C Note that these routines are the EXPLICIT formulations of the partial 
C derivatives as indicated by the modelling approach used (eg., GSM01
C for the water model).  If we change a particular model then we must
C recode the derivatives in this section.

C **********************************************************************

C                   ***** Master derivative routine *****
C Any change to DFCNLSQ_GRAD Master Function will require modifications to exact 
C derivatives in respective called routines here.
        
C Computes partial derivatives D(RhoA+rhow)/D(mi) & D(RhoA+rhow)/D(mr),
C given computed partial derivatives in auxillary functions.

      SUBROUTINE DPDR_MIMR ( NBand, zwRho, pdZrho,
     &                       zTau865,pdZtau865,zTau,
     &                       zC865,pdZc865,
     &                       t0,zA0,zB0,pdZa0,pdZb0,
     &                       t ,zA ,zB ,pdZa ,pdZb ,
     &                       pdRmi,pdRmr )

c input: 
c        NBand - wavelength index [1::6]
c        zwRho - candidate rhow, as a function of candidates C,acdm,bbp (updates)
c       pdZrho - interpolated partial derivative rhoA   (as a function of candidate (optimal) mi,mr)
c      zTau865 - interpolated tau865                    (                   "                      )
c    pdZtau865 - interpolated partial derivative tau865 (                   "                      )
c         zTau - interpolated tau(lambda)               (                   "                      )
c        zC865 - interpolated tau(lambda)/tau(865)      (                   "                      )
c      pdZc865 - interpolated partial derivative tau(lambda)/tau865 (       "                      )
c           t0 - interpolated t(inbound)                (                   "                      )
c          zA0 - interpolated diff trans A coeff at theta0 (                "                      )               
c          zB0 - interpolated diff trans B coeff at theta0 (                "                      )           
c        pdZa0 - interpolated partial derivative diff trans A coeff at theta0 (          "         )     
c        pdZb0 - interpolated partial derivative diff trans B coeff at theta0 (          "         )     
c            t - t(outbound), interpolated for candidate mi,mr
c           zA - interpolated diff trans A coeff at theta  (                "                      )           
c           zB - interpolated diff trans B coeff at theta  (                "                      )           
c         pdZa - interpolated partial derivative diff trans A coeff at theta  (          "         )     
c         pdZb - interpolated partial derivative diff trans B coeff at theta  (          "         )     

c output:  
c        pdRmi - interpolated partial derivative rhoA+rhow (as a function of candidate (optimal) mi)
c        pdRmr - interpolated partial derivative rhoA+rhow (as a function of candidate (optimal) mr)


C ***** Master derivative routine *****
C Any change to Master Function will require modifications to exact derivatives
C in respective routines here        



C Calculations:
C ------------          NBand
C                       ____
C              D Rho    \\    D Slsq  |                  D Rho_j
C              ------ = /     ------  |                * -------
C              D m      ----  D Rho_j | Rho_j=Rho_j(m)   D m
C                       j=1
C  where,
C              Rho_j = RhoCand_j + trans0_j*RhoWat_j*trans_j (j=Lambda) 
C
C
C  For each j=Lambda and m = (mi or mr):
C
C   D Rho    D RhoCand            / D trans0                    D trans \
C   ------ = --------- + RhoWat * | -------- * trans + trans0 * ------- |
C   D m      D m                  \ D m                         D m     /
C
C  where,
C          D trans0   D a0     b0*tau            / D b0              D tau  \
C          -------- = ----- * e       + trans0 * | ---- * tau + b0 * ------ |
C          D m        D m                        \ D m               D m    /
C
C  and same for D_trans/D_m.
C
C  Assuming Z(x,y) is a bivariate interpolating function, such that
C
C  D RhoCand   D Z(RhoCand)
C  --------- = ------------
C  D m         D m
C
C  D tau       D Z(tau)                    D Z(cc865)
C  ------    = -------- * cc865 + tau865 * ----------
C  D m         D m                         D m
C
C  D a0        D Z(a0)
C  -----     = ------- , and the same for b0,a,b - diffuse trans coeffs.
C  D m         D m
C
C
      IMPLICIT NONE
      INTEGER    NLAMBDA
      PARAMETER (NLAMBDA=8)

      INTEGER  NBand
      REAL*8   pdZrho(2,NLAMBDA),zwRho(NLAMBDA),       ! [1,2]=[mi,mr]=x-y reference of mi,mr interpolation
     &         zTau865,pdZtau865(2),zTau(NLAMBDA),     ! [1,2]=[mi,mr]
     &         zC865(NLAMBDA),pdZc865(2,NLAMBDA),      ! [1,2]=[mi,mr]
     &         t0(NLAMBDA),
     &         zA0(NLAMBDA),zB0(NLAMBDA),              ! [1,2]=[mr,mi]
     &         pdZa0(2,NLAMBDA),pdZb0(2,NLAMBDA),      ! [1,2]=[mr,mi]
     &         t(NLAMBDA), 
     &         zA(NLAMBDA), zB(NLAMBDA),               ! [1,2]=[mr,mi]
     &         pdZa(2,NLAMBDA), pdZb(2,NLAMBDA)        ! [1,2]=[mr,mi]

      REAL*8   pdRmi(NLAMBDA),pdRmr(NLAMBDA)

      INTEGER  jl
      REAL*8   dtau_dm,
     &         da_dm,db_dm,da0_dm,db0_dm,dt0_dm,dt_dm,
     &         btau,b0tau

c code below duplicates math given above (for candidates mi,mr)

C --------------------------------------------- PARTIAL DERIVATIVE IN mi
      DO jl=1,NBand

        dtau_dm = pdZtau865(1)*zC865(jl) + zTau865*pdZc865(1,jl)
        da0_dm  = pdZa0(2,jl)
        db0_dm  = pdZb0(2,jl)
        da_dm   = pdZa(2,jl)
        db_dm   = pdZb(2,jl)
        b0tau   = zB0(jl)*zTau(jl)
        btau    = zB(jl)*zTau(jl)

        dt0_dm  = da0_dm * DEXP(b0tau)  +
     &            t0(jl) * (db0_dm*zTau(jl) + zB0(jl)*dtau_dm)
        dt_dm   = da_dm  * DEXP(btau)   +
     &            t(jl)  * (db_dm*zTau(jl)  + zB(jl)*dtau_dm)

        pdRmi(jl) = pdZrho(1,jl) + zwRho(jl)*(dt0_dm*t(jl)+t0(jl)*dt_dm)

      ENDDO  ! jl
C --------------------------------------------- PARTIAL DERIVATIVE IN mr
      DO jl=1,NBand

        dtau_dm = pdZtau865(2)*zC865(jl) + zTau865*pdZc865(2,jl)
        da0_dm  = pdZa0(1,jl)
        db0_dm  = pdZb0(1,jl)
        da_dm   = pdZa(1,jl)
        db_dm   = pdZb(1,jl)
        b0tau   = zB0(jl)*zTau(jl)
        btau    = zB(jl)*zTau(jl)

        dt0_dm  = da0_dm * DEXP(b0tau)  +
     &            t0(jl) * (db0_dm*zTau(jl) + zB0(jl)*dtau_dm)
        dt_dm   = da_dm  * DEXP(btau)   +
     &            t(jl)  * (db_dm*zTau(jl)  + zB(jl)*dtau_dm)

        pdRmr(jl) = pdZrho(2,jl) + zwRho(jl)*(dt0_dm*t(jl)+t0(jl)*dt_dm)

      ENDDO  ! jl

      RETURN
      END

C **********************************************************************

C                     ***** Master derivative routine *****
C Any change to DFCNLSQ_GRAD Master Function will require modifications to exact 
C derivatives in respective called routines here.   
            
C Purpose:  Computes partial derivatives D(LwN)/D(C,acdm,bbp) based
C           on GS97 water reflectance model.

      SUBROUTINE DPDLWN_GS97 ( NBand,
     &                         C,acdm,bbp,
     &                         pdLc,pdLacdm,pdLbbp )

c: input: 
c          NBand - wavelength index [1::6]
c     C,acdm,bbp - candidates from Lw

c output:
c           pdLc - partial derivative LwN as a function of C
c        pdLacdm -                   "                     acdm
c         pdLbbp -                   "                     bbp


C Calculations:
C ------------
C
C GS97 LwN Model (j=Lambda):
C
C         (j)         (j)                        2
C        L    = g0 * F    * [ g1 * x   +  g2 * x    ]
C         wN          0             (j)         (j)
C where,
C
C                         (j)          bbp_s
C                        b    + (443/j)     * BBP
C                         bw
C  x   = -----------------------------------------------------------------------
C   (j)    (j)         (j)       -adm_s (j-443)          (j)          bbp_s
C        aw   + aph_lin   * C + e              * Acdm + b    + (443/j)     * BBP
C                                                        bw
C
C Partial Derivatives:
C -------------------
C                       (j)    (j)         (j)
C      assuming:      CC   = aw   + aph_lin   * C
C
C                       (j)     (j)         bbp_s
C                     BB   = bbw   + (443/j)     * BBP
C
C                       (j)   -adm_s (j-443)
C                     AA   = e              * Acdm
C      obtaining:
C
C            (j)                                        (j)
C       D LwN              (j)                       D x
C       --------- = g0 * F0   * [ g1 + g2 * 2x   ] * ---------
C       D (C,a,b)                             (j)    D (C,a,b)
C
C      where,
C                (j)                (j)    (j)
C             D x            aph_lin   * BB
C             ------ = - ------------------------
C             D C            (j)    (j)    (j)  2
C                        ( CC   + AA   + BB    )
C
C
C                (j)       -adm_s*(j-443)    (j)
C             D x         e              * BB
C             ------ = - ------------------------
C             D a            (j)    (j)    (j)  2
C                        ( CC   + AA   + BB    )
C
C
C                (j)          bbp_s  /           (j)        \
C             D x      / 443 \       |      1 - x           |
C             ------ = | --- |     * | -------------------- |
C             D b      \  j  /       |   (j)    (j)    (j)  |
C                                    \ CC   + AA   + BB     /
C

      IMPLICIT NONE
      INTEGER    NLAMBDA
      PARAMETER (NLAMBDA=8)

      INTEGER  NBand
      REAL*8   C,acdm,bbp

      REAL*8   pdLc(NLAMBDA),pdLacdm(NLAMBDA),pdLbbp(NLAMBDA)

      INTEGER  jl
      REAL*8   lambda(NLAMBDA),Fobar(NLAMBDA)

      REAL*8   admstar,bbpstar,
     &         g0,g1,g2,x,
     &         AA,BB,CC,DD,GG,
     &         dxda,dxdb,dxdc

      DATA g0/0.5238D0/, g1/0.0949D0/, g2/0.0794D0/

      DATA lambda  /412.D0, 443.D0, 490.D0, 510.D0,
     &              555.D0, 670.D0, 765.D0, 865.D0 /

      REAL*8 d_adm_s,d_bbp_s,
     &       d_aw(NLAMBDA),d_bbw(NLAMBDA),d_aph_lin(NLAMBDA)
      COMMON /D_WATERPAR/ d_adm_s,d_bbp_s,d_aw,d_bbw,d_aph_lin

c      DATA Fobar   /170.79D0, 189.45D0, 193.66D0, 188.35D0,
c     &              185.33D0, 153.D0  , 0.D0,     0.D0     /

      DATA Fobar   /173.004D0, 190.154D0, 196.473D0, 188.158D0,
     &              183.01D0,  151.143D0  , 0.D0,     0.D0     /

      

      DO jl=1,NBand

        admstar = DEXP(-d_adm_s*(lambda(jl) - 443.D0))
        bbpstar = (443.D0/lambda(jl))**d_bbp_s

        AA = admstar*acdm
        BB = d_bbw(jl) + bbpstar*bbp 
        CC = d_aw(jl) + d_aph_lin(jl)*C
        DD = AA+BB+CC
        x  = BB/DD
        GG = g0*Fobar(jl)*(g1+2.D0*g2*x)

        dxdc = -d_aph_lin(jl)*BB/DD**2
        dxda = -admstar*BB/DD**2
        dxdb = bbpstar*(1.D0-x)/DD

        pdLc(jl)    = GG*dxdc
        pdLacdm(jl) = GG*dxda
        pdLbbp(jl)  = GG*dxdb

      ENDDO

      RETURN
      END
C **********************************************************************

C THIS SUBROUTINE DOES NOT NEED TO BE ALTERED - UNIVERSAL!!!
C Estimates three partial derivatives, zx, zy, and zxy, of a bivariate 
C function, z(x,y), on a rectangular grid in
C the x-y plane.  It is based on the revised Akima method that has
C the accuracy of a bicubic polynomial. 
C Called by main optimisation SUBROUTINE DRV_SPECOPT_LBFGSB to calculate
C partial derivatives of rhoA(lambda),tau(lambda/tau(865),diff_trans_coeffs(lambda),
C tau865, ssalbedo and nu, all as a function of 'candidates' mi,mr.

      SUBROUTINE DRGPD3P(NXD,NYD,XD,YD,ZD, PDD)
*
* Partial derivatives of a bivariate function on a rectangular grid
* (a supporting subroutine of the RGBI3P/RGSF3P subroutine package)
*
* Hiroshi Akima
* U.S. Department of Commerce, NTIA/ITS
* Version of 1995/08

* The input arguments are
*   NXD = number of the input-grid data points in the x
*         coordinate (must be 2 or greater),
*   NYD = number of the input-grid data points in the y
*         coordinate (must be 2 or greater),
*   XD  = array of dimension NXD containing the x coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   YD  = array of dimension NYD containing the y coordinates
*         of the input-grid data points (must be in a
*         monotonic increasing order),
*   ZD  = two-dimensional array of dimension NXD*NYD
*         containing the z(x,y) values at the input-grid data
*         points.
*
* The output argument is
*   PDD = three-dimensional array of dimension 3*NXD*NYD,
*         where the estimated zx, zy, and zxy values at the
*         input-grid data points are to be stored.
*
*
* Specification statements
*     .. Scalar Arguments ..
      INTEGER          NXD,NYD
*     ..
*     .. Array Arguments ..
      REAL*8           PDD(3,NXD,NYD),XD(NXD),YD(NYD),ZD(NXD,NYD)
*     ..
*     .. Local Scalars ..
      REAL*8           B00,B00X,B00Y,B01,B10,B11,CX1,CX2,CX3,CY1,CY2,
     +                 CY3,DISF,DNM,DZ00,DZ01,DZ02,DZ03,DZ10,DZ11,DZ12,
     +                 DZ13,DZ20,DZ21,DZ22,DZ23,DZ30,DZ31,DZ32,DZ33,
     +                 DZX10,DZX20,DZX30,DZXY11,DZXY12,DZXY13,DZXY21,
     +                 DZXY22,DZXY23,DZXY31,DZXY32,DZXY33,DZY01,DZY02,
     +                 DZY03,EPSLN,PEZX,PEZXY,PEZY,SMPEF,SMPEI,SMWTF,
     +                 SMWTI,SX,SXX,SXXY,SXXYY,SXY,SXYY,SXYZ,SXZ,SY,SYY,
     +                 SYZ,SZ,VOLF,WT,X0,X1,X2,X3,XX1,XX2,XX3,Y0,Y1,Y2,
     +                 Y3,Z00,Z01,Z02,Z03,Z10,Z11,Z12,Z13,Z20,Z21,Z22,
     +                 Z23,Z30,Z31,Z32,Z33,ZXDI,ZXYDI,ZYDI,ZZ0,ZZ1,ZZ2
      INTEGER          IPEX,IPEY,IX0,IX1,IX2,IX3,IY0,IY1,IY2,IY3,JPEXY,
     +                 JXY,NX0,NY0
*     ..
*     .. Local Arrays ..
      REAL*8           B00XA(4),B00YA(4),B01A(4),B10A(4),CXA(3,4),
     +                 CYA(3,4),SXA(4),SXXA(4),SYA(4),SYYA(4),XA(3,4),
     +                 YA(3,4),Z0IA(3,4),ZI0A(3,4)
      INTEGER          IDLT(3,4)
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        MAX
*     ..
*     .. Statement Functions ..
      REAL*8           Z2F,Z3F
*     ..
* Data statements 
      DATA             ((IDLT(JXY,JPEXY),JPEXY=1,4),JXY=1,3)/-3,-2,-1,1,
     +                 -2,-1,1,2,-1,1,2,3/
*     ..
* Statement Function definitions 
      Z2F(XX1,XX2,ZZ0,ZZ1) = (ZZ1-ZZ0)*XX2/XX1 + ZZ0
      Z3F(XX1,XX2,XX3,ZZ0,ZZ1,ZZ2) = ((ZZ2-ZZ0)* (XX3-XX1)/XX2-
     +                               (ZZ1-ZZ0)* (XX3-XX2)/XX1)*
     +                               (XX3/ (XX2-XX1)) + ZZ0
*     ..
* Calculation
* Initial setting of some local variables
      NX0 = MAX(4,NXD)
      NY0 = MAX(4,NYD)
* Double DO-loop with respect to the input grid points
      DO 60 IY0 = 1,NYD
          DO 50 IX0 = 1,NXD
              X0 = XD(IX0)
              Y0 = YD(IY0)
              Z00 = ZD(IX0,IY0)
* Part 1.  Estimation of ZXDI
* Initial setting
              SMPEF = 0.0
              SMWTF = 0.0
              SMPEI = 0.0
              SMWTI = 0.0
* DO-loop with respect to the primary estimate
              DO 10 IPEX = 1,4
* Selects necessary grid points in the x direction.
                  IX1 = IX0 + IDLT(1,IPEX)
                  IX2 = IX0 + IDLT(2,IPEX)
                  IX3 = IX0 + IDLT(3,IPEX)
                  IF ((IX1.LT.1) .OR. (IX2.LT.1) .OR. (IX3.LT.1) .OR.
     +                (IX1.GT.NX0) .OR. (IX2.GT.NX0) .OR.
     +                (IX3.GT.NX0)) GO TO 10
* Selects and/or supplements the x and z values.
                  X1 = XD(IX1) - X0
                  Z10 = ZD(IX1,IY0)
                  IF (NXD.GE.4) THEN
                      X2 = XD(IX2) - X0
                      X3 = XD(IX3) - X0
                      Z20 = ZD(IX2,IY0)
                      Z30 = ZD(IX3,IY0)
                  ELSE IF (NXD.EQ.3) THEN
                      X2 = XD(IX2) - X0
                      Z20 = ZD(IX2,IY0)
                      X3 = 2*XD(3) - XD(2) - X0
                      Z30 = Z3F(X1,X2,X3,Z00,Z10,Z20)
                  ELSE IF (NXD.EQ.2) THEN
                      X2 = 2*XD(2) - XD(1) - X0
                      Z20 = Z2F(X1,X2,Z00,Z10)
                      X3 = 2*XD(1) - XD(2) - X0
                      Z30 = Z2F(X1,X3,Z00,Z10)
                  END IF
                  DZX10 = (Z10-Z00)/X1
                  DZX20 = (Z20-Z00)/X2
                  DZX30 = (Z30-Z00)/X3
* Calculates the primary estimate of partial derivative zx as
* the coefficient of the bicubic polynomial.
                  CX1 = X2*X3/ ((X1-X2)* (X1-X3))
                  CX2 = X3*X1/ ((X2-X3)* (X2-X1))
                  CX3 = X1*X2/ ((X3-X1)* (X3-X2))
                  PEZX = CX1*DZX10 + CX2*DZX20 + CX3*DZX30
* Calculates the volatility factor and distance factor in the x
* direction for the primary estimate of zx.
                  SX = X1 + X2 + X3
                  SZ = Z00 + Z10 + Z20 + Z30
                  SXX = X1*X1 + X2*X2 + X3*X3
                  SXZ = X1*Z10 + X2*Z20 + X3*Z30
                  DNM = 4.0*SXX - SX*SX
                  B00 = (SXX*SZ-SX*SXZ)/DNM
                  B10 = (4.0*SXZ-SX*SZ)/DNM
                  DZ00 = Z00 - B00
                  DZ10 = Z10 - (B00+B10*X1)
                  DZ20 = Z20 - (B00+B10*X2)
                  DZ30 = Z30 - (B00+B10*X3)
                  VOLF = DZ00**2 + DZ10**2 + DZ20**2 + DZ30**2
                  DISF = SXX
* Calculates the EPSLN value, which is used to decide whether or
* not the volatility factor is essentially zero.
                  EPSLN = (Z00**2+Z10**2+Z20**2+Z30**2)*1.0E-12
* Accumulates the weighted primary estimates of zx and their
* weights.
                  IF (VOLF.GT.EPSLN) THEN
* - For a finite weight.
                      WT = 1.0/ (VOLF*DISF)
                      SMPEF = SMPEF + WT*PEZX
                      SMWTF = SMWTF + WT
                  ELSE
* - For an infinite weight.
                      SMPEI = SMPEI + PEZX
                      SMWTI = SMWTI + 1.0
                  END IF
* Saves the necessary values for estimating zxy
                  XA(1,IPEX) = X1
                  XA(2,IPEX) = X2
                  XA(3,IPEX) = X3
                  ZI0A(1,IPEX) = Z10
                  ZI0A(2,IPEX) = Z20
                  ZI0A(3,IPEX) = Z30
                  CXA(1,IPEX) = CX1
                  CXA(2,IPEX) = CX2
                  CXA(3,IPEX) = CX3
                  SXA(IPEX) = SX
                  SXXA(IPEX) = SXX
                  B00XA(IPEX) = B00
                  B10A(IPEX) = B10
   10         CONTINUE
* Calculates the final estimate of zx.
              IF (SMWTI.LT.0.5) THEN
* - When no infinite weights exist.
                  ZXDI = SMPEF/SMWTF
              ELSE
* - When infinite weights exist.
                  ZXDI = SMPEI/SMWTI
              END IF
* End of Part 1.
* Part 2.  Estimation of ZYDI
* Initial setting
              SMPEF = 0.0
              SMWTF = 0.0
              SMPEI = 0.0
              SMWTI = 0.0
* DO-loop with respect to the primary estimate
              DO 20 IPEY = 1,4
* Selects necessary grid points in the y direction.
                  IY1 = IY0 + IDLT(1,IPEY)
                  IY2 = IY0 + IDLT(2,IPEY)
                  IY3 = IY0 + IDLT(3,IPEY)
                  IF ((IY1.LT.1) .OR. (IY2.LT.1) .OR. (IY3.LT.1) .OR.
     +                (IY1.GT.NY0) .OR. (IY2.GT.NY0) .OR.
     +                (IY3.GT.NY0)) GO TO 20
* Selects and/or supplements the y and z values.
                  Y1 = YD(IY1) - Y0
                  Z01 = ZD(IX0,IY1)
                  IF (NYD.GE.4) THEN
                      Y2 = YD(IY2) - Y0
                      Y3 = YD(IY3) - Y0
                      Z02 = ZD(IX0,IY2)
                      Z03 = ZD(IX0,IY3)
                  ELSE IF (NYD.EQ.3) THEN
                      Y2 = YD(IY2) - Y0
                      Z02 = ZD(IX0,IY2)
                      Y3 = 2*YD(3) - YD(2) - Y0
                      Z03 = Z3F(Y1,Y2,Y3,Z00,Z01,Z02)
                  ELSE IF (NYD.EQ.2) THEN
                      Y2 = 2*YD(2) - YD(1) - Y0
                      Z02 = Z2F(Y1,Y2,Z00,Z01)
                      Y3 = 2*YD(1) - YD(2) - Y0
                      Z03 = Z2F(Y1,Y3,Z00,Z01)
                  END IF
                  DZY01 = (Z01-Z00)/Y1
                  DZY02 = (Z02-Z00)/Y2
                  DZY03 = (Z03-Z00)/Y3
* Calculates the primary estimate of partial derivative zy as
* the coefficient of the bicubic polynomial.
                  CY1 = Y2*Y3/ ((Y1-Y2)* (Y1-Y3))
                  CY2 = Y3*Y1/ ((Y2-Y3)* (Y2-Y1))
                  CY3 = Y1*Y2/ ((Y3-Y1)* (Y3-Y2))
                  PEZY = CY1*DZY01 + CY2*DZY02 + CY3*DZY03
* Calculates the volatility factor and distance factor in the y
* direction for the primary estimate of zy.
                  SY = Y1 + Y2 + Y3
                  SZ = Z00 + Z01 + Z02 + Z03
                  SYY = Y1*Y1 + Y2*Y2 + Y3*Y3
                  SYZ = Y1*Z01 + Y2*Z02 + Y3*Z03
                  DNM = 4.0*SYY - SY*SY
                  B00 = (SYY*SZ-SY*SYZ)/DNM
                  B01 = (4.0*SYZ-SY*SZ)/DNM
                  DZ00 = Z00 - B00
                  DZ01 = Z01 - (B00+B01*Y1)
                  DZ02 = Z02 - (B00+B01*Y2)
                  DZ03 = Z03 - (B00+B01*Y3)
                  VOLF = DZ00**2 + DZ01**2 + DZ02**2 + DZ03**2
                  DISF = SYY
* Calculates the EPSLN value, which is used to decide whether or
* not the volatility factor is essentially zero.
                  EPSLN = (Z00**2+Z01**2+Z02**2+Z03**2)*1.0E-12
* Accumulates the weighted primary estimates of zy and their
* weights.
                  IF (VOLF.GT.EPSLN) THEN
* - For a finite weight.
                      WT = 1.0/ (VOLF*DISF)
                      SMPEF = SMPEF + WT*PEZY
                      SMWTF = SMWTF + WT
                  ELSE
* - For an infinite weight.
                      SMPEI = SMPEI + PEZY
                      SMWTI = SMWTI + 1.0
                  END IF
* Saves the necessary values for estimating zxy
                  YA(1,IPEY) = Y1
                  YA(2,IPEY) = Y2
                  YA(3,IPEY) = Y3
                  Z0IA(1,IPEY) = Z01
                  Z0IA(2,IPEY) = Z02
                  Z0IA(3,IPEY) = Z03
                  CYA(1,IPEY) = CY1
                  CYA(2,IPEY) = CY2
                  CYA(3,IPEY) = CY3
                  SYA(IPEY) = SY
                  SYYA(IPEY) = SYY
                  B00YA(IPEY) = B00
                  B01A(IPEY) = B01
   20         CONTINUE
* Calculates the final estimate of zy.
              IF (SMWTI.LT.0.5) THEN
* - When no infinite weights exist.
                  ZYDI = SMPEF/SMWTF
              ELSE
* - When infinite weights exist.
                  ZYDI = SMPEI/SMWTI
              END IF
* End of Part 2.
* Part 3.  Estimation of ZXYDI
* Initial setting
              SMPEF = 0.0
              SMWTF = 0.0
              SMPEI = 0.0
              SMWTI = 0.0
* Outer DO-loops with respect to the primary estimates in the x
* direction
              DO 40 IPEX = 1,4
                  IX1 = IX0 + IDLT(1,IPEX)
                  IX2 = IX0 + IDLT(2,IPEX)
                  IX3 = IX0 + IDLT(3,IPEX)
                  IF ((IX1.LT.1) .OR. (IX2.LT.1) .OR. (IX3.LT.1) .OR.
     +                (IX1.GT.NX0) .OR. (IX2.GT.NX0) .OR.
     +                (IX3.GT.NX0)) GO TO 40
* Retrieves the necessary values for estimating zxy in the x
* direction.
                  X1 = XA(1,IPEX)
                  X2 = XA(2,IPEX)
                  X3 = XA(3,IPEX)
                  Z10 = ZI0A(1,IPEX)
                  Z20 = ZI0A(2,IPEX)
                  Z30 = ZI0A(3,IPEX)
                  CX1 = CXA(1,IPEX)
                  CX2 = CXA(2,IPEX)
                  CX3 = CXA(3,IPEX)
                  SX = SXA(IPEX)
                  SXX = SXXA(IPEX)
                  B00X = B00XA(IPEX)
                  B10 = B10A(IPEX)
* Inner DO-loops with respect to the primary estimates in the y
* direction
                  DO 30 IPEY = 1,4
                      IY1 = IY0 + IDLT(1,IPEY)
                      IY2 = IY0 + IDLT(2,IPEY)
                      IY3 = IY0 + IDLT(3,IPEY)
                      IF ((IY1.LT.1) .OR. (IY2.LT.1) .OR.
     +                    (IY3.LT.1) .OR. (IY1.GT.NY0) .OR.
     +                    (IY2.GT.NY0) .OR. (IY3.GT.NY0)) GO TO 30
* Retrieves the necessary values for estimating zxy in the y
* direction.
                      Y1 = YA(1,IPEY)
                      Y2 = YA(2,IPEY)
                      Y3 = YA(3,IPEY)
                      Z01 = Z0IA(1,IPEY)
                      Z02 = Z0IA(2,IPEY)
                      Z03 = Z0IA(3,IPEY)
                      CY1 = CYA(1,IPEY)
                      CY2 = CYA(2,IPEY)
                      CY3 = CYA(3,IPEY)
                      SY = SYA(IPEY)
                      SYY = SYYA(IPEY)
                      B00Y = B00YA(IPEY)
                      B01 = B01A(IPEY)
* Selects and/or supplements the z values.
                      IF (NYD.GE.4) THEN
                          Z11 = ZD(IX1,IY1)
                          Z12 = ZD(IX1,IY2)
                          Z13 = ZD(IX1,IY3)
                          IF (NXD.GE.4) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z23 = ZD(IX2,IY3)
                              Z31 = ZD(IX3,IY1)
                              Z32 = ZD(IX3,IY2)
                              Z33 = ZD(IX3,IY3)
                          ELSE IF (NXD.EQ.3) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z23 = ZD(IX2,IY3)
                              Z31 = Z3F(X1,X2,X3,Z01,Z11,Z21)
                              Z32 = Z3F(X1,X2,X3,Z02,Z12,Z22)
                              Z33 = Z3F(X1,X2,X3,Z03,Z13,Z23)
                          ELSE IF (NXD.EQ.2) THEN
                              Z21 = Z2F(X1,X2,Z01,Z11)
                              Z22 = Z2F(X1,X2,Z02,Z12)
                              Z23 = Z2F(X1,X2,Z03,Z13)
                              Z31 = Z2F(X1,X3,Z01,Z11)
                              Z32 = Z2F(X1,X3,Z02,Z12)
                              Z33 = Z2F(X1,X3,Z03,Z13)
                          END IF
                      ELSE IF (NYD.EQ.3) THEN
                          Z11 = ZD(IX1,IY1)
                          Z12 = ZD(IX1,IY2)
                          Z13 = Z3F(Y1,Y2,Y3,Z10,Z11,Z12)
                          IF (NXD.GE.4) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z31 = ZD(IX3,IY1)
                              Z32 = ZD(IX3,IY2)
                          ELSE IF (NXD.EQ.3) THEN
                              Z21 = ZD(IX2,IY1)
                              Z22 = ZD(IX2,IY2)
                              Z31 = Z3F(X1,X2,X3,Z01,Z11,Z21)
                              Z32 = Z3F(X1,X2,X3,Z02,Z12,Z22)
                          ELSE IF (NXD.EQ.2) THEN
                              Z21 = Z2F(X1,X2,Z01,Z11)
                              Z22 = Z2F(X1,X2,Z02,Z12)
                              Z31 = Z2F(X1,X3,Z01,Z11)
                              Z32 = Z2F(X1,X3,Z02,Z12)
                          END IF
                          Z23 = Z3F(Y1,Y2,Y3,Z20,Z21,Z22)
                          Z33 = Z3F(Y1,Y2,Y3,Z30,Z31,Z32)
                      ELSE IF (NYD.EQ.2) THEN
                          Z11 = ZD(IX1,IY1)
                          Z12 = Z2F(Y1,Y2,Z10,Z11)
                          Z13 = Z2F(Y1,Y3,Z10,Z11)
                          IF (NXD.GE.4) THEN
                              Z21 = ZD(IX2,IY1)
                              Z31 = ZD(IX3,IY1)
                          ELSE IF (NXD.EQ.3) THEN
                              Z21 = ZD(IX2,IY1)
                              Z31 = Z3F(X1,X2,X3,Z01,Z11,Z21)
                          ELSE IF (NXD.EQ.2) THEN
                              Z21 = Z2F(X1,X2,Z01,Z11)
                              Z31 = Z2F(X1,X3,Z01,Z11)
                          END IF
                          Z22 = Z2F(Y1,Y2,Z20,Z21)
                          Z23 = Z2F(Y1,Y3,Z20,Z21)
                          Z32 = Z2F(Y1,Y2,Z30,Z31)
                          Z33 = Z2F(Y1,Y3,Z30,Z31)
                      END IF
* Calculates the primary estimate of partial derivative zxy as
* the coefficient of the bicubic polynomial.
                      DZXY11 = (Z11-Z10-Z01+Z00)/ (X1*Y1)
                      DZXY12 = (Z12-Z10-Z02+Z00)/ (X1*Y2)
                      DZXY13 = (Z13-Z10-Z03+Z00)/ (X1*Y3)
                      DZXY21 = (Z21-Z20-Z01+Z00)/ (X2*Y1)
                      DZXY22 = (Z22-Z20-Z02+Z00)/ (X2*Y2)
                      DZXY23 = (Z23-Z20-Z03+Z00)/ (X2*Y3)
                      DZXY31 = (Z31-Z30-Z01+Z00)/ (X3*Y1)
                      DZXY32 = (Z32-Z30-Z02+Z00)/ (X3*Y2)
                      DZXY33 = (Z33-Z30-Z03+Z00)/ (X3*Y3)
                      PEZXY = CX1* (CY1*DZXY11+CY2*DZXY12+CY3*DZXY13) +
     +                        CX2* (CY1*DZXY21+CY2*DZXY22+CY3*DZXY23) +
     +                        CX3* (CY1*DZXY31+CY2*DZXY32+CY3*DZXY33)
* Calculates the volatility factor and distance factor in the x
* and y directions for the primary estimate of zxy.
                      B00 = (B00X+B00Y)/2.0
                      SXY = SX*SY
                      SXXY = SXX*SY
                      SXYY = SX*SYY
                      SXXYY = SXX*SYY
                      SXYZ = X1* (Y1*Z11+Y2*Z12+Y3*Z13) +
     +                       X2* (Y1*Z21+Y2*Z22+Y3*Z23) +
     +                       X3* (Y1*Z31+Y2*Z32+Y3*Z33)
                      B11 = (SXYZ-B00*SXY-B10*SXXY-B01*SXYY)/SXXYY
                      DZ00 = Z00 - B00
                      DZ01 = Z01 - (B00+B01*Y1)
                      DZ02 = Z02 - (B00+B01*Y2)
                      DZ03 = Z03 - (B00+B01*Y3)
                      DZ10 = Z10 - (B00+B10*X1)
                      DZ11 = Z11 - (B00+B01*Y1+X1* (B10+B11*Y1))
                      DZ12 = Z12 - (B00+B01*Y2+X1* (B10+B11*Y2))
                      DZ13 = Z13 - (B00+B01*Y3+X1* (B10+B11*Y3))
                      DZ20 = Z20 - (B00+B10*X2)
                      DZ21 = Z21 - (B00+B01*Y1+X2* (B10+B11*Y1))
                      DZ22 = Z22 - (B00+B01*Y2+X2* (B10+B11*Y2))
                      DZ23 = Z23 - (B00+B01*Y3+X2* (B10+B11*Y3))
                      DZ30 = Z30 - (B00+B10*X3)
                      DZ31 = Z31 - (B00+B01*Y1+X3* (B10+B11*Y1))
                      DZ32 = Z32 - (B00+B01*Y2+X3* (B10+B11*Y2))
                      DZ33 = Z33 - (B00+B01*Y3+X3* (B10+B11*Y3))
                      VOLF = DZ00**2 + DZ01**2 + DZ02**2 + DZ03**2 +
     +                       DZ10**2 + DZ11**2 + DZ12**2 + DZ13**2 +
     +                       DZ20**2 + DZ21**2 + DZ22**2 + DZ23**2 +
     +                       DZ30**2 + DZ31**2 + DZ32**2 + DZ33**2
                      DISF = SXX*SYY
* Calculates EPSLN.
                      EPSLN = (Z00**2+Z01**2+Z02**2+Z03**2+Z10**2+
     +                        Z11**2+Z12**2+Z13**2+Z20**2+Z21**2+Z22**2+
     +                        Z23**2+Z30**2+Z31**2+Z32**2+Z33**2)*
     +                        1.0E-12
* Accumulates the weighted primary estimates of zxy and their
* weights.
                      IF (VOLF.GT.EPSLN) THEN
* - For a finite weight.
                          WT = 1.0/ (VOLF*DISF)
                          SMPEF = SMPEF + WT*PEZXY
                          SMWTF = SMWTF + WT
                      ELSE
* - For an infinite weight.
                          SMPEI = SMPEI + PEZXY
                          SMWTI = SMWTI + 1.0
                      END IF
   30             CONTINUE
   40         CONTINUE
* Calculates the final estimate of zxy.
              IF (SMWTI.LT.0.5) THEN
* - When no infinite weights exist.
                  ZXYDI = SMPEF/SMWTF
              ELSE
* - When infinite weights exist.
                  ZXYDI = SMPEI/SMWTI
              END IF
* End of Part 3
              PDD(1,IX0,IY0) = ZXDI
              PDD(2,IX0,IY0) = ZYDI
              PDD(3,IX0,IY0) = ZXYDI
   50     CONTINUE
   60 CONTINUE

      RETURN
      END

C ############################################################################
C ############################################################################

C Section 6:

C Forward WATER-model routines

C Limits at 443nm:
C                   Chl:  0.02   - 10.0 [mg/m3]
C                   acdm: 0.003  - 0.3  [1/m]
C                   bbp:  0.0003 - 0.01 [1/m]

C *****************************************************************************

C called by SPEC1GEOFIT case 2 loop

        SUBROUTINE GS97_CASE2 (C,acdm0,bbp0,LwN)

        INTEGER    NLAMBDA
        PARAMETER (NLAMBDA=8)
        real LwN(NLAMBDA),lambda(NLAMBDA), aw(NLAMBDA),
     &       bbw(NLAMBDA), aph_lin(NLAMBDA), Fobar(NLAMBDA)
        REAL*4 adm_s,bbp_s

        INTEGER i
        REAL*4  abs_coef,admstar,bb,bbpstar,
     &          grd1,grd2,x,F

        COMMON /WATERPAR/ adm_s,bbp_s,aw,bbw,aph_lin

        data lambda /412.0, 443.0, 490.0, 510.0,
     &               555.0, 670.0, 765.0, 865.0 /

c       data Fobar /170.79,189.45,193.66,188.35,
c     &              185.33,153.,122.24, 98.82/

        DATA Fobar   /173.004, 190.154, 196.473, 188.158,
     &                183.01,  151.143, 122.316, 96.302  /

C       Gordon constants
        grd1 = 0.0949 
        grd2 = 0.0794                                                   

        do i = 1, NLAMBDA

C       acdm spectral dependence
C       ************************
        admstar = exp( - adm_s * (lambda(i) - 443.0))

C       total absorption coefficient ( = aw + aph + acdm)
C       *************************************************
        abs_coef = aw(i) + aph_lin(i)*C + acdm0*admstar

C       bbp spectral dependence
C       ************************
        bbpstar = (443.0/lambda(i)) ** bbp_s

C       total backscattering coefficient ( = bbw + bbp)
C       ***********************************************
        bb = bbw(i) + bbp0*bbpstar

C       bb/(bb+a) ratio
C       ****************
        x = bb/(abs_coef + bb)

C       R/Q (= ~Lu[0-]/Ed[0-]) as in Gordon et al., 1988
C       *************************************************
        F = grd1*x + grd2*x**2          

C       Conversion to LwN (assumes flat surface, high sza)
C       ***************************************************
        LwN(i) = F*0.54*0.97*Fobar(i)
        
        enddo

c        LwN(7) = 0.0
c        LwN(8) = 0.0

        return
        end

C *****************************************************************************

C called by DFCNLSQ_GRAD

        SUBROUTINE DGS97(C, acdm0, bbp0, LwN)

C       Forward model (Garver & Siegel, 1997)
C       The model uses the quadratic formulation R/Q = f(bb/a+bb)
C       from Gordon et al., 1988
C       It has been optimized through simulated annealing using a 
C       1075 stations Rrs/Kd/Chl data set
C
C       Inputs :
C       C       : [Chl]
C       acdm0   : acdm (i.e. ay + ad) at a reference wavelength (443 nm)
C       bbp0    : particulate backscattering coeff. at a reference wavelength 
C                 (443 nm)
C
C       Output : LwN at the SeaWIFS wavebands

        IMPLICIT NONE
        INTEGER    NLAMBDA
        PARAMETER (NLAMBDA=8)
        REAL*8  C, acdm0, bbp0
        real*8  LwN(NLAMBDA),lambda(NLAMBDA),Fobar(NLAMBDA)

        INTEGER i
        REAL*8  abs_coef,admstar,bb,bbpstar,
     &          grd1,grd2,x,F

        REAL*8 d_adm_s,d_bbp_s,
     &         d_aw(NLAMBDA),d_bbw(NLAMBDA),d_aph_lin(NLAMBDA)
        COMMON /D_WATERPAR/ d_adm_s,d_bbp_s,d_aw,d_bbw,d_aph_lin

        data lambda /412.D0, 443.D0, 490.D0, 510.D0,
     &               555.D0, 670.D0, 765.D0, 865.D0 /

c       data aw /0.00456,0.00707,0.015,0.0325,0.0596,0.439,0.,0./
c       data bbw /0.00333,0.00239,0.00155,0.0013,0.000925,0.0004,0.,0./
c       data aph_lin /0.00665, 0.05582, 0.02055, 0.01910, 
c     &                0.01015, 0.01424, 0.0, 0.0/

c       data Fobar /170.79D0,189.45D0,193.66D0,188.35D0,
c     &              185.33D0,153.D0,0.D0,0.D0/

        DATA Fobar   /173.004D0, 190.154D0, 196.473D0, 188.158D0,
     &                183.01D0,  151.143D0, 0.D0, 0.D0  /

C       slope of the CDM ( i.e. ay + ad) absorption
C       adm_s = 0.0206
C       exponent of particulate backscattering
C       bbp_s = 1.03

C       Gordon constants
        grd1 = 0.0949D0 
        grd2 = 0.0794D0                                                 

        do i = 1, 6

C       acdm spectral dependence
C       ************************
        admstar = dexp( - d_adm_s * (lambda(i) - 443.D0))

C       total absorption coefficient ( = aw + aph + acdm)
C       *************************************************
        abs_coef = d_aw(i) + d_aph_lin(i)*C + acdm0*admstar

C       bbp spectral dependence
C       ************************        
        bbpstar = (443.D0/lambda(i)) ** d_bbp_s

C       total backscattering coefficient ( = bbw + bbp)
C       ***********************************************
        bb = d_bbw(i) + bbp0*bbpstar

C       bb/(bb+a) ratio
C       ****************
        x = bb/(abs_coef + bb)

C       R/Q (= ~Lu[0-]/Ed[0-]) as in Gordon et al., 1988
C       *************************************************
        F = grd1*x + grd2*x**2          

C       Conversion to LwN (assumes flat surface, high sza)
C       ***************************************************
        LwN(i) = F*0.54D0*0.97D0*Fobar(i)
        
        enddo

        LwN(7) = 0.D0
        LwN(8) = 0.D0

        return
        end

C *****************************************************************************

C called by SPECOPT_POSTPROC

        SUBROUTINE DGS97_CASE2 (C,acdm0,bbp0,LwN,absCoefnw)

        INTEGER    NLAMBDA
        PARAMETER (NLAMBDA=8)
        REAL*8 C,acdm0,bbp0
        REAL*8 LwN(NLAMBDA),lambda(NLAMBDA), Fobar(NLAMBDA)

        DATA lambda /412.D0, 443.D0, 490.D0, 510.D0,
     &               555.D0, 670.D0, 765.D0, 865.D0 /

        REAL*8  absCoefnw(NLAMBDA)

        INTEGER i
        REAL*8  abs_coef,admstar,bb,bbpstar,
     &          grd1,grd2,x,F

        REAL*8 d_adm_s,d_bbp_s,
     &         d_aw(NLAMBDA),d_bbw(NLAMBDA),d_aph_lin(NLAMBDA)
        COMMON /D_WATERPAR/ d_adm_s,d_bbp_s,d_aw,d_bbw,d_aph_lin

c       DATA Fobar /170.79D0,189.45D0,193.66D0,188.35D0,
c     &              185.33D0,153.D0,122.24D0, 98.82D0/

        DATA Fobar   /173.004D0, 190.154D0, 196.473D0, 188.158D0,
     &                183.01D0,  151.143D0, 122.316D0, 96.302D0  /

C       Gordon constants
        grd1 = 0.0949D0
        grd2 = 0.0794D0

        do i = 1, NLAMBDA

C       acdm spectral dependence
C       ************************
        admstar = dexp( - d_adm_s * (lambda(i) - 443.D0))

C       total absorption coefficient ( = aw + aph + acdm)
C       *************************************************
        abs_coef = d_aw(i) + d_aph_lin(i)*C + acdm0*admstar

C       non-water absorption coefficient (= aph + acdm)
C       *************************************************
        absCoefnw(i) = d_aph_lin(i)*C + acdm0*admstar

C       bbp spectral dependence
C       ************************
        bbpstar = (443.D0/lambda(i)) ** d_bbp_s

C       total backscattering coefficient ( = bbw + bbp)
C       ***********************************************
        bb = d_bbw(i) + bbp0*bbpstar

C       bb/(bb+a) ratio
C       ****************
        x = bb/(abs_coef + bb)

C       R/Q (= ~Lu[0-]/Ed[0-]) as in Gordon et al., 1988
C       *************************************************
        F = grd1*x + grd2*x**2          

C       Conversion to LwN (assumes flat surface, high sza)
C       ***************************************************
        LwN(i) = F*0.54D0*0.97D0*Fobar(i)
        
        enddo

c       LwN(7) = 0.D0
c       LwN(8) = 0.D0

        return
        end

C *****************************************************************************

        subroutine get_chl_oc4v4(nLw,chl)

        INTEGER    NLAMBDA, len
        PARAMETER (NLAMBDA=8, len=5)
        
        INTEGER i443,i490,i510,i555
        PARAMETER (i443=2,i490=3,i510=4,i555=5)

        REAL nLw(NLAMBDA)
        REAL Fobar(NLAMBDA)
        REAL a(len)

        REAL Rrs443,Rrs490,Rrs510,Rrs555
        REAL minRrs

        REAL ratio,rat
  
        REAL chl

c        DATA Fobar /170.79,189.45,193.66,188.35,
c     &              185.33,153.41,122.24, 98.82/

c       Thuillier band-averaged (square block) values
        DATA Fobar /173.004,190.154,196.473,188.158,
     &              183.010,151.143,122.316, 96.302/

        DATA a /0.366,-3.067,1.930,0.649,-1.532/

        Rrs443 = nLw(i443)/Fobar(i443)
        Rrs490 = nLw(i490)/Fobar(i490)
        Rrs510 = nLw(i510)/Fobar(i510)
        Rrs555 = nLw(i555)/Fobar(i555)

c        print*, Rrs443/Rrs555
c        print*, Rrs490/Rrs555
c        print*, Rrs510/Rrs555
c        print*, Rrs555/Rrs555



        minRrs = MIN(Rrs443,Rrs490)

c       We require Rrs555 to be positive, and we require that if any band
c       goes negative, it must occur in order of wavelength
c       note: this is Franz 'standard algorith' requirement.  Not relavent to SOA
c       but mentioned for clarity

        IF (Rrs555 .gt. 0.0 .and. Rrs510 .gt. 0.0 .and. 
     &           (Rrs490 .gt. 0.0 .or. Rrs443*Rrs490 .gt. 0.0) .and.
     &                                 minRrs .gt. -0.001) THEN

          ratio = MAX(MAX(Rrs443,Rrs490),Rrs510)/Rrs555

c         Fail if ratio is unphysical (Rat=0.21 -> Chl=640)
          IF (ratio .gt. 0.21) THEN
            rat = log10(ratio)
            chl = 10.0**(a(1)+rat*(a(2)+rat*(a(3)+rat*(a(4)+rat*a(5)))))
          ELSE
            chl = -1.   !dummy
          ENDIF

        ELSE
          chl = -2.     !dummy
        ENDIF
        
        return
        end

C *****************************************************************************

      REAL FUNCTION OC4_CHL4S(rhow)

      INTEGER NLAMBDA
      PARAMETER (NLAMBDA=8)

      REAL rhow(NLAMBDA)

      INTEGER i443,i490,i510,i555
      PARAMETER (i443=2,i490=3,i510=4,i555=5)
      REAL r4s,rs443,rs490,rs510

      rs443 = rhow(i443)/rhow(i555)
      rs490 = rhow(i490)/rhow(i555)
      rs510 = rhow(i510)/rhow(i555)

      IF(rs443.GT.rs490) THEN
          r4s = LOG10(rs443)      
      ELSE IF(rs490.GT.rs510) THEN
          r4s = LOG10(rs490)      
      ELSE
          r4s = LOG10(rs510)      
      ENDIF

      OC4_CHL4S = 10.0**(0.366 - 3.067*r4s + 1.93*r4s**2 +
     &                  0.649*r4s**3 - 1.532*r4s**4)

      RETURN
      END


C ############################################################################
C ############################################################################

C SECTION 7:

C Atmosphere LUT and and water initialisation data

C Subroutine READ_OCEAN_DATA not required now.  LUT's for SeaWiFS water reflectance
C contain x=b0*10/3 where 0.4 < x < 1.5.  bbp is now calculated in optimization
C and hence water LUT's not required.

C *****************************************************************************

      SUBROUTINE READ_MODEL_DATA (sensorNm)

c for each model jmodel [mr,mi,nu] this subroutine reads in the single scattering albedo,
c extinction coefficient and the 75 elements (scattering directions) of the Scat Muller
c matrix, all at each wavelength (within each model 'jmodel').
c ccTab865 is then calculated as the ratio of TAU(lambda)/TAU(865nm) = C(lambda)/C(865nm)
c for each wavelength (within each model 'jmodel').

c LUT values of the rhoA v Tau fourier coefficients a,b,c,d are also read.
c They have been computed as a quartic fit for each of 72 models (mi,mr,v) 
c at each sun angle, lambda and phi (33*8*15*72 = 285 120 values).
c They are stored in luttau_a/b/c/d :
c rhoa (params) = a(params)tau + b(params)tau^2 + c(params)tau^3 + d(params)tau^4
c They are required by subroutine funct_a_b_c_d_fourier_interpolate.

c Finally, DIFFUSE TANSMITTANCE COEFFICIENTS are stored in LUT's (diffTrans_a / b).

      INTEGER    NRAD,MPHI,NUM
      PARAMETER (NRAD=100,MPHI=15,NUM=75)

C -- INPUT/PARAMETERS=SPECMATCH_GENERIC
      INTEGER    NMI,NMR,NV,
     &           NLAMBDA,NLAMBD0
C -----------------------------------------------------------------------------
C File: SPECMATCH_GENERIC.PAR
C -----------------------------------------------------------------------------
      PARAMETER (NMI=6,NMR=2,NV=6)
      PARAMETER (NLAMBDA=8,NLAMBD0=2)
C -----------------------------------------------------------------------------
C -- SEADAS-SPECIFIC
      INTEGER    NDIR,NDIR_OUT,NMODEL,NSUN
      PARAMETER (NDIR=35,NDIR_OUT=26,NMODEL=NV*NMR*NMI,NSUN=33)

C -- INTERNAL/MISC
      INTEGER    jmodel,jl,jt,jord,jdir
C
C -- PHASE FUNCTION DATA ----------------------------------------------------
C
C       aerW0     Single Scat Albedo W0
C       aerC      Extinction Coeff   C
C       s11       S11 elements of Scat Muller Matrix
C       anglPhs   scat angles (not required)
C
      REAL    aerW0(NLAMBDA,NMI,NMR,NV), aerC(NLAMBDA,NMI,NMR,NV),
     &        cc865tab(NLAMBDA,NMI,NMR,NV),
     &        s11dummy(NUM)
C
C -- Haze-C MODEL PARAMETER DATA --------------------------------------------
C
C --> MODEL
C#1
        INTEGER*4 iv(NV), jv
        DATA iv
     &            /20, 25, 30, 35, 40, 45/    ! actually 2.0,2.5,...,4.5
C#2
        INTEGER*4 imr(NMR), jr
        DATA imr
     &            /133, 150/                  ! actually 1.33, 1.50
C#3
        INTEGER*4 imi(NMI), ji
        DATA imi
     &            /0, 1, 3, 10, 30, 40/       ! actually 0.000,0.001,...,0.040
C
C    Look-Up Reflectance LUT
      REAL    luttau_a(NDIR,MPHI,NSUN,NLAMBDA,NMODEL),
     &        luttau_b(NDIR,MPHI,NSUN,NLAMBDA,NMODEL),
     &        luttau_c(NDIR,MPHI,NSUN,NLAMBDA,NMODEL),
     &        luttau_d(NDIR,MPHI,NSUN,NLAMBDA,NMODEL)

C    Look-Up Diffuse Transmittance LUT
      REAL    diffTrans_a(NSUN,NLAMBDA,NMODEL),
     &        diffTrans_b(NSUN,NLAMBDA,NMODEL)
C
CBAF
      CHARACTER sensorNm*20
      INTEGER   flen
      INTEGER   slen
      CHARACTER filedir*255
      CHARACTER cSEADASlut*255

C -- COMMON
      COMMON /COMPHASE/    aerW0,aerC,cc865tab,s11dummy
      COMMON /COMRTELUT/   luttau_a,luttau_b,luttau_c,luttau_d
      COMMON /COMDIFF/     diffTrans_a,diffTrans_b
cc      SAVE   /COMPHASE/,/COMRTELUT/,/COMDIFF/

C:START$RUN
CBAF

      call getenv('OCDATAROOT',filedir)
      if (lenstr(filedir) .eq. 0) then
          write(*,*)
     &    '-E- : Environment variable OCDATAROOT undefined'
          call exit(1)
      endif
      flen = lenstr(filedir)
      slen = lenstr(sensorNm)

      DO jmodel=1,NMODEL
C ---------------------- LOOKING FOR Parameters_Index ACCORDING TO MODEL INDEX
        CALL IDX2MODEL(jmodel,NV,NMR,NMI,jv,jr,ji)

C ---------------------- MAKING OUTPUT FILE_NAME
        WRITE(unit=cSEADASlut,fmt=10) iv(jv),imr(jr),imi(ji)
10      FORMAT(
     &    'hzc_v',i2.2,'_r',i3.3,'_i',i3.3,'_z20_s00.dat')

        cSEADASlut = filedir(1:flen)//'/'//sensorNm(1:slen)//
     &        '/aerosol/'//sensorNm(1:slen)//'_'//cSEADASlut  
        WRITE(*,*) cSEADASlut

CBAF        WRITE(unit=cSEADASlut,fmt=10) iv(jv),imr(jr),imi(ji)
CBAF10      FORMAT(
CBAF     &    '/home/kuchinke/MySOA/SeaWiFS_HazeC_LUT_orig/seawifs_',
CBAF     &       'hzc_v',i2.2,'_r',i3.3,'_i',i3.3,'_z20_s00.dat')

C ---------------------- READING DATA
        OPEN(unit=30, file=cSEADASlut, status='UNKNOWN')

C -- PHASE FUNCTION
c read in ssalbedo, extinction coefficient and elements of Scat Muller matrix
        DO jl=1,NLAMBDA
          READ(30,21)
          READ(30, * ) aerW0(jl,ji,jr,jv), aerC(jl,ji,jr,jv)
          READ(30,21)
          READ(30,20) (s11dummy(jj), jj = 1,NUM)
        ENDDO ! jl

C -- TAU(lambda)/TAU(865nm) = C(lambda)/C(865nm)
        DO jl=1,NLAMBDA

         IF(jl.EQ.NLAMBDA) THEN
          cc865tab(jl,ji,jr,jv)=1.0   !ratio is 1.0 at 865nm
         ELSE
          ! exact values of ratio at other wavelengths
          cc865tab(jl,ji,jr,jv)=aerC(jl,ji,jr,jv)/aerC(NLAMBDA,ji,jr,jv)
         ENDIF
        ENDDO ! jl

C -- RTELUT
c LUT values of Tau a,b,c,d coefficients for each model (mi,mr,v) at each
c sun angle, lambda and phi (33*8*15) = 285 120 values.  Stored in luttau_a/b/c/d
        DO jt=1,NSUN
          DO jl=1,NLAMBDA
           DO jord=1,MPHI
            READ(30,21)
            READ(30,20) (luttau_a(jdir,jord,jt,jl,jmodel), 
     &                   jdir=1,NDIR_OUT) 
            READ(30,21)
            READ(30,20) (luttau_b(jdir,jord,jt,jl,jmodel), 
     &                   jdir=1,NDIR_OUT) 
            READ(30,21)
            READ(30,20) (luttau_c(jdir,jord,jt,jl,jmodel), 
     &                   jdir=1,NDIR_OUT) 
            READ(30,21)
            READ(30,20) (luttau_d(jdir,jord,jt,jl,jmodel), 
     &                   jdir=1,NDIR_OUT) 
           ENDDO  ! MPHI
          ENDDO  ! NLAMBDA

C -- DIFFUSE TANSMITTANCE COEFFICIENTS
c LUT values stored in diffTrans_a / b
          READ(30,21)
          DO jl = 1, nlambda    ! wavelength
            READ(30,*) diffTrans_a(jt,jl,jmodel), 
     &                 diffTrans_b(jt,jl,jmodel) 
          ENDDO    ! NLAMBDA
        ENDDO    ! NSUN

        CLOSE(30)

      ENDDO ! NMODEL

20    FORMAT( 5(2x,E13.6) )
21    FORMAT( )

      RETURN
      END

C *****************************************************************************
c GSM01 coeffs -various

      SUBROUTINE INIT_WATER_PAR1()

      IMPLICIT NONE

      INTEGER NLAMBDA
      PARAMETER (NLAMBDA=8)

      INTEGER jl

      INTEGER iCase1,iSantaBarb,iMaineOct,iMaineDec,iTanaka,
     & iToratani,iAndreaINsh,  iAndreaMLsm, iAndreaUPfa
      PARAMETER (iCase1=1,iSantaBarb=2,iMaineOct=3,iMaineDec=4,
     & iTanaka=5,iToratani=6,iAndreaINsh=7,iAndreaMLsm=8,
     & iAndreaUPfa=9)
      INTEGER iWhichWater

      REAL*4 adm_s,bbp_s,aw(NLAMBDA),bbw(NLAMBDA),aph_lin(NLAMBDA)
      COMMON /WATERPAR/ adm_s,bbp_s,aw,bbw,aph_lin

      REAL*8 d_adm_s,d_bbp_s,
     &       d_aw(NLAMBDA),d_bbw(NLAMBDA),d_aph_lin(NLAMBDA)
      COMMON /D_WATERPAR/ d_adm_s,d_bbp_s,d_aw,d_bbw,d_aph_lin

      iWhichWater = iCase1
c      iWhichWater = iSantaBarb
c      iWhichWater = iMaineOct
c      iWhichWater = iTanaka
c      iWhichWater = iToratani
c      iWhichWater = iAndreaINsh
c      iWhichWater = iAndreaMLsm
c      iWhichWater = iAndreaUPfa   

      IF (iWhichWater.EQ.iCase1) THEN

           adm_s = 0.0206
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

c           aph_lin(1) = 0.00665   !original
           aph_lin(1) = 0.02       !increased
           aph_lin(2) = 0.05582
           aph_lin(3) = 0.02055
           aph_lin(4) = 0.0191
           aph_lin(5) = 0.01015
           aph_lin(6) = 0.01424
           aph_lin(7) = 0.0
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iSantaBarb) THEN 

           adm_s = 0.0190473077  
           bbp_s = 1.2291496351 

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.0409560862
           aph_lin(2) = 0.0555838757
           aph_lin(3) = 0.0184040713
           aph_lin(4) = 0.0234303520
           aph_lin(5) = 0.0078382969
           aph_lin(6) = 0.02024         
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iMaineOct) THEN

           adm_s = 0.01821
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.038238392   !actual
           aph_lin(2) = 0.047316183
           aph_lin(3) = 0.031646183
           aph_lin(4) = 0.02057424
           aph_lin(5) = 0.0082816764
           aph_lin(6) = 0.019338809
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iTanaka) THEN

           adm_s = 0.014
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.03880149
           aph_lin(2) = 0.0472593
           aph_lin(3) = 0.02982772
           aph_lin(4) = 0.01772824
           aph_lin(5) = 0.0057802
           aph_lin(6) = 0.01898189
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iToratani) THEN

           adm_s = 0.014
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.0559
           aph_lin(2) = 0.0628
           aph_lin(3) = 0.0364
           aph_lin(4) = 0.0175
           aph_lin(5) = 0.0075
           aph_lin(6) = 0.0245
           aph_lin(7) = 0.0
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iAndreaINsh) THEN

           adm_s = 0.014
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.07952
           aph_lin(2) = 0.09755
           aph_lin(3) = 0.06507
           aph_lin(4) = 0.04336
           aph_lin(5) = 0.01814
           aph_lin(6) = 0.03885
           aph_lin(7) = 0.0
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iAndreaMLsm) THEN

           adm_s = 0.014
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.03349
           aph_lin(2) = 0.03883
           aph_lin(3) = 0.02305
           aph_lin(4) = 0.01661
           aph_lin(5) = 0.00614
           aph_lin(6) = 0.02277
           aph_lin(7) = 0.0
           aph_lin(8) = 0.0


      ELSE IF (iWhichWater.EQ.iAndreaUPfa) THEN

           adm_s = 0.014
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.05206
           aph_lin(2) = 0.04812
           aph_lin(3) = 0.02835
           aph_lin(4) = 0.02287
           aph_lin(5) = 0.01134
           aph_lin(6) = 0.02318
           aph_lin(7) = 0.0
           aph_lin(8) = 0.0

      ENDIF

      d_adm_s = DBLE(adm_s)
      d_bbp_s = DBLE(bbp_s)
      DO jl=1,NLAMBDA
          d_aw(jl)      = DBLE(aw(jl))
          d_bbw(jl)     = DBLE(bbw(jl))
          d_aph_lin(jl) = DBLE(aph_lin(jl))
      ENDDO

      RETURN
      END

C *****************************************************************************
c GSM01 coeffs -Chesapeake Bay (Magnuson, Harding, Mallonee, Adolf 2004)

      SUBROUTINE INIT_WATER_PAR2(lat,lon,iseason)

      IMPLICIT NONE

      INTEGER NLAMBDA
      PARAMETER (NLAMBDA=8)

      INTEGER jl

      INTEGER iupperBay_spring,iupperBay_summer,iupperBay_fall,
     &        imidBay_spring,imidBay_summer,imidBay_fall,
     &        ilowerBay_spring,ilowerBay_summer,ilowerBay_fall,
     &        inshore,ioffshore
      PARAMETER (iupperBay_spring=1,iupperBay_summer=2,iupperBay_fall=3,
     &        imidBay_spring=4,imidBay_summer=5,imidBay_fall=6,
     &        ilowerBay_spring=7,ilowerBay_summer=8,ilowerBay_fall=9,
     &        inshore=10,
     &        ioffshore=11)
      INTEGER iWhichWater

      REAL*4 adm_s,bbp_s,aw(NLAMBDA),bbw(NLAMBDA),aph_lin(NLAMBDA)
      COMMON /WATERPAR/ adm_s,bbp_s,aw,bbw,aph_lin

      REAL*8 d_adm_s,d_bbp_s,
     &       d_aw(NLAMBDA),d_bbw(NLAMBDA),d_aph_lin(NLAMBDA)
      COMMON /D_WATERPAR/ d_adm_s,d_bbp_s,d_aw,d_bbw,d_aph_lin

      REAL*4    lat,lon
      INTEGER   iseason

      INTEGER   ispring,isummer,ifall
      PARAMETER (ispring=1,isummer=2,ifall=3)
      INTEGER   iwhichseason
      
      iwhichseason = iseason

      IF (lat .gt. 38.60 .and. lat .lt. 39.56 
     &                .and. lon .gt. -77.00  .and. lon .lt. -75.70) THEN

        IF (iwhichseason .eq. ispring) THEN
          iWhichWater = iupperBay_spring
        ELSE IF (iwhichseason .eq. isummer) THEN
          iWhichWater = iupperBay_summer
        ELSE IF (iwhichseason .eq. ifall) THEN
          iWhichWater = iupperBay_fall
        ENDIF


      ELSE IF (lat .gt. 37.60 .and. lat .le. 38.60 
     &                .and. lon .gt. -77.00  .and. lon .lt. -75.50) THEN

        IF (iwhichseason .eq. ispring) THEN
          iWhichWater = imidBay_spring
        ELSE IF (iwhichseason .eq. isummer) THEN
          iWhichWater = imidBay_summer
        ELSE IF (iwhichseason .eq. ifall) THEN
          iWhichWater = imidBay_fall
        ENDIF

      ELSE IF (lat .gt. 37.11 .and. lat .le. 37.60 
     &                .and. lon .gt. -77.00  .and. lon .lt. -75.85) THEN

        IF (iwhichseason .eq. ispring) THEN
          iWhichWater = ilowerBay_spring
        ELSE IF (iwhichseason .eq. isummer) THEN
          iWhichWater = ilowerBay_summer
        ELSE IF (iwhichseason .eq. ifall) THEN
          iWhichWater = ilowerBay_fall
        ENDIF

      ELSE IF (lat .ge. 36.90 .and. lat .le. 37.11
     &             .and. lon .gt. -77.00 .and. lon .lt. -75.9) THEN

        iWhichWater = inshore

      ELSE

        iWhichWater = ioffshore

      ENDIF

      IF (iWhichWater.EQ. iupperBay_spring) THEN

           adm_s = 0.01218
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.05960
           aw(6) = 0.43900
           aw(7) = 2.570
           aw(8) = 4.4400

           bbw(1) = 0.003330
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.00130
           bbw(5) = 0.000925
           bbw(6) = 0.00040
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.02119
           aph_lin(2) = 0.02509
           aph_lin(3) = 0.01282
           aph_lin(4) = 0.00919 
           aph_lin(5) = 0.00427
           aph_lin(6) = 0.02087
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iupperBay_summer) THEN

           adm_s = 0.01218 
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.05960
           aw(6) = 0.43900
           aw(7) = 2.570
           aw(8) = 4.4400

           bbw(1) = 0.003330
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.00130
           bbw(5) = 0.000925
           bbw(6) = 0.00040
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.02653
           aph_lin(2) = 0.02979
           aph_lin(3) = 0.01655
           aph_lin(4) = 0.01208 
           aph_lin(5) = 0.00470
           aph_lin(6) = 0.02122
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.iupperBay_fall) THEN 

           adm_s = 0.01218
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.02259
           aph_lin(2) = 0.02573
           aph_lin(3) = 0.01372
           aph_lin(4) = 0.01019
           aph_lin(5) = 0.00376
           aph_lin(6) = 0.02066        
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.imidBay_spring) THEN

           adm_s = 0.01385
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.02001
           aph_lin(2) = 0.02212
           aph_lin(3) = 0.01279
           aph_lin(4) = 0.00974
           aph_lin(5) = 0.00449
           aph_lin(6) = 0.01588
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.imidBay_summer) THEN

           adm_s = 0.01385  
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.03345 
           aph_lin(2) = 0.039
           aph_lin(3) = 0.02318
           aph_lin(4) = 0.01664
           aph_lin(5) = 0.00609
           aph_lin(6) = 0.02285
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0


      ELSE IF (iWhichWater.EQ.imidBay_fall) THEN

           adm_s = 0.01385
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.02758
           aph_lin(2) = 0.0308
           aph_lin(3) = 0.01826
           aph_lin(4) = 0.01438
           aph_lin(5) = 0.00691
           aph_lin(6) = 0.02112
           aph_lin(7) = 0.0
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.ilowerBay_spring) THEN

            adm_s = 0.0133
            bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.02001
           aph_lin(2) = 0.02212
           aph_lin(3) = 0.01279
           aph_lin(4) = 0.00974
           aph_lin(5) = 0.00449
           aph_lin(6) = 0.01588
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0


      ELSE IF (iWhichWater.EQ.ilowerBay_summer) THEN

           adm_s = 0.0133
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.57
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.03345
           aph_lin(2) = 0.039
           aph_lin(3) = 0.02318
           aph_lin(4) = 0.01664
           aph_lin(5) = 0.00609
           aph_lin(6) = 0.02285
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.ilowerBay_fall) THEN

           adm_s = 0.0133
           bbp_s = 1.1

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.570
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

           aph_lin(1) = 0.02758
           aph_lin(2) = 0.0308
           aph_lin(3) = 0.01826
           aph_lin(4) = 0.01438
           aph_lin(5) = 0.00691
           aph_lin(6) = 0.02112
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0

      ELSE IF (iWhichWater.EQ.inshore) THEN

           adm_s = 0.01236
           bbp_s = 1.1
 
           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.570
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016
      
           aph_lin(1) = 0.07123
           aph_lin(2) = 0.08843
           aph_lin(3) = 0.06024
           aph_lin(4) = 0.04072
           aph_lin(5) = 0.01693
           aph_lin(6) = 0.03815
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0
       
      ELSE IF (iWhichWater.EQ.ioffshore) THEN

c           adm_s = 0.01646    !Magnuson et al. 2005
           adm_s = 0.0185      !for testing

c          Case 1
c           adm_s = 0.0206
           bbp_s = 1.03

           aw(1) = 0.00456
           aw(2) = 0.00707
           aw(3) = 0.015
           aw(4) = 0.0325
           aw(5) = 0.0596
           aw(6) = 0.439
           aw(7) = 2.570
           aw(8) = 4.44

           bbw(1) = 0.00333
           bbw(2) = 0.00239
           bbw(3) = 0.00155
           bbw(4) = 0.0013
           bbw(5) = 0.000925
           bbw(6) = 0.0004
           bbw(7) = 0.00026
           bbw(8) = 0.00016

c          Magnuson et al. 2005
           aph_lin(1) = 0.11331
           aph_lin(2) = 0.14678
           aph_lin(3) = 0.09832
           aph_lin(4) = 0.06048
           aph_lin(5) = 0.0192
           aph_lin(6) = 0.04349
           aph_lin(7) = 0.000071
           aph_lin(8) = 0.0
 
c          Case 1
c           aph_lin(1) = 0.01   !0.00665  
c           aph_lin(2) = 0.05582    
c           aph_lin(3) = 0.02055   
c           aph_lin(4) = 0.01910    
c           aph_lin(5) = 0.01015    
c           aph_lin(6) = 0.01424
c           aph_lin(7) = 0.0  
c           aph_lin(8) = 0.0

      ENDIF

      d_adm_s = DBLE(adm_s)
      d_bbp_s = DBLE(bbp_s)
      DO jl=1,NLAMBDA
          d_aw(jl)      = DBLE(aw(jl))
          d_bbw(jl)     = DBLE(bbw(jl))
          d_aph_lin(jl) = DBLE(aph_lin(jl))
      ENDDO

      RETURN
      END

C ############################################################################
C ############################################################################

C SECTION 8:

c radiometric correction routines.

      SUBROUTINE INBAND (rho_band7,rho_band8,fi6_out,fi7_out,fi8_out)

c derivation of 'f' correction factor for SeaWiFS bands 6 to 8.  
c Gordon (1995) conversion of La band centered to band averaged values
c (or vice versa).

      INTEGER    NLAMBDA
      PARAMETER (NLAMBDA=8)

      REAL   r765,r865

      REAL   theta0,theta,delt_phi
      REAL   M,mu,mu0,pi
      REAL   water_vapor_c,w_vap               !g/cm^2

      REAL   rho_band7,rho_band8
      REAL   fi6_out,fi7_out,fi8_out

      REAL   epsilon_78(4),c78(4)

      REAL   fi_band6(3),fi_band7(3),fi_band8(3)
      REAL   fi6(12),fi7(12),fi8(12)

      DATA   fi6 /    0.9986, -0.0007046,     2.459,  0.002545,
     &             -0.001644,  -0.001188,  -0.01015, -0.008021,
     &             0.0001378,  0.0001079,  0.001233, 0.0004105   /

      DATA   fi7 /    0.9983, -0.0008214,   -0.4094,   0.03732,
     &             -0.003537,  -0.001303,    0.1767,  0.008578,
     &             0.0003686,  0.0001534,  -0.02471, -0.002145   /

      DATA   fi8 /    0.9958,  -0.001561,     6.442,  -0.01894,
     &             -0.006337,  -0.002679,  -0.01037,  -0.03583, 
     &             0.0006157,  0.0003080, -0.002428,  0.003628   /


      COMMON /SPEC_GEOM/   theta0,theta,delt_phi
      COMMON /watv/        water_vapor_c

      pi       = 4.0*atan(1.0)
      mu = COS(theta * pi / 180.0)
      mu0 = COS(theta0 * pi / 180.0)
      M  = 1./mu0 + 1./mu

      w_vap=water_vapor_c

      r765 = 765.
      r865 = 865.

c     make sure epsilon does not go negative
      IF((rho_band7 .LE. 0.0) .OR. (rho_band8 .LE. 0.0)) then
         epsilon_78(1) = 1.0
      else
         epsilon_78(1) = rho_band7/rho_band8 !band centred
      endif

      c78(1) =  (ALOG(epsilon_78(1)))/(r865-r765)               !band centred 765 to 865

      do i = 2,4 

        fi_band6(i-1) = 
     &         (fi6(1)+fi6(2)*M)+(fi6(3)+fi6(4)*M)*c78(i-1)
     &      +( (fi6(5)+fi6(6)*M)+(fi6(7)+fi6(8)*M)*c78(i-1) )* w_vap
     &   +( (fi6(9)+fi6(10)*M)+(fi6(11)+fi6(12)*M)*c78(i-1) )*(w_vap**2)

        fi_band7(i-1) =  
     &         (fi7(1)+fi7(2)*M)+(fi7(3)+fi7(4)*M)*c78(i-1)
     &      +( (fi7(5)+fi7(6)*M)+(fi7(7)+fi7(8)*M)*c78(i-1) )* w_vap
     &   +( (fi7(9)+fi7(10)*M)+(fi7(11)+fi7(12)*M)*c78(i-1) )*(w_vap**2)

        fi_band8(i-1) =  
     &         (fi8(1)+fi8(2)*M)+(fi8(3)+fi8(4)*M)*c78(i-1)
     &      +( (fi8(5)+fi8(6)*M)+(fi8(7)+fi8(8)*M)*c78(i-1) )* w_vap
     &   +( (fi8(9)+fi8(10)*M)+(fi8(11)+fi8(12)*M)*c78(i-1) )*(w_vap**2)

        epsilon_78(i) = (fi_band7(i-1)*epsilon_78(i-1))/fi_band8(i-1)     !band averaged
        
c     make sure epsilon does not go negative
        IF( epsilon_78(i) .LE. 0.0) THEN
           epsilon_78(i) = 1.0
        ENDIF
        
        c78(i) =  (ALOG(epsilon_78(i)))/(r865-r765) !band averaged 765 to 865

      enddo

c      OPEN(unit=76,file='fi.txt',status='UNKNOWN',
c     &   form='FORMATTED',access='APPEND')
c         write(76,*) w_vap,epsilon_78(3),c78(3),
c     &                              fi_band6(3),fi_band7(3),fi_band8(3)
c      CLOSE (76)

      fi6_out = fi_band6(3)
      fi7_out = fi_band7(3)
      fi8_out = fi_band8(3)

      RETURN
      END


C ############################################################################
C ############################################################################


C SECTION 9: 

C optimisation routines for NIR

C **********************************************************************
C $RCSfile: zxmwd.f,v $
C
C                      P R O J E C T:  Z X M W D
C
C ----------------------------------------------------------------------
C   PURPOSE:            - Extension to IMSL routine ZXMWD
C ----------------------------------------------------------------------
CVS:
C 
C $RCSfile: zxmwd.f,v $
C $Revision: 2.18 $
C $Date: 2000/04/14 23:25:51 $
C $Author: chomko $
C
C $Header: /usr/users/chomko/Project/ZXMWD/zxmwd.f,v 2.18 2000/04/14 23:25:51 chomko Exp $
C
CVS-END
C ----------------------------------------------------------------------
C
C   COMPUTER            - DEC Alpha OSF1 / SINGLE
C
C   VERSION             - 2
C   LATEST RELEASE      - 2.1           (tag: rel-2_1)
C
C   DATE                - April 14, 2000
C
C   PURPOSE             - GLOBAL MINIMUM (w/ BOUND-REGION CONSTRAINTS) OF
C                           FUNCTION OF N VARIABLES, w/ CONTROLLED HESSIAN
C                           SETUP, "TRUE" NUMBER OF SEARCH-POINTS AND
C                           ADDED RESTARTS.
C
C   USAGE               - CALL ZXMWD(FCN,IOPT,N,NSIG,NRST,A,B,NSRCH,X,F,
C                                    WORK,IWORK,IER)
C
C   ARGUMENTS    FCN    - A USER SUPPLIED <LOGICAL FUNCTION> WHICH CALCULATES
C                           F GIVEN X(1),X(2),...,X(N).
C                           FCN IS REFERENCED AS FOLLOWS,
C                             lResult = FCN(N,X,F) (LOGICAL lResult)
C                             WHERE X IS A VECTOR OF LENGTH N
C                             FCN MUST APPEAR IN AN EXTERNAL STATEMENT
C                               IN THE CALLING PROGRAM.
C                             FCN MUST NOT ALTER THE VALUES OF
C                               X(I),I=1,...,N,  OR N.
C
C                         (LOGICAL RESULT):
C                                          .TRUE.  - success
C                                          .FALSE. - cannot be computed
C
C                IOPT   - INTEGER VARIABLE which indicates how to set-up
C                         Hessian Matrix before starting iterations. To use
C                         the original ZXMWD set IOPT = 0,1,2,3 (any)
C
C                                NOTE: Restarts ARE NOT available
C                                      for IOPT=0,...,3
C                         IOPT=0
C                         IOPT=1
C                         IOPT=2
C                         IOPT=3   Integer VARIABLE (INTERNAL USE ONLY) 
C                                  NOTE: some (any) of them MUST be set if
C                                        the original ZXMWD is to be used. 
C
C                         IOPT=4   Same in spirit as IOPT=0,...,3 BUT:
C
C                                  IF NSRCH <  5, final stage of optimization
C                                                 is carried out w/ NSRCH #of
C                                                 starting points;
C                                  IF NSRCH >= 5, final stage of optimization
C                                                 is carried out w/ NSRCH = 5
C                                                 #of starting points;
C                                  IF NSRCH = any, BEST HESSIAN is saved in
C                                                             WORK(NHH)-array
C                                                  AND RESTARTS are performed
C                                                  if IER=129,130,135 (i.e.,
C                                                  saddle,non-convergence,stuck.)
C
C                         IOPT=5   Optimization is carried out w/ a start.pt.
C                                  given in X vector AND INITIAL HESSIAN is
C                                  set to UNITY. Final (best) HESSIAN is saved.
C                                  NOTE:  if IER=129,130,135 as a result, restarts
C                                         are performed w/ the #of NSRCH start.pt's
C                                         This is to say,
C                                               NSRCH must not necessarily be 1
C
C                         IOPT=6   same as IOPT=5, but the INITIAL HESSIAN is
C                                                             actually computed.
C
C                         IOPT=7   (REQUIRES some care from the user, but pays-off)
C
C                                  Optimization is carried out w/ a start.pt.
C                                  given in X vector AND INITIAL HESSIAN taken
C                                  from the LATEST BEST HESSIAN (saved during a
C                                  previous call to ZXMWD).
C                                  NOTE: if there were no previous calls EVER to
C                                        ZXMWD, then it will INTERNALLY proceed
C                                        w/ IOPT=4, NSRCH=NSRCH.
C                                        Otherwise, if, in the process, you screw up
C                                                   HESSIAN saved in WORK array,
C                                                   YOYO (you're on your own) and
C                                                   expect miracles (e.g., core dump)
C
C                N      - THE NUMBER OF UNKNOWN PARAMETERS. (INPUT)
C
C                NSIG   - CONVERGENCE CRITERION. (INPUT)  NSIG IS THE
C                           NUMBER OF DIGITS OF ACCURACY REQUIRED IN
C                           THE PARAMETER ESTIMATES.
C                           NOTE: the actual accuracy of computation is provided
C                                 in COMMON /ZXMWD_STIX/ Block in the REAL variable
C                                 RSIG (might be used to reduce NSIG for acceleration),
C                                 BUT, reducing NSIG might lead to IER=129,130, thus
C                                 increasing the number of IFCN_TOT (see below).
C
C                NRST   - INTEGER Number of RESTARTS (to use w/ IOPT=4,...,7)
C                         NOTE: if IOPT=0,...,3, set NRST = whatever_integer >= 0
C
C                A,B    - CONSTRAINT VECTORS OF LENGTH N. (INPUT)
C                           X(I) IS REQUIRED TO SATISFY -
C                                A(I) .LE. X(I) .LE. B(I)
C
C                NSRCH  - NUMBER OF STARTING POINTS TO BE GENERATED.
C                           (INPUT) SUGGESTED VALUE = MIN(2**N+5,100)
C
C                X      - VECTOR OF LENGTH N CONTAINING THE FINAL
C                           PARAMETER ESTIMATES. (OUTPUT)
C
C                F      - VALUE OF THE FUNCTION AT THE FINAL
C                           PARAMETER ESTIMATES. (OUTPUT)
C
C                WORK   - REAL WORK VECTOR OF LENGTH (N+12)*N
C
C                IWORK  - INTEGER WORK VECTOR OF LENGTH N
C
C                IER    - ERROR PARAMETER (OUTPUT)
C                         TERMINAL ERROR:
C
C                           IER = 129 INDICATES THAT THE ALGORITHM
C                             HAS CONVERGED TO A POINT WHICH MAY
C                             ONLY BE A SADDLE POINT.
C
C                           IER = 130 INDICATES THAT IT WAS NOT
C                             POSSIBLE TO CALCULATE THE SOLUTION
C                             TO NSIG DIGITS.  SEE REMARKS.
C
C                           IER = 131 INDICATES THAT THE ITERATION
C                             WAS TERMINATED AFTER 200*(N+1)
C                             FUNCTION EVALUATIONS.  SEE REMARKS.
C
C                           IER = 132 INDICATES THAT A(I).GE.B(I)
C                             FOR SOME I=1,...,N. NO ATTEMPT IS MADE
C                             TO FIND THE MINIMUM IN THIS CASE.
C
C                           IER = 133 INDICATES THAT either
C                             X(I) < A(I) .OR. X(I) > B(I)
C                             .OR. B(I)-A(I) <= REPS (rmach-zero)
C
C                           IER = 134 INDICATES THAT NSIG.LT.1
C
C                           IER = 135 INDICATES THAT while iterating
C                             some point: ITN >= 2*IHUPD (i.e., the
C                             number of iterations exceeded the double
C                             of actual Hessian updates meaning some
C                             problems with gradients and step lengths,
C                             e.g.,"STUCK-CRITERION").
C
C                           IER = 136 INDICATES THAT LOGICAL FUNCTION FCN
C                             could not be computed for a given and/or
C                             computed set of points.
C
C
C   COMMON BLOCKS        - /ZXMWD_STIX/
C
C =============================================== COMMON BLOCK: ZXMWD-STIX
C         DESCRIPTION: to be included in USER routines (if needed)
C
C                   IFCN_TOT   -  actual NUMBER OF FUNCTION EVALUATIONS
C                   ITER_TOT   -  actual NUMBER OF ITERATIONS
C                   IHUPD_TOT  -  actual NUMBER OF HESSIAN UPDATES
C                   IRST_TOT   -  actual NUMBER OF RESTARTS
C                   RSIG       -  actual PRECISION (same as NSIG but REAL)
C ------------------------------------------------------------------------
C      REAL                RSIG
C      INTEGER             IFCN_TOT,ITER_TOT,IHUPD_TOT,IRST_TOT
C      COMMON /ZXMWD_STIX/ IFCN_TOT,ITER_TOT,IHUPD_TOT,IRST_TOT,RSIG
C ========================================================================
C
C   PRECISION/HARDWARE  - SINGLE/H64
C
C   REQD. ROUTINES      - ZXMWE, ZSRCH (included)
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      WHEN IER IS RETURNED AS 130 OR 131, THE PARAMETER
C                ESTIMATES IN X MAY NOT BE RELIABLE.  FURTHER CHECKING
C                SHOULD BE PERFORMED.  USE OF A LARGER NSRCH VALUE MAY
C                PRODUCE MORE RELIABLE PARAMETER ESTIMATES OR TRY
C                USING IOPT=4,...,7
C
C   COPYRIGHT           - 1984 BY IMSL, INC. ALL RIGHTS RESERVED.
C                         2000 R.M.Chomko U.Miami-PHYSICS, All Rights Reserved
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C                         W/ IOPT=0,...,3 GOTO IMSL WARRANTY, otherwise,
C                           there is no warranty for illegal redistribution.
C
C-----------------------------------------------------------------------------
CVS:
C $Id: zxmwd.f,v 2.18 2000/04/14 23:25:51 chomko Exp $
C
C                ========== HISTORY ==========
C $Log: zxmwd.f,v $
C Revision 2.18  2000/04/14 23:25:51  chomko
C    Update for NW5 variable (using PARAMETER NFINSRCH).
C
C Revision 2.17  2000/04/14 23:13:32  chomko
C    Extreme case fix: in RESTART loop, IER=136 iff F >= BIG. Otherwise,
C    it might have blown good X from IOPT=5,6,7 processing (if any).
C
C Revision 2.16  2000/04/14 21:50:44  chomko
C    Looks good. Preparing for Release-Version 2.1 (tag: rel-2_1).
C
C Revision 2.15  2000/04/14 21:45:11  chomko
C    Erased the original Header for ZXMWD from IMSL (i.e., v1.0)
C
C Revision 2.14  2000/04/14 21:41:04  chomko
C    Cleaned out of unnecessary oldNSRCH variable.
C
C Revision 2.13  2000/04/14 21:12:21  chomko
C    Bug fix: variable F always contains the value for the best X, while
C    array FAST(1:J2) contains best FQ's for A GIVEN RESTART, i.e.,
C    immaterial of previous restarts (if any). F changes only when
C    FQ < F. Those things were not quite followed well.
C
C Revision 2.12  2000/04/14 20:50:52  chomko
C    Extension to Revision 2.7:
C    NSRCHMAX is now properly handled for different IOPT's.
C    A reminder: ZSRCH routine bases its search on NSRCHMAX that
C    was incorrect previously (though, functional most of
C    the time). ZSRCH covers the hyper-rectangle completely in a
C    pre-determined manner (e.g., not randomly).
C
C Revision 2.11  2000/04/14 20:17:59  chomko
C    The number of SRCH points to look for in ZSRCH is now properly
C    handled for IOPT<4 (original IMSL ZXMWD).
C
C Revision 2.10  2000/04/14 18:48:02  chomko
C    In the loop IOPT<=4, if J2==0, set IER=136 and GO TO RESTART (if any).
C
C Revision 2.9  2000/04/14 18:19:00  chomko
C    When IOPT>4, the results were not prepared in FAST(J), thus, damaging
C    possibly the best result (sometimes).
C
C Revision 2.8  2000/04/14 17:40:01  chomko
C    Hard-wired number 5 is now a PARAMETER (NFINSRCH=5). Note: not variable.
C
C Revision 2.7  2000/04/14 17:27:25  chomko
C    1.Search routine ZSRCH now properly handles the MAX number of
C      starting points it should serve (NSRCH5 >> NSRCHMAX=NSRCH*NRST);
C    2.IER=0 has been moved out of the search loop for IOPT<=4, thus,
C    reserving IER for whatever best answer found.
C
C Revision 2.6  2000/04/14 00:25:54  chomko
C    Correction to Revision 2.5: ZXMWD actually should be ZXMWE, sorry.
C
C Revision 2.5  2000/04/14 00:21:45  chomko
C    Inconsistency Fix and Speedup (related to processing IER=135):
C    Before:
C      w/ IOPT<=4, for all NSRCH, if after ZXMWE(IOPT_in=0) FQ < F, then
C      it would do another but ZXMWE(IOPT_in=3), and keep results as final.
C    Now:
C      it will do ZXMWE(IOPT_in=3) iff JER[ZXMWE(IOPT_in=0)] != 0.
C    Facts FOR:
C      1). After IOPT_in=0, solution XX is almost the correct one or
C          at least a very close one. Near this point FCN most probably
C          will get stuck, thus giving rise to IER=135 (unwanted);
C      2). Well, that was slowing the process down.
C    Facts AGAINST: N/A.
C
C Revision 2.4  2000/04/13 23:08:56  chomko
C    1. J2 (final #of optimization points) now respects JER136 out of NSRCH5;
C    2. Restarts are made available for IER=136 (w/ IOPT >= 4).
C
C Revision 2.3  2000/04/13 22:11:01  chomko
C    1. IER=132 (A(j)>B(j)) is controlled just once per call;
C    2. Modified help for FCN and IER=136.
C
C Revision 2.2  2000/04/13 21:29:22  chomko
C    Assuming that it's not always possible to compute EXTERNAL FCN
C    function. So, that the following changes are made:
C    1. SUBROUTINE FCN  ->>> LOGICAL FUNCTION FCN (.TRUE. for SUCCESS);
C    2. If FCN==.FALSE., then IER=136.
C
C Revision 2.1  2000/04/12 22:41:58  chomko
C    Bringing ALL up to Revision 2.1
C
C Revision 2.0.0.1  2000/04/12 22:39:57  chomko
C    Making CVS support, v2.0.0
C
C $Header: /usr/users/chomko/Project/ZXMWD/zxmwd.f,v 2.18 2000/04/14 23:25:51 chomko Exp $
C
CVS-END
C
CSRC: **********************************************************************

      SUBROUTINE ZXMWD (FCN,IOPT,N,NSIG,NRST,A,B,NSRCH,X,F,
     &                  WORK,IWORK,IER)

c input:  fcn,iopt,n,nsig,nrst,a,b,nsrch
c output: x,f,work,iwork,ier

C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NFINSRCH
      PARAMETER         (NFINSRCH=5)
      INTEGER            N,NSIG,NSRCH,IER,IWORK(1)
      REAL               A(N),B(N),X(N),F,WORK(1)

C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I1,I2,I5,IOPT,IOPT_in,IP,I,
     *                   J1,J2,JER136,JER,J,
     *                   KK,K,NG,NH,NSRCH5,NS,NW,NX,IW(9)
      REAL               BIG,FAST(NFINSRCH),FQ,PI,TEMP
      REAL               REPS
      DATA               PI /3.141593/
      DATA               BIG /1.7E+38/, REPS /1.7E-37/
      SAVE               PI,BIG,REPS

      LOGICAL            FCN     
      EXTERNAL           FCN 

      LOGICAL       goodHess, wasErr
      DATA          goodHess /.FALSE./
      SAVE          goodHess
      INTEGER       RESTART,NRST,NHsz,NSRCHMAX
      INTEGER       oldIOPT
      REAL          RANGE,SIGEPS,XJ
C ------------------------------------------------------------- ZXMWD-STATISTICS
      REAL                RSIG
      INTEGER             IFCN_TOT,ITER_TOT,IHUPD_TOT,IRST_TOT
      COMMON /ZXMWD_STIX/ IFCN_TOT,ITER_TOT,IHUPD_TOT,IRST_TOT,RSIG
C ------------------------------------------------------------------------------

C                                  FIRST EXECUTABLE STATEMENT
      IF (NSRCH.LT.1) THEN
        PRINT*, "WARNING-ZXMWD: provide  NSRCH > 0"
        PRINT*, "               Assumed: NSRCH=5"
        NSRCH = NFINSRCH
      ENDIF
      IF (NRST.LT.0) THEN
        PRINT*, "WARNING-ZXMWD: provide  NRST >= 0"
        PRINT*, "               Assumed: NRST=2"
        NRST = 2
      ENDIF
      IF(IOPT.GT.7.OR.IOPT.LT.0) THEN
        PRINT*, "WARNING-ZXMWD: provide  0 <= IOPT <= 7"
        PRINT*, "               Assumed: IOPT=0"
        IOPT=0
      ENDIF

      IFCN_TOT  = 0         ! #of FCN function evaluations
      ITER_TOT  = 0         ! #of Iterations
      IHUPD_TOT = 0         ! #of Hessian Updates
      IRST_TOT  = 0         ! #of Restarts

      RSIG     = 0.0
      oldIOPT  = IOPT
      IP       = 0
      IER      = 0

      DO J=1,N
        IF (A(J).GE.B(J)) THEN
          PRINT*, "ERROR-ZXMWD: provide region A(J) < B(J)"
          IER=132
          GO TO 911
        ENDIF
      ENDDO

      RESTART  = -1          ! on the first run it's not a restart

      SIGEPS   = 10.0**(-NSIG)
      IF(SIGEPS.GE.1) THEN
        IER=134
        GO TO 911
      ENDIF

      IF(IOPT.EQ.7.AND..NOT.goodHess) IOPT=4

      IF (IOPT.LT.4) THEN
        NSRCHMAX = MAX0(NSRCH,NFINSRCH)    ! for ZSRCH routine
      ELSE IF (IOPT.EQ.4.) THEN
        NSRCHMAX = NSRCH*(NRST+1)
      ELSE
        NSRCHMAX = NSRCH*NRST
      ENDIF
C                          <<<<<<< RESTART >>>>>>>

    1 RESTART = RESTART+1

      IF (IOPT.GE.4) THEN
        J2 = MIN0(NSRCH,NFINSRCH)
      ELSE
        J2 = NFINSRCH
      ENDIF

      DO I=1,J2
         FAST(I) = BIG
      ENDDO

      NH   = N+1
      NG   = NH+N*(N+1)/2
      NX   = NG+N
      NS   = NX+N
      NW5  = NS+NFINSRCH*N
      NW   = NS+J2*N
      NHH  = NW5+3*N
      NHsz = N*(N+1)/2

      DO 10 I=NS,NW
         WORK(I) = 0.0
   10 CONTINUE

      IF (IOPT.LE.4) THEN

       NSRCH5 = MAX0(NSRCH,J2)
       JER136 = 0               ! the #of non-computable FCN's

       DO 35 I=1,NSRCH5
C                                  GENERATE STARTING POINTS
         CALL ZSRCH(A,B,N,NSRCHMAX,IP,WORK,IWORK,IW,JER)

         DO 15 J=1,N
            XJ = WORK(J)
            RANGE = SIGEPS*(B(J)-A(J))
            IF ((B(J)-A(J)).GT.REPS) THEN
              IF ((XJ-A(J)).LT.RANGE) WORK(J)=A(J)+RANGE  ! move off LeftBoundary
              IF ((B(J)-XJ).LT.RANGE) WORK(J)=B(J)-RANGE  ! move off RightBoundary
            ELSE
              IER=133
              GO TO 911
            ENDIF
            WORK(J) = PI/2.*(WORK(J)-A(J))/(B(J)-A(J))
   15    CONTINUE

C                                  DO 4 ITERATIONS WITH EACH
crc         CALL ZXMWE(FCN,N,NSIG,7*(N+1)+5*N+1,2,WORK,WORK(NH),WORK(NG),

         CALL ZXMWE(FCN,N,NSIG,4*(N+1)+2*N+1,2,WORK,WORK(NH),WORK(NG),
     *              FQ,WORK(NW),A,B,WORK(NX),WORK(NHH),JER)

         IF(JER.EQ.136) THEN
           JER136 = JER136 + 1
           GO TO 35
         ENDIF
         IF (FQ.GE.FAST(J2)) GO TO 35
         FAST(J2) = FQ
         DO 20 K=1,N
            I5 = NS+(J2-1)*N+K-1
            WORK(I5) = AMOD(WORK(K),PI)
   20    CONTINUE
         DO 30 J1=1,J2-1
            IF (FAST(J2).GE.FAST(J1)) GO TO 30
            TEMP = FAST(J1)
            FAST(J1) = FAST(J2)
            FAST(J2) = TEMP
            DO 25 K=1,N
               I1 = NS+(J1-1)*N+K-1
               I2 = NS+(J2-1)*N+K-1
               TEMP = WORK(I1)
               WORK(I1) = WORK(I2)
               WORK(I2) = TEMP
   25       CONTINUE
   30    CONTINUE
   35  CONTINUE

       IF(RESTART.EQ.0) F = BIG      ! keep the best from previous restarts
       J2  = MIN0(NSRCH5-JER136,J2)

       IF (J2.EQ.0) THEN  ! all srch points led to FCN=.FALSE., hence, restart.
         IF (F.GE.BIG) IER = 136
         GO TO 9010
       ENDIF

       DO 60 I=1,J2
         DO 40 K=1,N
            KK = NS+(I-1)*N+K-1
            WORK(K) = WORK(KK)
   40    CONTINUE
C                                  DO 200 ITERATIONS WITH THE 5
C                                    STARTING POINTS WHICH GAVE
C                                    THE SMALLEST SUM OF SQUARES
C                                    AFTER 4 ITERATIONS
         IOPT_in = 0
         CALL ZXMWE(FCN,N,NSIG,200*(N+1),IOPT_in,WORK,WORK(NH),
     *              WORK(NG),FQ,WORK(NW),A,B,WORK(NX),
     &              WORK(NHH),JER)
         IF (FQ.GE.F.OR.JER.EQ.136) GO TO 60
         IF (JER.NE.0) THEN
           IOPT_in = 3
           CALL ZXMWE(FCN,N,NSIG,200*(N+1),IOPT_in,WORK,WORK(NH),
     *                WORK(NG),FQ,WORK(NW),A,B,WORK(NX),
     &                WORK(NHH),JER)
         ENDIF
         IER = JER
         F = FQ
         DO 55 J=1,N
            X(J) = A(J)+(B(J)-A(J))*SIN(WORK(J))**2
   55    CONTINUE
         IF (oldIOPT.GE.4) THEN
           DO 56 J=1,NHsz
              K  = NH+J-1
              KK = NHH+J-1
              WORK(KK) = WORK(K)     ! save the BEST HESSIAN
   56      CONTINUE
         ENDIF
C                                  SAVE NSIG ESTIMATE FOR BEST X
         RSIG = WORK(NW+2)
   60  CONTINUE
C                                  RETURN NSIG ESTIMATE IN WORK(1)
 9000   CONTINUE
        IF(oldIOPT.LT.4) WORK(1) = RSIG

      ELSE  ! IOPT=5,6,7

        IER = 132
        DO J=1,N
          IF (A(J).GE.B(J)) GO TO 911
          XJ = X(J)
          IF (XJ.LT.A(J).OR.XJ.GT.B(J)) THEN
             IER=133
             GO TO 911
          ENDIF
          RANGE = SIGEPS*(B(J)-A(J))
          IF ((B(J)-A(J)).GT.REPS) THEN
            IF ((XJ-A(J)).LT.RANGE) X(J)=A(J)+RANGE  ! move off LeftBoundary
            IF ((B(J)-XJ).LT.RANGE) X(J)=B(J)-RANGE  ! move off RightBoundary
          ELSE
            IER=133
            GO TO 911
          ENDIF
        ENDDO

        DO 70 J=1,N
          WORK(J) = ASIN(SQRT((X(J)-A(J))/(B(J)-A(J))))
   70   CONTINUE

        IER=0
        IF     (IOPT.EQ.7) THEN
            IOPT_in = 7
        ELSE IF(IOPT.EQ.6) THEN      ! Honestly compute Hessian
            IOPT_in = 3
        ELSE IF(IOPT.EQ.5) THEN      ! Set Hessian to Unity
            IOPT_in = 0
        ENDIF
        CALL ZXMWE(FCN,N,NSIG,200*(N+1),IOPT_in,WORK,WORK(NH),
     *               WORK(NG),FQ,WORK(NW),A,B,WORK(NX),
     &               WORK(NHH),IER)
        F = FQ
        DO 75 J=1,N
          X(J) = A(J)+(B(J)-A(J))*SIN(WORK(J))**2
   75   CONTINUE
        DO J=1,NHsz
          K  = NH+J-1
          KK = NHH+J-1
          WORK(KK) = WORK(K)     ! save HESSIAN
        ENDDO
        RSIG = WORK(NW+2)        ! actual precision (REAL)

      ENDIF ! IOPT

 9010 CONTINUE
      IF (oldIOPT.GE.4) THEN
        wasErr=.FALSE.
        IF (IER.EQ.129.OR.
     &      IER.EQ.130.OR.
     &      IER.EQ.135.OR.
     &      IER.EQ.136    ) wasErr=.TRUE.
        IF (wasErr.AND.RESTART.LT.NRST) THEN
          IOPT  = 4
          GOTO 1                  ! RESTART
        ELSE                      ! FINISH
          IF (wasErr) THEN
            goodHess = .FALSE.  ! bad Hessian
          ELSE
            goodHess = .TRUE.   ! tolerable Hessian
          ENDIF
          IRST_TOT = IRST_TOT + RESTART
        ENDIF
      ENDIF

  911 CONTINUE
      IOPT  = oldIOPT

      RETURN
      END

CSRC-END
C ********************************************************************
 
C **********************************************************************
C $RCSfile: zxmwe.f,v $
C
C                      P R O J E C T:  Z X M W D
C
C ----------------------------------------------------------------------
C   PURPOSE             - Extension to IMSL routine ZXMWD
C ----------------------------------------------------------------------
CVS:
C 
C $RCSfile: zxmwe.f,v $
C $Revision: 2.2 $
C $Date: 2000/04/13 21:31:12 $
C $Author: chomko $
C
C $Header: /usr/users/chomko/Project/ZXMWD/zxmwe.f,v 2.2 2000/04/13 21:31:12 chomko Exp $
C ----------------------------------------------------------------------
C
C   IMSL ROUTINE NAME   - ZXMWE
C
C   COMPUTER            - VAX/SINGLE
C
C   LAST IMSL REVISION  - JUNE 1, 1982
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY ZXMWD
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO,ZXMJN
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1982 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
CVS:
C $Id: zxmwe.f,v 2.2 2000/04/13 21:31:12 chomko Exp $
C
C                ========== HISTORY ==========
C $Log: zxmwe.f,v $
C Revision 2.2  2000/04/13 21:31:12  chomko
C Changed: SUBROUTINE FCN ->>> LOGICAL FUNCTION FCN, such that,
C if FCN==.FALSE., then function cannot be computed for given or
C computed point and IER=136 as a result.
C
C Revision 2.1  2000/04/12 22:41:58  chomko
C Bringing ALL up to Revision 2.1
C
C Revision 2.0.0.1  2000/04/12 22:39:57  chomko
C Making CVS support, v2.0.0
C
C $Header: /usr/users/chomko/Project/ZXMWD/zxmwe.f,v 2.2 2000/04/13 21:31:12 chomko Exp $
C **********************************************************************

      SUBROUTINE ZXMWE (FUNCT,N,NSIG,MAXFN,IOPT,
     &                  X,H,G,F,W,A,B,XX,HOLD,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,NSIG,MAXFN,IOPT,IER
      REAL               X(N),H(1),G(N),F,W(3),A(1),B(1),XX(1),HOLD(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IDIFF,IFN,IGG,IG,II,IJ,IM1,IR,IS,ITN,I,JB,JJ,
     *                   JNT,JP1,J,KJ,K,LINK,L,NJ,NM1,NP1,
     &                   NHsz,IHUPD
      REAL               AEPS,ALPHA,AX,DF,DGS,DIFF,EPS,F11,F12,F1,F21,
     *                   F22,F2,FF,FIVE,GHH,GNRM,GS0,GYS,H2,HALF,HHH,HH,
     *                   HJJ,HMAX,HMIN,ONE,P1,RELX,REPS,SEVEN,SIG,TEN,
     *                   TOT,TWELVE,V,ZERO,ZZ,Z
      DATA               REPS /1.1921E-07/,AX /0.1/
      DATA               ZERO /0.0/,ONE /1.0/,HALF /0.5/,SEVEN
     *                   /7.0/,FIVE /5.0/,TWELVE /12.0/,
     *                   TEN /10.0/,P1 /0.1/

      LOGICAL   FUNCT 
      EXTERNAL  FUNCT    
C ------------------------------------------------------------- ZXMWD-STATISTICS
      REAL                RSIG
      INTEGER             IFCN_TOT,ITER_TOT,IHUPD_TOT,IRST_TOT
      COMMON /ZXMWD_STIX/ IFCN_TOT,ITER_TOT,IHUPD_TOT,IRST_TOT,RSIG
C ------------------------------------------------------------------------------

      INTEGER      uf           ! history file unit
      LOGICAL      lHistory     ! do output for history?
      CHARACTER*20 fHistName    ! history file name

      REAL         XX0

      lHistory  = .FALSE.       
      fHistName = 'zxmwh.out'

      IF (lHistory) CALL ZXMWH_INI (uf,fHistName)
C
C                                  INITIALIZATION
C                                  FIRST EXECUTABLE STATEMENT
      IER   = 0
      HH    = SQRT(REPS)
      H2    = SQRT(HH)
      EPS   = TEN**(-NSIG)
      IG    = N
      IGG   = N+N
      IS    = IGG
      IDIFF = 1
      IR    = N
      IHUPD = 0
      W(1)  = -ONE
      W(2)  = ZERO
      W(3)  = ZERO
      NHsz  = N*(N+1)/2
C                                  EVALUATE FUNCTION AT STARTING POINT
      DO 5 I=1,N
         G(I) = X(I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
    5 CONTINUE
      IF(.NOT.FUNCT(N,XX,F)) GO TO 9010  !first call to FCN_TAUVV_DIFFTRAN (N,x,f)
      IFN = 1                            !N,x sent to this routine, f retrieved from this routine
      IF (IOPT.EQ.1) GO TO 50
C                                  SET OFF-DIAGONAL ELEMENTS OF H TO 0.0
      IF (N.EQ.1) GO TO 20
C                                  Set OLD (PREVIOUS) Hessian
      IF (IOPT.EQ.7) THEN
         DO IJ=1,NHsz
           H(IJ) = HOLD(IJ)
         ENDDO
         GOTO 61                  ! H=L*D*L' decomposition
      ENDIF

      IJ = 2
      DO 15 I=2,N
         DO 10 J=2,I
            H(IJ) = ZERO
            IJ = IJ+1
   10    CONTINUE
         IJ = IJ+1
   15 CONTINUE
   20 IF (IOPT.NE.0) GO TO 30
C                                  SET DIAGONAL ELEMENTS OF H TO ONE
      IJ = 0
      DO 25 I=1,N
         IJ = IJ+I
         H(IJ) = ONE
   25 CONTINUE
      GO TO 95
C                                  GET DIAGONAL ELEMENTS OF HESSIAN
   30 IM1 = 1
      NM1 = 1
      NP1 = N+1
      DO 35 I=2,NP1
         HHH = H2*AMAX1(ABS(X(IM1)),AX)
         G(IM1) = X(IM1)+HHH
         XX(IM1) = A(IM1)+(B(IM1)-A(IM1))*SIN(G(IM1))**2
         IF(.NOT.FUNCT(N,XX,F2)) GO TO 9010
         G(IM1) = X(IM1)-HHH
         XX(IM1) = A(IM1)+(B(IM1)-A(IM1))*SIN(G(IM1))**2
         IF(.NOT.FUNCT(N,XX,FF)) GO TO 9010
         H(NM1) = (FF-F+F2-F)/(HHH*HHH)
         G(IM1) = X(IM1)
         XX(IM1) = A(IM1)+(B(IM1)-A(IM1))*SIN(G(IM1))**2
         IM1 = I
         NM1 = I+NM1
   35 CONTINUE
      IFN = IFN+N+N
      IF (IOPT.NE.3 .OR. N.EQ.1) GO TO 50
C                                  GET THE REST OF THE HESSIAN
      JJ = 1
      II = 2
      DO 45 I=2,N
         GHH = H2*AMAX1(ABS(X(I)),AX)
         DO 40 J=1,JJ
            HHH = H2*AMAX1(ABS(X(J)),AX)
            G(I) = X(I)+GHH
            XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
            G(J) = X(J)+HHH
            XX(J) = A(J)+(B(J)-A(J))*SIN(G(J))**2
            IF(.NOT.FUNCT(N,XX,F22)) GO TO 9010
            G(I) = X(I)-GHH
            XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
            IF(.NOT.FUNCT(N,XX,F12)) GO TO 9010
            G(J) = X(J)-HHH
            XX(J) = A(J)+(B(J)-A(J))*SIN(G(J))**2
            IF(.NOT.FUNCT(N,XX,F11)) GO TO 9010
            G(I) = X(I)+GHH
            XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
            IF(.NOT.FUNCT(N,XX,F21)) GO TO 9010
            H(II) = (F22-F21-F12+F11)/(4.*HHH*GHH)
            G(J) = X(J)
            XX(J) = A(J)+(B(J)-A(J))*SIN(G(J))**2
            II = II+1
   40    CONTINUE
         G(I) = X(I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(G(I))**2
         JJ = JJ+1
         II = II+1
   45 CONTINUE
      IFN = IFN+((N*N-N)*2)
C                                  ADD MULTIPLE OF IDENTITY TO
C                                  MAKE DIAGONAL ELEMENTS POSITIVE
   50 HMIN = H(1)
      HMAX = H(1)
      NM1 = 1
      DO 55 I=1,N
         HMIN = AMIN1(HMIN,H(NM1))
         HMAX = AMAX1(HMAX,H(NM1))
         NM1 = NM1+I+1
   55 CONTINUE
      HMIN = AMAX1(0.01*(ABS(HMAX)+ABS(HMIN))-HMIN,0.0)
      NM1 = 1
      DO 60 I=1,N
         H(NM1) = H(NM1)+HMIN
         NM1 = NM1+I+1
   60 CONTINUE
C                                  FACTOR H TO L*D*L-TRANSPOSE
   61 IR = N
      IF (N.GT.1) GO TO 65
      IF (H(1).GT.ZERO) GO TO 95
      H(1) = ZERO
      IR = 0
      GO TO 90
   65 NM1 = N-1
      JJ = 0
      DO 85 J=1,N
         JP1 = J+1
         JJ = JJ+J
         HJJ = H(JJ)
         IF (HJJ.GT.ZERO) GO TO 70
         H(JJ) = ZERO
         IR = IR-1
         GO TO 85
   70    IF (J.EQ.N) GO TO 85
         IJ = JJ
         L = 0
         DO 80 I=JP1,N
            L = L+1
            IJ = IJ+I-1
            V = H(IJ)/HJJ
            KJ = IJ
            DO 75 K=I,N
               H(KJ+L) = H(KJ+L)-H(KJ)*V
               KJ = KJ+K
   75       CONTINUE
            H(IJ) = V
   80    CONTINUE
   85 CONTINUE
   90 IF (IR.EQ.N) GO TO 95
      IER = 129
      GO TO 9000

   95 ITN = 0
      DF = -ONE
C                                  EVALUATE GRADIENT W(IG+I),I=1,...,N
  100 LINK = 1
      GO TO 275

  105 CONTINUE
      IF (lHistory) CALL ZXMWH_WRITE (uf,N,ITN,F,XX,W(IG),IFN)

C                                  BEGIN ITERATION LOOP
      IF (IFN.GE.MAXFN) GO TO 240
      ITN = ITN+1
      DO 110 I=1,N
         W(I) = -W(IG+I)
  110 CONTINUE
C                                  DETERMINE SEARCH DIRECTION W
C                                    BY SOLVING H*W = -G WHERE
C                                    H = L*D*L-TRANSPOSE
      IF (IR.LT.N) GO TO 140
C                                  N .EQ. 1
      G(1) = W(1)
      IF (N.GT.1) GO TO 115
      W(1) = W(1)/H(1)
      GO TO 140
C                                  N .GT. 1
  115 II = 1
C                                  SOLVE L*W = -G
      DO 125 I=2,N
         IJ = II
         II = II+I
         V = W(I)
         IM1 = I-1
         DO 120 J=1,IM1
            IJ = IJ+1
            V = V-H(IJ)*W(J)
  120    CONTINUE
         G(I) = V
         W(I) = V
  125 CONTINUE
C                                  SOLVE (D*LT)*Z = W WHERE
C                                  LT = L-TRANSPOSE
      W(N) = W(N)/H(II)
      JJ = II
      NM1 = N-1
      DO 135 NJ=1,NM1
C                                  J = N-1,N-2,...,1
         J = N-NJ
         JP1 = J+1
         JJ = JJ-JP1
         V = W(J)/H(JJ)
         IJ = JJ
         DO 130 I=JP1,N
            IJ = IJ+I-1
            V = V-H(IJ)*W(I)
  130    CONTINUE
         W(J) = V
  135 CONTINUE
C                                  DETERMINE STEP LENGTH ALPHA
  140 RELX = ZERO
      GS0 = ZERO
      DO 145 I=1,N
         W(IS+I) = W(I)
         DIFF = ABS(W(I))/AMAX1(ABS(X(I)),AX)
         RELX = AMAX1(RELX,DIFF)
         GS0 = GS0+W(IG+I)*W(I)
  145 CONTINUE
      IF (RELX.EQ.ZERO) GO TO 245
      AEPS = EPS/RELX
      IER = 130
      IF (GS0.GE.ZERO) GO TO 245
      IF (DF.EQ.ZERO) GO TO 245
      IER = 0
      ALPHA = (-DF-DF)/GS0
      IF (ALPHA.LE.ZERO) ALPHA = ONE
      ALPHA = AMIN1(ALPHA,ONE)
      IF (IDIFF.EQ.2) ALPHA = AMAX1(P1,ALPHA)
      FF = F
      TOT = ZERO
      JNT = 0
C                                  SEARCH ALONG  X+ALPHA*W
  150 IF (IFN.GE.MAXFN) GO TO 240
      DO 155 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(W(I))**2
  155 CONTINUE
      IF(.NOT.FUNCT(N,XX,F1)) GO TO 9010
      IFN = IFN+1
      IF (F1.GE.F) GO TO 180
      F2 = F
      TOT = TOT+ALPHA
  160 IER = 0
      F = F1
      DO 165 I=1,N
         X(I) = W(I)
  165 CONTINUE
      IF (JNT-1) 170, 200, 205
  170 IF (IFN.GE.MAXFN) GO TO 240
      DO 175 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(W(I))**2
  175 CONTINUE
      IF(.NOT.FUNCT(N,XX,F1)) GO TO 9010
      IFN = IFN+1
      IF (F1.GE.F) GO TO 205
      IF (F1+F2.GE.F+F .AND. SEVEN*F1+FIVE*F2.GT.TWELVE*F) JNT = 2
      TOT = TOT+ALPHA
      ALPHA = ALPHA+ALPHA
      GO TO 160
  180 CONTINUE
      IF (F.EQ.FF .AND. IDIFF.EQ.2 .AND. RELX.GT.EPS) IER = 130
      IF (ALPHA.LT.AEPS) GO TO 245
      IF (IFN.GE.MAXFN) GO TO 240
      ALPHA = HALF*ALPHA
      DO 185 I=1,N
         W(I) = X(I)+ALPHA*W(IS+I)
         XX(I) = A(I)+(B(I)-A(I))*SIN(W(I))**2
  185 CONTINUE
      IF(.NOT.FUNCT(N,XX,F2)) GO TO 9010
      IFN = IFN+1
      IF (F2.GE.F) GO TO 195
      TOT = TOT+ALPHA
      IER = 0
      F = F2
      DO 190 I=1,N
         X(I) = W(I)
  190 CONTINUE
      GO TO 200
  195 Z = P1
      IF (F1+F.GT.F2+F2) Z = ONE+HALF*(F-F1)/(F+F1-F2-F2)
      Z = AMAX1(P1,Z)
      ALPHA = Z*ALPHA
      JNT = 1
      GO TO 150
  200 IF (TOT.LT.AEPS) GO TO 245
  205 ALPHA = TOT
C                                  SAVE OLD GRADIENT
      DO 210 I=1,N
         W(I) = W(IG+I)
  210 CONTINUE
C                                  EVALUATE GRADIENT W(IG+I), I=1,...,N
      LINK = 2
      GO TO 275
  215 IF (IFN.GE.MAXFN) GO TO 240
      GYS = ZERO
      DO 220 I=1,N
         GYS = GYS+W(IG+I)*W(IS+I)
         W(IGG+I) = W(I)
  220 CONTINUE
      DF = FF-F
      DGS = GYS-GS0
      IF (DGS.LE.ZERO) GO TO 105
      IF (DGS+ALPHA*GS0.GT.ZERO) GO TO 230
C                                  UPDATE HESSIAN H USING
C                                    COMPLEMENTARY DFP FORMULA
      SIG = ONE/GS0
      IR = -IR
      CALL ZXMJN(H,N,W,SIG,G,IR,0,ZERO)
      DO 225 I=1,N
         G(I) = W(IG+I)-W(IGG+I)
  225 CONTINUE
      SIG = ONE/(ALPHA*DGS)
      IR = -IR
      CALL ZXMJN(H,N,G,SIG,W,IR,0,ZERO)
      IHUPD = IHUPD+1
      GO TO 105
C                                  UPDATE HESSIAN USING
C                                    DFP FORMULA
  230 ZZ = ALPHA/(DGS-ALPHA*GS0)
      SIG = -ZZ
      CALL ZXMJN(H,N,W,SIG,G,IR,0,REPS)
      Z = DGS*ZZ-ONE
      DO 235 I=1,N
         G(I) = W(IG+I)+Z*W(IGG+I)
  235 CONTINUE
      SIG = ONE/(ZZ*DGS*DGS)
      CALL ZXMJN(H,N,G,SIG,W,IR,0,ZERO)
      IHUPD = IHUPD+1
      GO TO 105
  240 IER = 131
C                                  MAXFN FUNCTION EVALUATIONS
      GO TO 250
  245 IF (IDIFF.EQ.2) GO TO 250
C                                  CHANGE TO CENTRAL DIFFERENCES
      IDIFF = 2
      GO TO 100
  250 IF (RELX.GT.EPS .AND. IER.EQ.0) GO TO 100
C                                  MOVE GRADIENT TO G AND RETURN
      GNRM = ZERO
      DO 255 I=1,N
         G(I) = W(IG+I)
         GNRM = GNRM+G(I)*G(I)
  255 CONTINUE
      GNRM = SQRT(GNRM)
      W(1) = GNRM
      W(2) = IFN
      W(3) = -ALOG10(AMAX1(REPS,RELX))
C                                  COMPUTE H = L*D*L-TRANSPOSE
      IF (N.EQ.1) GO TO 9000
      NP1 = N+1
      NM1 = N-1
      JJ = (N*(NP1))/2
      DO 270 JB=1,NM1
         JP1 = NP1-JB
         JJ = JJ-JP1
         HJJ = H(JJ)
         IJ = JJ
         L = 0
         DO 265 I=JP1,N
            L = L+1
            IJ = IJ+I-1
            V = H(IJ)*HJJ
            KJ = IJ
            DO 260 K=I,N
               H(KJ+L) = H(KJ+L)+H(KJ)*V
               KJ = KJ+K
  260       CONTINUE
            H(IJ) = V
  265    CONTINUE
         HJJ = H(JJ)
  270 CONTINUE
      GO TO 9000
C                                  EVALUATE GRADIENT
  275 IF (IDIFF.EQ.2) GO TO 290
C                                  FORWARD DIFFERENCES
C                                    GRADIENT = W(IG+I), I=1,...,N
      DO 280 I=1,N
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  280 CONTINUE
      DO 285 I=1,N
         Z = HH*AMAX1(ABS(X(I)),AX)
c         Z = 0.000001
         ZZ = X(I)
         X(I) = ZZ+Z
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
         IF(.NOT.FUNCT(N,XX,F1)) GO TO 9010
         W(IG+I) = (F1-F)/Z
         X(I) = ZZ
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  285 CONTINUE
      IFN = IFN+N
      GO TO (105, 215), LINK
C                                  CENTRAL DIFFERENCES
C                                    GRADIENT = W(IG+I), I=1,...,N
  290 DO 295 I=1,N
         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  295 CONTINUE
      DO 300 I=1,N
c         Z = HH*AMAX1(ABS(X(I)),AX)

         Z = 0.00001
         XX0 = A(I)+(B(I)-A(I))*SIN(X(I))**2
         IF ( (XX0-A(I)).LT.Z ) THEN
             XX0 = A(I)+Z
         ELSE IF ( (B(I)-XX0).LT.Z ) THEN
             XX0 = B(I)-Z
         ENDIF

         ZZ = XX0
         XX0 = ZZ+Z
         XX(I) = XX0
         IF(.NOT.FUNCT(N,XX,F1)) GO TO 9010

         XX0 = ZZ-Z
         XX(I) = XX0
         IF(.NOT.FUNCT(N,XX,F2)) GO TO 9010

         W(IG+I) = (B(I)-A(I))*SIN(2.0*X(I))*(F1-F2)/(Z+Z)
c         X(I) = ZZ
c         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
         XX(I) = ZZ

c         Z = HH*AMAX1(ABS(X(I)),AX)
cc         Z = 0.000001
c         ZZ = X(I)
c         X(I) = ZZ+Z
c         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
c         IF(.NOT.FUNCT(N,XX,F1)) GO TO 9010
c         X(I) = ZZ-Z
c         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
c         IF(.NOT.FUNCT(N,XX,F2)) GO TO 9010
c         W(IG+I) = (F1-F2)/(Z+Z)
c         X(I) = ZZ
c         XX(I) = A(I)+(B(I)-A(I))*SIN(X(I))**2
  300 CONTINUE
      IFN = IFN+N+N
      GO TO (105, 215), LINK
 9000 CONTINUE

      IF (IER.EQ.0.AND.ITN.GE.2*IHUPD) IER=135

      IFCN_TOT  =  IFCN_TOT  + IFN
      ITER_TOT  =  ITER_TOT  + ITN
      IHUPD_TOT =  IHUPD_TOT + IHUPD

      IF (lHistory) CALL ZXMWH_END (uf)
 9005 RETURN

 9010 IER = 136
      DO J=1,N
        X(J) = ASIN(SQRT((XX(J)-A(J))/(B(J)-A(J))))
      ENDDO
      IF (lHistory) CALL ZXMWH_END (uf)

      RETURN
      END

C ***********************************************************************
 
C **********************************************************************
C $RCSfile: zxmjn.f,v $
C
C                      P R O J E C T:  Z X M W D
C
C ----------------------------------------------------------------------
C   PURPOSE             - Extension to IMSL routine ZXMWD
C ----------------------------------------------------------------------
CVS:
C 
C $RCSfile: zxmjn.f,v $
C $Revision: 2.1 $
C $Date: 2000/04/12 22:41:58 $
C $Author: chomko $
C
C $Header: /usr/users/chomko/Project/ZXMWD/zxmjn.f,v 2.1 2000/04/12 22:41:58 chomko Exp $
C-----------------------------------------------------------------------
C
C   IMSL ROUTINE NAME   - ZXMJN
C
C   COMPUTER            - VAX/SINGLE
C
C   LAST IMSL REVISION  - NOVEMBER 1, 1984
C
C   PURPOSE             - NUCLEUS CALLED BY IMSL ROUTINES ZXMIN AND
C                           ZXMWD
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - DOUBLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - NONE REQUIRED
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
CVS:
C $Id: zxmjn.f,v 2.1 2000/04/12 22:41:58 chomko Exp $
C
C                ========== HISTORY ==========
C $Log: zxmjn.f,v $
C Revision 2.1  2000/04/12 22:41:58  chomko
C Bringing ALL up to Revision 2.1
C
C Revision 2.0.0.1  2000/04/12 22:39:57  chomko
C Making CVS support, v2.0.0
C
C $Header: /usr/users/chomko/Project/ZXMWD/zxmjn.f,v 2.1 2000/04/12 22:41:58 chomko Exp $
C **********************************************************************
C
      SUBROUTINE ZXMJN(A,N,Z,SIG,W,IR,MK,EPS)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,IR,MK
      REAL               A(1),Z(N),SIG,W(N),EPS
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            J,JJ,IJ,JP1,I,II,MM
      REAL               ZERO,ONE,FOUR,TI,V,TIM,AL,R,B,GM,Y,RINF,SQRINF
      DATA               ZERO/0.0/,ONE/1.0/,FOUR/4.0/
      DATA               RINF/1.7E+38/
C                                  UPDATE FACTORS GIVEN IN A
C                                    SIG*Z*Z-TRANSPOSE IS ADDED
C                                  FIRST EXECUTABLE STATEMENT
      SQRINF = SQRT(RINF)
      IF (N.GT.1) GO TO 5
C                                  N .EQ. 1
      A(1) = A(1) + SIG*Z(1)*Z(1)
      IR = 1
      IF (A(1).GT.ZERO) GO TO 9005
      A(1) = ZERO
      IR = 0
      GO TO 9005
C                                  N .GT. 1
    5 IF (SIG.GT.ZERO) GO TO 65
      IF (SIG.EQ.ZERO .OR. IR.EQ.0) GO TO 9005
      TI = ONE/SIG
      JJ = 0
      IF (MK.EQ.0) GO TO 15
C                                  L*W = Z ON INPUT
      DO 10 J = 1, N
        JJ = JJ + J
        IF (A(JJ).NE.ZERO) TI = TI + (W(J)*W(J))/A(JJ)
   10 CONTINUE
      GO TO 40
C                                  SOLVE L*W = Z
   15 DO 20 J = 1, N
        W(J) = Z(J)
   20 CONTINUE
      DO 35 J = 1, N
        JJ = JJ + J
        V = W(J)
        IF (A(JJ).GT.ZERO) GO TO 25
        W(J) = ZERO
        GO TO 35
   25   TI = TI + (V*V)/A(JJ)
        IF (J.EQ.N) GO TO 35
        IJ = JJ
        JP1 = J + 1
        DO 30 I = JP1, N
          IJ = IJ + I - 1
          W(I) = W(I) - V*A(IJ)
   30   CONTINUE
   35 CONTINUE
C                                  SET TI, TIM AND W
   40 IF (IR.LE.0) GO TO 45
      IF (TI.GT.ZERO) GO TO 50
      IF (MK-1) 65, 65, 55
   45 TI = ZERO
      IR = -IR - 1
      GO TO 55
   50 TI = EPS/SIG
      IF (EPS.EQ.ZERO) IR = IR - 1
   55 TIM = TI
      II = JJ
      I = N
      DO 60 J = 1, N
        IF (A(II).NE.ZERO) TIM = TI - (W(I)*W(I))/A(II)
        W(I) = TI
        TI = TIM
        II = II - I
        I = I - 1
   60 CONTINUE
      MM = 1
      GO TO 70
   65 MM = 0
      TIM = ONE/SIG
   70 JJ = 0
C                                  UPDATE A
      DO 120 J = 1, N
        JJ = JJ + J
        IJ = JJ
        JP1 = J + 1
C                                  UPDATE A(J,J)
        V = Z(J)
        IF (A(JJ).GT.ZERO) GO TO 95
C                                  A(J,J) .EQ. ZERO
        IF (IR.GT.0 .OR. SIG.LT.ZERO .OR. V.EQ.ZERO) GO TO 90
        IR = 1 - IR
        IF (ABS(V).GE.SQRINF) GO TO 75
        A(JJ) = (V*V)/TIM
        GO TO 80
   75   A(JJ) = RINF/TIM
   80   IF (J.EQ.N) GO TO 9005
        DO 85 I = JP1, N
          IJ = IJ + I - 1
          A(IJ) = Z(I)/V
   85   CONTINUE
        GO TO 9005
   90   TI = TIM
        GO TO 120
C                                  A(J,J) .GT. ZERO
   95   AL = V/A(JJ)
        TI = W(J)
        IF (MM.EQ.0) TI = TIM + V*AL
        R = TI/TIM
        A(JJ) = R*A(JJ)
        IF (R.EQ.ZERO) GO TO 125
        IF (J.EQ.N) GO TO 125
C                                  UPDATE REMAINDER OF COLUMN J
        B = AL/TI
        IF (R.GT.FOUR) GO TO 105
        DO 100 I = JP1, N
          IJ = IJ + I - 1
          Z(I) = Z(I) - V*A(IJ)
          A(IJ) = A(IJ) + B*Z(I)
  100   CONTINUE
        GO TO 115
  105   GM = TIM/TI
        DO 110 I = JP1, N
          IJ = IJ + I - 1
          Y = A(IJ)
          A(IJ) = B*Z(I) + Y*GM
          Z(I) = Z(I) - V*Y
  110   CONTINUE
  115   TIM = TI
  120 CONTINUE
  125 IF (IR.LT.0) IR = -IR
 9005 CONTINUE

      RETURN
      END

C **************************************************************************

C $RCSfile: zsrch.f,v $
C
C                      P R O J E C T:  Z X M W D
C
C ----------------------------------------------------------------------
C   PURPOSE             - Extension to IMSL routine ZXMWD
C ----------------------------------------------------------------------
CVS:
C 
C $RCSfile: zsrch.f,v $
C $Revision: 2.1 $
C $Date: 2000/04/12 22:41:57 $
C $Author: chomko $
C
C $Header: /usr/users/chomko/Project/ZXMWD/zsrch.f,v 2.1 2000/04/12 22:41:57 chomko Exp $
C-----------------------------------------------------------------------
C
C   IMSL ROUTINE NAME   - ZSRCH
C
C   COMPUTER            - VAX/SINGLE
C
C   LAST IMSL REVISION  - JANUARY 1, 1978
C
C   PURPOSE             - GENERATE POINTS IN AN N DIMENSIONAL SPACE
C
C   USAGE               - CALL ZSRCH (A,B,N,K,IP,S,M,IW,IER)
C
C   ARGUMENTS    A,B,N  - PARAMETERS WHICH DEFINE THE RECTANGULAR REGION
C                           IN N DIMENSIONAL SPACE. (INPUT)
C                           A AND B ARE VECTORS OF LENGTH N.
C                           GENERATED POINTS SATISFY
C                           A(I) .LT. S(I) .LT. B(I) FOR I=1,2,...,N.
C                           NOTE THAT IF B(I) .LT. A(I), THEN
C                           B(I) .LT. S(I) .LT. A(I).
C                K      - NUMBER OF POINTS TO BE GENERATED. (INPUT)
C                IP     - INITIALIZATION PARAMETER. (INPUT)
C                           IP MUST BE SET TO 0 FOR THE FIRST CALL.
C                           ZSRCH RESETS IP TO 1 AND RETURNS THE
C                           FIRST GENERATED POINT IN S.
C                           SUBSEQUENT CALLS SHOULD BE MADE WITH
C                           IP = 1.
C                S      - VECTOR OF LENGTH N CONTAINING THE
C                           GENERATED POINT. (OUTPUT)
C                           EACH CALL RESULTS IN THE NEXT GENERATED
C                           POINT BEING STORED IN S
C                           (THAT IS S(1),S(2),...,S(N)).
C                M      - WORK VECTOR OF LENGTH N.
C                           M MUST BE PRESERVED BETWEEN CALLS TO ZSRCH.
C                IW     - WORK VECTOR OF LENGTH 9.
C                           IW MUST BE PRESERVED BETWEEN CALLS TO ZSRCH.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER=129, ATTEMPT TO GENERATE MORE THAN
C                             K POINTS.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - UERTST,UGETIO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS      ZSRCH MAY BE USED WITH ANY NONLINEAR OPTIMIZATION
C                ROUTINE THAT REQUIRES STARTING POINTS. THE RECTANGLE
C                TO BE SEARCHED (DEFINED BY PARAMETERS A, B, AND N)
C                MUST BE DETERMINED AND THE NUMBER OF STARTING POINTS,
C                K, MUST BE CHOSEN. ONE POSSIBLE USE FOR ZSRCH WOULD
C                BE TO CALL ZSRCH TO GENERATE A POINT IN THE CHOSEN
C                RECTANGLE. THEN CALL THE NONLINEAR OPTIMIZATION
C                ROUTINE USING THIS POINT AS AN INITIAL GUESS FOR THE
C                SOLUTION. REPEAT THIS PROCESS K TIMES. THE NUMBER OF
C                ITERATIONS THAT THE OPTIMIZATION ROUTINE IS ALLOWED
C                TO PERFORM SHOULD BE QUITE SMALL (5 TO 10) DURING
C                THIS SEARCH PROCESS. THE BEST (OR BEST SEVERAL)
C                POINT(S) FOUND DURING THE SEARCH MAY BE USED AS AN
C                INITIAL GUESS TO ALLOW THE OPTIMIZATION ROUTINE TO
C                DETERMINE THE OPTIMUM MORE ACCURATELY. IN THIS
C                MANNER, AN N DIMENSIONAL RECTANGLE MAY BE
C                EFFECTIVELY SEARCHED FOR A GLOBAL OPTIMUM OF A
C                NONLINEAR FUNCTION. THE CHOICE OF K DEPENDS UPON
C                THE NONLINEARITY OF THE FUNCTION BEING OPTIMIZED.
C                ONE WITH MANY LOCAL OPTIMA REQUIRES A LARGER VALUE
C                THAN ONE WITH ONLY A FEW LOCAL OPTIMA.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
CVS:
C $Id: zsrch.f,v 2.1 2000/04/12 22:41:57 chomko Exp $
C
C                ========== HISTORY ==========
C $Log: zsrch.f,v $
C Revision 2.1  2000/04/12 22:41:57  chomko
C Bringing ALL up to Revision 2.1
C
C Revision 2.0.0.1  2000/04/12 22:39:57  chomko
C Making CVS support, v2.0.0
C
C $Header: /usr/users/chomko/Project/ZXMWD/zsrch.f,v 2.1 2000/04/12 22:41:57 chomko Exp $
C **********************************************************************
C
      SUBROUTINE ZSRCH  (A,B,N,K,IP,S,M,IW,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,K,IP,M(N),IW(9),IER
      REAL               A(N),B(N),S(N)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            L,I,NS,J,JP1,L1,NML1,L2,I1,I2,NJP1,MS,JG,JM,
     1                   NCHALF
      REAL               E
      REAL               C,RJP1,RMI,ZERO,HALF,ONE,RK
      DATA               ZERO,HALF,ONE/0.0,0.5,1.0/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  IP .EQ. 0 SIGNALS INITIAL ENTRY TO
C                                    INITIALIZE AND GENERATE FIRST
C                                    POINT
      IF (IP.GT.0) GO TO 15
      E = ONE/N
C                                  COMPUTE J, L1 AND L2 SO THAT
C                                    K = J**(N-L1) * (J+1)**L1 + L2
      RK = K
      J = RK**E
      JP1 = J+1
      IF (JP1**N.LE.K) J = JP1
      JP1 = J+1
      DO 5 L=1,N
         L1 = L
         J4PLUS=J
         NL=N-L
         J4PLUS=J4PLUS**NL
         J5PLUS=JP1
         J5PLUS=J5PLUS**L
         IDK=J4PLUS*J5PLUS
         IF (IDK .GT. K) GO TO 10
    5 CONTINUE
   10 L1 = L1-1
      NML1 = N-L1
      L2 = K-J**NML1*JP1**L1
C                                  I1 IS THE NUMBER OF POINTS THAT CAN
C                                    BE GENERATED WITHOUT ESTABLISHING
C                                    A NEW BASE FOR M
      I1 = J**(NML1-1)*JP1**L1
      IF (L2.EQ.0.AND.L1.GT.0) I1 = (J*I1)/JP1
C
C                                  I2 COUNTS(MODULO I1) THE POINTS AS
C                                    THEY ARE GENERATED
      I2 = 0
      NJP1 = NML1
      IF (L2.EQ.0) NJP1 = NML1+1
      MS = 0
      NS = MIN0(N,NJP1)
C                                  M IS MOVED TO THE NEXT DIAGONAL
C                                    AFTER EACH GROUP OF JG POINTS
      JG = JP1
      IF (NS.EQ.N) JG = J
C                                  PREVENT REJECTION OF POINTS WHEN
C                                    L2 .EQ. 0
      IF (L2.EQ.0) L2 = I1
      IP = 1
      GO TO 20
   15 CONTINUE
C                                  ENTRY TO GENERATE NEXT POINT
C                                    RESTORE LOCAL VARIABLES
      J = IW(1)
      L1 = IW(2)
      L2 = IW(3)
      I1 = IW(4)
      I2 = IW(5)
      NJP1 = IW(6)
      MS = IW(7)
      NS = IW(8)
      JG = IW(9)
      JP1 = J+1
      NML1 = N-L1
   20 CONTINUE
C                                  INCREMENT I2
      I2 = MOD(I2,I1)+1
      IF (I2.GT.1) GO TO 30
C                                  ESTABLISH BASE VALUE FOR M
      DO 25 I=1,N
   25 M(I) = 1
      MS = MS+1
      JM = J
      IF (L1.GT.0.OR.L2.GT.0) JM = JP1
C                                  IF MS .GT. JM ALL K POINTS HAVE BEEN
C                                    GENERATED
      IF (MS.GT.JM) GO TO 9000
      M(NS) = MS
      GO TO 60
C                                  MOVE M ALONG CURRENT DIAGONAL
   30 JM = J
      DO 35 I=1,N
         IF (I.EQ.NJP1) JM = JP1
         M(I) = MOD(M(I),JM)+1
   35 CONTINUE
C                                  CHECK FOR MOVE TO NEW DIAGONAL
      IF (MOD(I2,JG).NE.1) GO TO 60
      IF (NS.EQ.N) GO TO 45
C                                  MOVE M BACK TO LAST DIAGONAL
      DO 40 I=1,NS
         M(I) = M(I)-1
         IF (M(I).EQ.0) M(I) = J
   40 CONTINUE
      M(NS) = MS
C                                  MOVE M TO NEXT DIAGONAL
   45 JM = J
      DO 55 I=1,N
         IF (I.NE.NS) GO TO 50
         JM = JP1
         GO TO 55
   50    M(I) = MOD(M(I),JM)+1
         IF (M(I).GT.1) GO TO 60
   55 CONTINUE
C                                  IF THIS LOOP TERMINATES ALL K
C                                    POINTS HAVE BEEN GENERATED
      GO TO 9000
C                                  REJECT M(NS) .EQ. J+1 WHEN
C                                    I2 .GT. L2
C
   60 IF (I2.GT.L2.AND.M(NS).EQ.JP1) GO TO 20
C
C                                  MAP M INTO S
      C = ZERO
      NCHALF = NML1+1
      IF (I2.LE.L2) NCHALF = NJP1
      RJP1 = ONE/JP1
      DO 65 I=1,N
         IF (I.EQ.NCHALF) C = HALF
         RMI = M(I)
         S(I) = A(I)+(RMI-C)*(B(I)-A(I))*RJP1
   65 CONTINUE
C                                  SAVE LOCAL VARIABLES
      IW(1) = J
      IW(2) = L1
      IW(3) = L2
      IW(4) = I1
      IW(5) = I2
      IW(6) = NJP1
      IW(7) = MS
      IW(8) = NS
      IW(9) = JG
      GO TO 9005
 9000 CONTINUE
      IP = 0
      IER = 129

 9005 RETURN
      END

C ***************************************************

C ***************************************************

      SUBROUTINE ZXMWH_INI (uf,fName)

      IMPLICIT NONE
      CHARACTER*(*) fName
      INTEGER       uf

      uf=67
      OPEN(unit=uf,file=fName,status='UNKNOWN')
      WRITE(uf,100)
      WRITE(uf,101)

 100  FORMAT ('itn  f  x(5)  g(5)  ifn  f_lambda(6)  w0')
 101  FORMAT (//)

      RETURN
      END

C ***************************************************

C ***************************************************

      SUBROUTINE ZXMWH_WRITE (uf,N,itn,f,x,g,ifn)

      IMPLICIT NONE
      INTEGER  uf,N,itn,ifn
      REAL     f,x(N),g(N)

      IF (N.EQ.5) THEN
        WRITE(uf,105) itn,f,
     &                x(1),x(2),x(3),x(4),x(5),
     &                g(1),g(2),g(3),g(4),g(5),
     &                ifn
      ELSE IF (N.EQ.2) THEN
        WRITE(uf,102) itn,f,
     &                x(1),x(2),
     &                g(1),g(2),
     &                ifn
      ENDIF

 102  FORMAT(I3,' ',F8.4,' ',2(' ',F8.5),' ',2(' ',F9.5),' ',I4)
 105  FORMAT(I3,' ',F8.4,' ',5(' ',F8.5),' ',5(' ',F9.5),' ',I4)

      RETURN
      END

C ***************************************************

C ***************************************************

      SUBROUTINE DZXMWH_WRITE (lat,uf,N,itn,f,x,g,f_lam,
     &          tau865_i,v_i,w0_i,
     &          awr1,awr2,awr3,awr4,awr5,awr6,start_num,start_best)

      IMPLICIT NONE
      INTEGER  uf,N,itn,ifn,NLAMBDA,start_num,start_best
      PARAMETER (NLAMBDA=8)
      REAL*8   f,x(N),g(N),f_lam(NLAMBDA)
      REAL     tau865_i,v_i,w0_i,awr1,awr2,awr3,awr4,awr5,awr6,lat

      IF (N.EQ.5) THEN
        WRITE(uf,105) lat,itn,f,
     &     REAL(x(1)),REAL(x(2)),REAL(x(3)),REAL(x(4)),REAL(x(5)),
     &     REAL(g(1)),REAL(g(2)),REAL(g(3)),REAL(g(4)),REAL(g(5)),
     &     REAL(f_lam(1)),REAL(f_lam(2)),REAL(f_lam(3)),
     &     REAL(f_lam(4)),REAL(f_lam(5)),REAL(f_lam(6)),
     &     tau865_i,v_i,w0_i,
     &     awr1,awr2,awr3,awr4,awr5,awr6,start_num,start_best
      ELSE IF (N.EQ.2) THEN
        WRITE(uf,102) itn,f,
     &                REAL(x(1)),REAL(x(2)),
     &                REAL(g(1)),REAL(g(2)),
     &                ifn,
     &     REAL(f_lam(1)),REAL(f_lam(2)),REAL(f_lam(3)),
     &     REAL(f_lam(4)),REAL(f_lam(5)),REAL(f_lam(6))
      ENDIF

c 102  FORMAT(
c     & I3,' ',E11.4,' ',2(' ',E11.5),' ',2(' ',E12.5),' ',I4)
c 105  FORMAT(I3,' ',E11.4,' ',5(' ',E11.5),' ',5(' ',E12.5),' ',I4)

 102  FORMAT(I3,' ',F8.4,' ',2(' ',F8.5),' ',2(' ',F12.5),' ',I4,
     &                                               ' ',6(' ',F9.5))
 105  FORMAT(f8.3,I4,(' ',F12.3),' ',
     &  F16.11,3(' ',F9.4),' ',F9.6,'      ',
     &  F11.6,4(' ',F9.3),'   ',
     &  F12.3,5(' ',F9.3),'   ',
     &  3(' ',F9.5),' ',
     &  6(' ',F9.5),' ',
     &  2(' ',I4)         )

      RETURN
      END

       
C ***************************************************

C ***************************************************

      SUBROUTINE ZXMWH_END (uf)

      IMPLICIT NONE
      INTEGER       uf

      CLOSE (uf)

      RETURN
      END

C ############################################################################
C ############################################################################

C SECTION 10:

C optimisation routines for visible


c================    L-BFGS-B (version 2.1)   ==========================
 
      subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,
     +                 task, iprint,  csave, lsave, isave, dsave)
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, 
     +                 nbd(n), iwa(3*n), isave(44)
      double precision f, factr, pgtol, x(n), l(n), u(n), g(n),
     +                 wa(2*m*n+4*n+12*m*m+12*m), dsave(29)
 
c     ************
c
c     Subroutine setulb
c
c     This subroutine partitions the working arrays wa and iwa, and 
c       then uses the limited memory BFGS method to solve the bound
c       constrained optimization problem by calling mainlb.
c       (The direct method will be used in the subspace minimization.)
c
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x is an approximation to the solution.
c       On exit x is the current approximation.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound on x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound on x.
c       On exit u is unchanged.
c
c     nbd is an integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds, and
c                3 if x(i) has only an upper bound.
c       On exit nbd is unchanged.
c
c     f is a double precision variable.
c       On first entry f is unspecified.
c       On final exit f is the value of the function at x.
c
c     g is a double precision array of dimension n.
c       On first entry g is unspecified.
c       On final exit g is the value of the gradient at x.
c
c     factr is a double precision variable.
c       On entry factr >= 0 is specified by the user.  The iteration
c         will stop when
c
c         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c         where epsmch is the machine precision, which is automatically
c         generated by the code. Typical values for factr: 1.d+12 for
c         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
c         high accuracy.
c       On exit factr is unchanged.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.   
c       On exit pgtol is unchanged.
c
c     wa is a double precision working array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax.
c
c     iwa is an integer working array of length 3nmax.
c
c     task is a working string of characters of length 60 indicating
c       the current job when entering and quitting this subroutine.
c
c     iprint is an integer variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave is a working string of characters of length 60.
c
c     lsave is a logical working array of dimension 4.
c       On exit with 'task' = NEW_X, the following information is 
c                                                             available:
c         If lsave(1) = .true.  then  the initial X has been replaced by
c                                     its projection in the feasible set;
c         If lsave(2) = .true.  then  the problem is constrained;
c         If lsave(3) = .true.  then  each variable has upper and lower
c                                     bounds;
c
c     isave is an integer working array of dimension 44.
c       On exit with 'task' = NEW_X, the following information is 
c                                                             available:
c         isave(22) = the total number of intervals explored in the 
c                         search of Cauchy points;
c         isave(26) = the total number of skipped BFGS updates before 
c                         the current iteration;
c         isave(30) = the number of current iteration;
c         isave(31) = the total number of BFGS updates prior the current
c                         iteration;
c         isave(33) = the number of intervals explored in the search of
c                         Cauchy point in the current iteration;
c         isave(34) = the total number of function and gradient 
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         if isave(37) = 0  then the subspace argmin is within the box;
c         if isave(37) = 1  then the subspace argmin is beyond the box;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints in the current
c                         iteration;
c         n + 1 - isave(40) = the number of variables leaving the set of
c                           active constraints in the current iteration;
c         isave(41) = the number of variables entering the set of active
c                         constraints in the current iteration.
c
c     dsave is a double precision working array of dimension 29.
c       On exit with 'task' = NEW_X, the following information is
c                                                             available:
c         dsave(1) = current 'theta' in the BFGS matrix;
c         dsave(2) = f(x) in the previous iteration;
c         dsave(3) = factr*epsmch;
c         dsave(4) = 2-norm of the line search direction vector;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(7) = the accumulated time spent on searching for
c                                                         Cauchy points;
c         dsave(8) = the accumulated time spent on
c                                                 subspace minimization;
c         dsave(9) = the accumulated time spent on line search;
c         dsave(11) = the slope of the line search function at
c                                  the current point of line search;
c         dsave(12) = the maximum relative step length imposed in
c                                                           line search;
c         dsave(13) = the infinity norm of the projected gradient;
c         dsave(14) = the relative step length in the line search;
c         dsave(15) = the slope of the line search function at
c                                 the starting point of the line search;
c         dsave(16) = the square of the 2-norm of the line search
c                                                      direction vector.
c
c     Subprograms called:
c
c       L-BFGS-B Library ... mainlb.    
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
c       limited memory FORTRAN code for solving bound constrained
c       optimization problems'', Tech. Report, NAM-11, EECS Department,
c       Northwestern University, 1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer   l1,l2,l3,lws,lr,lz,lt,ld,lsg,lwa,lyg,
     +          lsgo,lwy,lsy,lss,lyy,lwt,lwn,lsnd,lygo

      if (task .eq. 'START') then
         isave(1)  = m*n
         isave(2)  = m**2
         isave(3)  = 4*m**2
         isave(4)  = 1
         isave(5)  = isave(4)  + isave(1)
         isave(6)  = isave(5)  + isave(1)
         isave(7)  = isave(6)  + isave(2)
         isave(8)  = isave(7)  + isave(2)
         isave(9)  = isave(8)  + isave(2)
         isave(10) = isave(9)  + isave(2)
         isave(11) = isave(10) + isave(3)
         isave(12) = isave(11) + isave(3)
         isave(13) = isave(12) + n
         isave(14) = isave(13) + n
         isave(15) = isave(14) + n
         isave(16) = isave(15) + n
         isave(17) = isave(16) + 8*m
         isave(18) = isave(17) + m
         isave(19) = isave(18) + m
         isave(20) = isave(19) + m   
      endif
      l1   = isave(1)
      l2   = isave(2)
      l3   = isave(3)
      lws  = isave(4)
      lwy  = isave(5)
      lsy  = isave(6)
      lss  = isave(7)
      lyy  = isave(8)
      lwt  = isave(9)
      lwn  = isave(10)
      lsnd = isave(11)
      lz   = isave(12)
      lr   = isave(13)
      ld   = isave(14)
      lt   = isave(15)
      lwa  = isave(16)
      lsg  = isave(17)
      lsgo = isave(18)
      lyg  = isave(19)
      lygo = isave(20)

      call mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol,
     +  wa(lws),wa(lwy),wa(lsy),wa(lss),wa(lyy),wa(lwt),
     +  wa(lwn),wa(lsnd),wa(lz),wa(lr),wa(ld),wa(lt),
     +  wa(lwa),wa(lsg),wa(lsgo),wa(lyg),wa(lygo),
     +  iwa(1),iwa(n+1),iwa(2*n+1),task,iprint,
     +  csave,lsave,isave(22),dsave)

      return

      end

c======================= The end of setulb =============================
 
      subroutine mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy,
     +                  sy, ss, yy, wt, wn, snd, z, r, d, t, wa, sg,
     +                  sgo, yg, ygo, index, iwhere, indx2, task,
     +                  iprint, csave, lsave, isave, dsave)
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, nbd(n), index(n),
     +                 iwhere(n), indx2(n), isave(23)
      double precision f, factr, pgtol,
     +                 x(n), l(n), u(n), g(n), z(n), r(n), d(n), t(n), 
     +                 wa(8*m), sg(m), sgo(m), yg(m), ygo(m), 
     +                 ws(n, m), wy(n, m), sy(m, m), ss(m, m), yy(m, m),
     +                 wt(m, m), wn(2*m, 2*m), snd(2*m, 2*m), dsave(29)

c     ************
c
c     Subroutine mainlb
c
c     This subroutine solves bound constrained optimization problems by
c       using the compact formula of the limited memory BFGS updates.
c       
c     n is an integer variable.
c       On entry n is the number of variables.
c       On exit n is unchanged.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric
c          corrections allowed in the limited memory matrix.
c       On exit m is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x is an approximation to the solution.
c       On exit x is the current approximation.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound of x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound of x.
c       On exit u is unchanged.
c
c     nbd is an integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds,
c                3 if x(i) has only an upper bound.
c       On exit nbd is unchanged.
c
c     f is a double precision variable.
c       On first entry f is unspecified.
c       On final exit f is the value of the function at x.
c
c     g is a double precision array of dimension n.
c       On first entry g is unspecified.
c       On final exit g is the value of the gradient at x.
c
c     factr is a double precision variable.
c       On entry factr >= 0 is specified by the user.  The iteration
c         will stop when
c
c         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c         where epsmch is the machine precision, which is automatically
c         generated by the code.
c       On exit factr is unchanged.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       On exit pgtol is unchanged.
c
c     ws, wy, sy, and wt are double precision working arrays used to
c       store the following information defining the limited memory
c          BFGS matrix:
c          ws, of dimension n x m, stores S, the matrix of s-vectors;
c          wy, of dimension n x m, stores Y, the matrix of y-vectors;
c          sy, of dimension m x m, stores S'Y;
c          ss, of dimension m x m, stores S'S;
c          yy, of dimension m x m, stores Y'Y;
c          wt, of dimension m x m, stores the Cholesky factorization
c                                  of (theta*S'S+LD^(-1)L'); see eq.
c                                  (2.26) in [3].
c
c     wn is a double precision working array of dimension 2m x 2m
c       used to store the LEL^T factorization of the indefinite matrix
c                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c
c       where     E = [-I  0]
c                     [ 0  I]
c
c     snd is a double precision working array of dimension 2m x 2m
c       used to store the lower triangular part of
c                 N = [Y' ZZ'Y   L_a'+R_z']
c                     [L_a +R_z  S'AA'S   ]
c            
c     z(n),r(n),d(n),t(n),wa(8*m) are double precision working arrays.
c       z is used at different times to store the Cauchy point and
c       the Newton point.
c
c     sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays. 
c
c     index is an integer working array of dimension n.
c       In subroutine freev, index is used to store the free and fixed
c          variables at the Generalized Cauchy Point (GCP).
c
c     iwhere is an integer working array of dimension n used to record
c       the status of the vector x for GCP computation.
c       iwhere(i)=0 or -3 if x(i) is free and has bounds,
c                 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
c                 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
c                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
c                -1       if x(i) is always free, i.e., no bounds on it.
c
c     indx2 is an integer working array of dimension n.
c       Within subroutine cauchy, indx2 corresponds to the array iorder.
c       In subroutine freev, a list of variables entering and leaving
c       the free set is stored in indx2, and it is passed on to
c       subroutine formk with this information.
c
c     task is a working string of characters of length 60 indicating
c       the current job when entering and leaving this subroutine.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave is a working string of characters of length 60.
c
c     lsave is a logical working array of dimension 4.
c
c     isave is an integer working array of dimension 23.
c
c     dsave is a double precision working array of dimension 29.
c
c
c     Subprograms called
c
c       L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk, 
c
c        errclb, prn1lb, prn2lb, prn3lb, active, projgr,
c
c        freev, cmprlb, matupd, formt.
c
c       Minpack2 Library ... timer, dpmeps.
c
c       Linpack Library ... dcopy, ddot.
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c       Subroutines for Large Scale Bound Constrained Optimization''
c       Tech. Report, NAM-11, EECS Department, Northwestern University,
c       1994.
c 
c       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
c       Quasi-Newton Matrices and their use in Limited Memory Methods'',
c       Mathematical Programming 63 (1994), no. 4, pp. 129-156.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      logical          prjctd,cnstnd,boxed,updatd,wrk
      character*3      word
      integer          i,k,nintol,itfile,iback,nskip,
     +                 head,col,iter,itail,iupdat,
     +                 nint,nfgv,info,ifun,
     +                 iword,nfree,nact,ileave,nenter,
     &                 iRCerr
      double precision theta,fold,ddot,dr,rr,tol,dpmeps,
     +                 xstep,sbgnrm,ddum,dnorm,dtd,epsmch,
     +                 cpu1,cpu2,cachyt,sbtime,lnscht,time1,time2,
     +                 gd,gdold,stp,stpmx,time
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
      
      if (task .eq. 'START') then

         call timer(time1)

c        Generate the current machine precision.

         epsmch = dpmeps()
c         print*,epsmch

c        Initialize counters and scalars when task='START'.

c           for the limited memory BFGS matrices:
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
 
c           for operation counts:
         iter   = 0
         nfgv   = 0
         nint   = 0
         nintol = 0
         nskip  = 0
         nfree  = n

c           for stopping tolerance:
         tol = factr*epsmch
c         print*,tol

c           for measuring running time:
         cachyt = 0
         sbtime = 0
         lnscht = 0
 
c           'word' records the status of subspace solutions.
         word = '---'

c           'info' records the termination information.
         info = 0

         if (iprint .ge. 1) then
c                                open a summary file 'iterate.dat'
            open (8, file = 'iterate.dat', status = 'unknown')
            itfile = 8
         endif            

c        Check the input arguments for errors.

         call errclb(n,m,factr,l,u,nbd,task,info,k)
         if (task(1:5) .eq. 'ERROR') then
            call prn3lb(n,x,f,task,iprint,info,itfile,
     +                  iter,nfgv,nintol,nskip,nact,sbgnrm,
     +                  zero,nint,word,iback,stp,xstep,k,
     +                  cachyt,sbtime,lnscht)
            return
         endif

         call prn1lb(n,m,l,u,x,iprint,itfile,epsmch)
 
c        Initialize iwhere & project x onto the feasible set.
 
         call active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed) 

c        The end of the initialization.

      else
c          restore local variables.

         prjctd = lsave(1)
         cnstnd = lsave(2)
         boxed  = lsave(3)
         updatd = lsave(4)

         nintol = isave(1)
         itfile = isave(3)
         iback  = isave(4)
         nskip  = isave(5)
         head   = isave(6)
         col    = isave(7)
         itail  = isave(8)
         iter   = isave(9)
         iupdat = isave(10)
         nint   = isave(12)
         nfgv   = isave(13)
         info   = isave(14)
         ifun   = isave(15)
         iword  = isave(16)
         nfree  = isave(17)
         nact   = isave(18)
         ileave = isave(19)
         nenter = isave(20)

         theta  = dsave(1)
         fold   = dsave(2)
         tol    = dsave(3)
         dnorm  = dsave(4)
         epsmch = dsave(5)
         cpu1   = dsave(6)
         cachyt = dsave(7)
         sbtime = dsave(8)
         lnscht = dsave(9)
         time1  = dsave(10)
         gd     = dsave(11)
         stpmx  = dsave(12)
         sbgnrm = dsave(13)
         stp    = dsave(14)
         gdold  = dsave(15)
         dtd    = dsave(16)
   
c        After returning from the driver go to the point where execution
c        is to resume.

         if (task(1:5) .eq. 'FG_LN') goto 666
         if (task(1:5) .eq. 'NEW_X') goto 777
         if (task(1:5) .eq. 'FG_ST') goto 111
         if (task(1:4) .eq. 'STOP') then
            if (task(7:9) .eq. 'CPU') then
c                                          restore the previous iterate.
               call dcopy(n,t,1,x,1)
               call dcopy(n,r,1,g,1)
               f = fold
            endif
            goto 999
         endif
      endif 

c     Compute f0 and g0.

      task = 'FG_START' 
c          return to the driver to calculate f and g; reenter at 111.
      goto 1000
 111  continue
      nfgv = 1
 
c     Compute the infinity norm of the (-) projected gradient.
 
      call projgr(n,l,u,nbd,x,g,sbgnrm)
  
      if (iprint .ge. 1) then
         write (6,1002) iter,f,sbgnrm
         write (itfile,1003) iter,nfgv,sbgnrm,f
      endif
      if (sbgnrm .le. pgtol) then
c                                terminate the algorithm.
         task = 'CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL'
         goto 999
      endif 
 
c ----------------- the beginning of the loop --------------------------
 
 222  continue
      if (iprint .ge. 99) write (6,1001) iter + 1
      iword = -1
c
      if (.not. cnstnd .and. col .gt. 0) then 
c                                            skip the search for GCP.
         call dcopy(n,x,1,z,1)
         wrk = updatd
         nint = 0
         goto 333
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the Generalized Cauchy Point (GCP).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call timer(cpu1) 
      call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z,
     +            m,wy,ws,sy,wt,theta,col,head,
     +            wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nint,
     +            sg,yg,iprint,sbgnrm,info,epsmch,iRCerr)
      if (info .ne. 0) then 
         IF(iRCerr.NE.0) GOTO 999
c         singular triangular system detected; refresh the lbfgs memory.
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         cachyt = cachyt + cpu2 - cpu1
         goto 222
      endif
      call timer(cpu2) 
      cachyt = cachyt + cpu2 - cpu1
      nintol = nintol + nint

c     Count the entering and leaving variables for iter > 0; 
c     find the index set of free and active variables at the GCP.

      call freev(n,nfree,index,nenter,ileave,indx2,
     +           iwhere,wrk,updatd,cnstnd,iprint,iter)

      nact = n - nfree
 
 333  continue
 
c     If there are no free variables or B=theta*I, then
c                                        skip the subspace minimization.
 
      if (nfree .eq. 0 .or. col .eq. 0) goto 555
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subspace minimization.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call timer(cpu1) 

c     Form  the LEL^T factorization of the indefinite
c       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c       where     E = [-I  0]
c                     [ 0  I]

      if (wrk) call formk(n,nfree,index,nenter,ileave,indx2,iupdat,
     +                 updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
      if (info .ne. 0) then
c          nonpositive definiteness in Cholesky factorization;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1006)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         sbtime = sbtime + cpu2 - cpu1 
         goto 222
      endif 

c        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
c                                                   from 'cauchy').
      call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index,
     +           theta,col,head,nfree,cnstnd,info)
      if (info .ne. 0) goto 444
c       call the direct method.
      call subsm(n,m,nfree,index,l,u,nbd,z,r,ws,wy,theta,
     +           col,head,iword,wa,wn,iprint,info)
 444  continue
      if (info .ne. 0) then 
c          singular triangular system detected;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         sbtime = sbtime + cpu2 - cpu1 
         goto 222
      endif
 
      call timer(cpu2) 
      sbtime = sbtime + cpu2 - cpu1 
 555  continue
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Line search and optimality tests.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c     Generate the search direction d:=z-x.

      do 40 i = 1, n
         d(i) = z(i) - x(i)
  40  continue
      call timer(cpu1) 
 666  continue
      call lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm,
     +            dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,
     +            boxed,cnstnd,csave,isave(22),dsave(17))
      if (info .ne. 0 .or. iback .ge. 20) then
c          restore the previous iterate.
         call dcopy(n,t,1,x,1)
         call dcopy(n,r,1,g,1)
         f = fold
         if (col .eq. 0) then
c             abnormal termination.
            if (info .eq. 0) then
               info = -9
c                restore the actual number of f and g evaluations etc.
               nfgv = nfgv - 1
               ifun = ifun - 1
               iback = iback - 1
            endif
            task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
            iter = iter + 1
            goto 999
         else
c             refresh the lbfgs memory and restart the iteration.
            if(iprint .ge. 1) write (6, 1008)
            if (info .eq. 0) nfgv = nfgv - 1
            info   = 0
            col    = 0
            head   = 1
            theta  = one
            iupdat = 0
            updatd = .false.
            task   = 'RESTART_FROM_LNSRCH'
            call timer(cpu2)
            lnscht = lnscht + cpu2 - cpu1
            goto 222
         endif
      else if (task(1:5) .eq. 'FG_LN') then
c          return to the driver for calculating f and g; reenter at 666.
         goto 1000
      else 
c          calculate and print out the quantities related to the new X.
         call timer(cpu2) 
         lnscht = lnscht + cpu2 - cpu1
         iter = iter + 1
 
c        Compute the infinity norm of the projected (-)gradient.
 
         call projgr(n,l,u,nbd,x,g,sbgnrm)
 
c        Print iteration information.

         call prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,
     +               sbgnrm,nint,word,iword,iback,stp,xstep)
         goto 1000
      endif
 777  continue

c     Test for termination.

      if (sbgnrm .le. pgtol) then
c                                terminate the algorithm.
         task = 'CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL'
         goto 999
      endif 

      ddum = max(abs(fold), abs(f), one)
      if ((fold - f) .le. tol*ddum) then
c                                        terminate the algorithm.
         task = 'CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH'
         if (iback .ge. 10) info = -5
c           i.e., to issue a warning if iback>10 in the line search.
         goto 999
      endif 

c     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
 
      do 42 i = 1, n
         r(i) = g(i) - r(i)
  42  continue
      rr = ddot(n,r,1,r,1)
      if (stp .eq. one) then  
         dr = gd - gdold
         ddum = -gdold
      else
         dr = (gd - gdold)*stp
         call dscal(n,stp,d,1)
         ddum = -gdold*stp
      endif
 
      if (dr .le. epsmch*ddum) then
c                            skip the L-BFGS update.
         nskip = nskip + 1
         updatd = .false.
         if (iprint .ge. 1) write (6,1004) dr, ddum
         goto 888
      endif 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update the L-BFGS matrix.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      updatd = .true.
      iupdat = iupdat + 1

c     Update matrices WS and WY and form the middle matrix in B.

      call matupd(n,m,ws,wy,sy,ss,d,r,itail,
     +            iupdat,col,head,theta,rr,dr,stp,dtd)

c     Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
c        Store T in the upper triangular of the array wt;
c        Cholesky factorize T to J*J' with
c           J' stored in the upper triangular of wt.

      call formt(m,wt,sy,ss,col,theta,info)
 
      if (info .ne. 0) then 
c          nonpositive definiteness in Cholesky factorization;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1007)
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         goto 222
      endif

c     Now the inverse of the middle matrix in B is

c       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
c       [ -L*D^(-1/2)   J ] [  0        J'          ]

 888  continue
 
c -------------------- the end of the loop -----------------------------
 
      goto 222
 999  continue
      call timer(time2)
      time = time2 - time1
      call prn3lb(n,x,f,task,iprint,info,itfile,
     +            iter,nfgv,nintol,nskip,nact,sbgnrm,
     +            time,nint,word,iback,stp,xstep,k,
     +            cachyt,sbtime,lnscht)
 1000 continue

c     Save local variables.

      lsave(1)  = prjctd
      lsave(2)  = cnstnd
      lsave(3)  = boxed
      lsave(4)  = updatd

      isave(1)  = nintol 
      isave(3)  = itfile 
      isave(4)  = iback 
      isave(5)  = nskip 
      isave(6)  = head 
      isave(7)  = col 
      isave(8)  = itail 
      isave(9)  = iter 
      isave(10) = iupdat 
      isave(12) = nint 
      isave(13) = nfgv 
      isave(14) = info 
      isave(15) = ifun 
      isave(16) = iword 
      isave(17) = nfree 
      isave(18) = nact 
      isave(19) = ileave 
      isave(20) = nenter 

      dsave(1)  = theta 
      dsave(2)  = fold 
      dsave(3)  = tol 
      dsave(4)  = dnorm 
      dsave(5)  = epsmch 
      dsave(6)  = cpu1 
      dsave(7)  = cachyt 
      dsave(8)  = sbtime 
      dsave(9)  = lnscht 
      dsave(10) = time1 
      dsave(11) = gd 
      dsave(12) = stpmx 
      dsave(13) = sbgnrm
      dsave(14) = stp
      dsave(15) = gdold
      dsave(16) = dtd  

 1001 format (//,'ITERATION ',i5)
 1002 format
     +  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 1003 format (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x,
     +        1p,2(1x,d10.3))
 1004 format ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,' BFGS update SKIPPED')
 1005 format (/, 
     +' Singular triangular system detected;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1006 format (/, 
     +' Nonpositive definiteness in Cholesky factorization in formk;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1007 format (/, 
     +' Nonpositive definiteness in Cholesky factorization in formt;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1008 format (/, 
     +' Bad direction in the line search;',/,
     +'   refresh the lbfgs memory and restart the iteration.')

      return   

      end
 
c======================= The end of mainlb =============================

      subroutine active(n, l, u, nbd, x, iwhere, iprint,
     +                  prjctd, cnstnd, boxed)

      logical          prjctd, cnstnd, boxed
      integer          n, iprint, nbd(n), iwhere(n)
      double precision x(n), l(n), u(n)

c     ************
c
c     Subroutine active
c
c     This subroutine initializes iwhere and projects the initial x to
c       the feasible set if necessary.
c
c     iwhere is an integer array of dimension n.
c       On entry iwhere is unspecified.
c       On exit iwhere(i)=-1  if x(i) has no bounds
c                         3   if l(i)=u(i)
c                         0   otherwise.
c       In cauchy, iwhere is given finer gradations.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          nbdd,i
      double precision zero
      parameter        (zero=0.0d0)

c     Initialize nbdd, prjctd, cnstnd and boxed.

      nbdd = 0
      prjctd = .false.
      cnstnd = .false.
      boxed = .true.

c     Project the initial x to the easible set if necessary.

      do 10 i = 1, n
         if (nbd(i) .gt. 0) then
            if (nbd(i) .le. 2 .and. x(i) .le. l(i)) then
               if (x(i) .lt. l(i)) then
                  prjctd = .true.
                  x(i) = l(i)
               endif
               nbdd = nbdd + 1
            else if (nbd(i) .ge. 2 .and. x(i) .ge. u(i)) then
               if (x(i) .gt. u(i)) then
                  prjctd = .true.
                  x(i) = u(i)
               endif
               nbdd = nbdd + 1
            endif
         endif
  10  continue

c     Initialize iwhere and assign values to cnstnd and boxed.

      do 20 i = 1, n
         if (nbd(i) .ne. 2) boxed = .false.
         if (nbd(i) .eq. 0) then
c                                this variable is always free
            iwhere(i) = -1

c           otherwise set x(i)=mid(x(i), u(i), l(i)).
         else
            cnstnd = .true.
            if (nbd(i) .eq. 2 .and. u(i) - l(i) .le. zero) then
c                   this variable is always fixed
               iwhere(i) = 3
            else 
               iwhere(i) = 0
            endif
         endif
  20  continue

      if (iprint .ge. 0) then
         if (prjctd) write (6,*)
     +   'The initial X is infeasible.  Restart with its projection.'
         if (.not. cnstnd)
     +      write (6,*) 'This problem is unconstrained.'
      endif

      if (iprint .gt. 0) write (6,1001) nbdd

 1001 format (/,'At X0 ',i9,' variables are exactly at the bounds') 

      return

      end

c======================= The end of active =============================
 
      subroutine bmv(m, sy, wt, col, v, p, info)

      integer m, col, info
      double precision sy(m, m), wt(m, m), v(2*col), p(2*col)

c     ************
c
c     Subroutine bmv
c
c     This subroutine computes the product of the 2m x 2m middle matrix 
c       in the compact L-BFGS formula of B and a 2m vector v;  
c       it returns the product in p.
c       
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     sy is a double precision array of dimension m x m.
c       On entry sy specifies the matrix S'Y.
c       On exit sy is unchanged.
c
c     wt is a double precision array of dimension m x m.
c       On entry wt specifies the upper triangular matrix J' which is 
c         the Cholesky factor of (thetaS'S+LD^(-1)L').
c       On exit wt is unchanged.
c
c     col is an integer variable.
c       On entry col specifies the number of s-vectors (or y-vectors)
c         stored in the compact L-BFGS formula.
c       On exit col is unchanged.
c
c     v is a double precision array of dimension 2col.
c       On entry v specifies vector v.
c       On exit v is unchanged.
c
c     p is a double precision array of dimension 2col.
c       On entry p is unspecified.
c       On exit p is the product Mv.
c
c     info is an integer variable.
c       On entry info is unspecified.
c       On exit info = 0       for normal return,
c                    = nonzero for abnormal return when the system
c                                to be solved by dtrsl is singular.
c
c     Subprograms called:
c
c       Linpack ... dtrsl.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer          i,k,i2
      double precision sum
 
      if (col .eq. 0) return
 
c     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
c                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

c       solve Jp2=v2+LD^(-1)v1.
      p(col + 1) = v(col + 1)
      do 20 i = 2, col
         i2 = col + i
         sum = 0.0d0
         do 10 k = 1, i - 1
            sum = sum + sy(i,k)*v(k)/sy(k,k)
  10     continue
         p(i2) = v(i2) + sum
  20  continue  
c     Solve the triangular system
      call dtrsl(wt,m,col,p(col+1),11,info)
      if (info .ne. 0) return
 
c       solve D^(1/2)p1=v1.
      do 30 i = 1, col
         p(i) = v(i)/sqrt(sy(i,i))
  30  continue 
 
c     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
c                    [  0         J'           ] [ p2 ]   [ p2 ]. 
 
c       solve J^Tp2=p2. 
      call dtrsl(wt,m,col,p(col+1),01,info)
      if (info .ne. 0) return
 
c       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
c                 =-D^(-1/2)p1+D^(-1)L'p2.  
      do 40 i = 1, col
         p(i) = -p(i)/sqrt(sy(i,i))
  40  continue
      do 60 i = 1, col
         sum = 0.d0
         do 50 k = i + 1, col
            sum = sum + sy(k,i)*p(col+k)/sy(i,i)
  50     continue
         p(i) = p(i) + sum
  60  continue

      return

      end

c======================== The end of bmv ===============================

      subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, 
     +                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, 
     +                  v, nint, sg, yg, iprint, sbgnrm, info, epsmch,
     &                  iRCerr)
      
      integer          n, m, head, col, nint, iprint, info, 
     +                 nbd(n), iorder(n), iwhere(n), iRCerr
      double precision theta, epsmch,
     +                 x(n), l(n), u(n), g(n), t(n), d(n), xcp(n),
     +                 sg(m), yg(m), wy(n, col), ws(n, col), sy(m, m),
     +                 wt(m, m), p(2*m), c(2*m), wbp(2*m), v(2*m)

c     ************
c
c     Subroutine cauchy
c
c     For given x, l, u, g (with sbgnrm > 0), and a limited memory
c       BFGS matrix B defined in terms of matrices WY, WS, WT, and
c       scalars head, col, and theta, this subroutine computes the
c       generalized Cauchy point (GCP), defined as the first local
c       minimizer of the quadratic
c
c                  Q(x + s) = g's + 1/2 s'Bs
c
c       along the projected gradient direction P(x-tg,l,u).
c       The routine returns the GCP in xcp. 
c       
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x is the starting point for the GCP computation.
c       On exit x is unchanged.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound of x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound of x.
c       On exit u is unchanged.
c
c     nbd is an integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds, and
c                3 if x(i) has only an upper bound. 
c       On exit nbd is unchanged.
c
c     g is a double precision array of dimension n.
c       On entry g is the gradient of f(x).  g must be a nonzero vector.
c       On exit g is unchanged.
c
c     iorder is an integer working array of dimension n.
c       iorder will be used to store the breakpoints in the piecewise
c       linear path and free variables encountered. On exit,
c         iorder(1),...,iorder(nleft) are indices of breakpoints
c                                which have not been encountered; 
c         iorder(nleft+1),...,iorder(nbreak) are indices of
c                                     encountered breakpoints; and
c         iorder(nfree),...,iorder(n) are indices of variables which
c                 have no bound constraits along the search direction.
c
c     iwhere is an integer array of dimension n.
c       On entry iwhere indicates only the permanently fixed (iwhere=3)
c       or free (iwhere= -1) components of x.
c       On exit iwhere records the status of the current x variables.
c       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
c                 0   if x(i) is free and has bounds, and is moved
c                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
c                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
c                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
c                 -1  if x(i) is always free, i.e., it has no bounds.
c
c     t is a double precision working array of dimension n. 
c       t will be used to store the break points.
c
c     d is a double precision array of dimension n used to store
c       the Cauchy direction P(x-tg)-x.
c
c     xcp is a double precision array of dimension n used to return the
c       GCP on exit.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections 
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     ws, wy, sy, and wt are double precision arrays.
c       On entry they store information that defines the
c                             limited memory BFGS matrix:
c         ws(n,m) stores S, a set of s-vectors;
c         wy(n,m) stores Y, a set of y-vectors;
c         sy(m,m) stores S'Y;
c         wt(m,m) stores the
c                 Cholesky factorization of (theta*S'S+LD^(-1)L').
c       On exit these arrays are unchanged.
c
c     theta is a double precision variable.
c       On entry theta is the scaling factor specifying B_0 = theta I.
c       On exit theta is unchanged.
c
c     col is an integer variable.
c       On entry col is the actual number of variable metric
c         corrections stored so far.
c       On exit col is unchanged.
c
c     head is an integer variable.
c       On entry head is the location of the first s-vector (or y-vector)
c         in S (or Y).
c       On exit col is unchanged.
c
c     p is a double precision working array of dimension 2m.
c       p will be used to store the vector p = W^(T)d.
c
c     c is a double precision working array of dimension 2m.
c       c will be used to store the vector c = W^(T)(xcp-x).
c
c     wbp is a double precision working array of dimension 2m.
c       wbp will be used to store the row of W corresponding
c         to a breakpoint.
c
c     v is a double precision working array of dimension 2m.
c
c     nint is an integer variable.
c       On exit nint records the number of quadratic segments explored
c         in searching for the GCP.
c
c     sg and yg are double precision arrays of dimension m.
c       On entry sg  and yg store S'g and Y'g correspondingly.
c       On exit they are unchanged. 
c 
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     sbgnrm is a double precision variable.
c       On entry sbgnrm is the norm of the projected gradient at x.
c       On exit sbgnrm is unchanged.
c
c     info is an integer variable.
c       On entry info is 0.
c       On exit info = 0       for normal return,
c                    = nonzero for abnormal return when the the system
c                              used in routine bmv is singular.
c
c     Subprograms called:
c 
c       L-BFGS-B Library ... hpsolb, bmv.
c
c       Linpack ... dscal dcopy, daxpy.
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c       Subroutines for Large Scale Bound Constrained Optimization''
c       Tech. Report, NAM-11, EECS Department, Northwestern University,
c       1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      logical          xlower,xupper,bnded
      integer          i,j,col2,nfree,nbreak,pointr,
     +                 ibp,nleft,ibkmin,iter
      double precision f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin,
     +                 tu,tl,wmc,wmp,wmw,ddot,tj,tj0,neggi,sbgnrm,
     +                 f2_org
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
 
c     Check the status of the variables, reset iwhere(i) if necessary;
c       compute the Cauchy direction d and the breakpoints t; initialize
c       the derivative f1 and the vector p = W'd (for theta = 1).
 
      if (sbgnrm .le. zero) then
         if (iprint .ge. 0) write (6,*) 'Subgnorm = 0.  GCP = X.'
         call dcopy(n,x,1,xcp,1)
         return
      endif 
      bnded = .true.
      nfree = n + 1
      nbreak = 0
      ibkmin = 0
      bkmin = zero
      col2 = 2*col
      f1 = zero
      if (iprint .ge. 99) write (6,3010)

c     We set p to zero and build it up as we determine d.

      do 20 i = 1, col2
         p(i) = zero
  20  continue 

c     In the following loop we determine for each variable its bound
c        status and its breakpoint, and update p accordingly.
c        Smallest breakpoint is identified.

      do 50 i = 1, n 
         neggi = -g(i)      
         if (iwhere(i) .ne. 3 .and. iwhere(i) .ne. -1) then
c             if x(i) is not a constant and has bounds,
c             compute the difference between x(i) and its bounds.
            if (nbd(i) .le. 2) tl = x(i) - l(i)
            if (nbd(i) .ge. 2) tu = u(i) - x(i)

c           If a variable is close enough to a bound
c             we treat it as at bound.
            xlower = nbd(i) .le. 2 .and. tl .le. zero
            xupper = nbd(i) .ge. 2 .and. tu .le. zero

c              reset iwhere(i).
            iwhere(i) = 0
            if (xlower) then
               if (neggi .le. zero) iwhere(i) = 1
            else if (xupper) then
               if (neggi .ge. zero) iwhere(i) = 2
            else
               if (abs(neggi) .le. zero) iwhere(i) = -3
            endif
         endif 
         pointr = head
         if (iwhere(i) .ne. 0 .and. iwhere(i) .ne. -1) then
            d(i) = zero
         else
            d(i) = neggi
            f1 = f1 - neggi*neggi
c             calculate p := p - W'e_i* (g_i).
            do 40 j = 1, col
               p(j) = p(j) +  wy(i,pointr)* neggi
               p(col + j) = p(col + j) + ws(i,pointr)*neggi
               pointr = mod(pointr,m) + 1
  40        continue 
            if (nbd(i) .le. 2 .and. nbd(i) .ne. 0
     +                        .and. neggi .lt. zero) then
c                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tl/(-neggi)
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else if (nbd(i) .ge. 2 .and. neggi .gt. zero) then
c                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tu/neggi
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else
c                x(i) + d(i) is not bounded.
               nfree = nfree - 1
               iorder(nfree) = i
               if (abs(neggi) .gt. zero) bnded = .false.
            endif
         endif
  50  continue 
 
c     The indices of the nonzero components of d are now stored
c       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
c       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
 
      if (theta .ne. one) then
c                   complete the initialization of p for theta not= one.
         call dscal(col,theta,p(col+1),1)
      endif
 
c     Initialize GCP xcp = x.

      call dcopy(n,x,1,xcp,1)

      if (nbreak .eq. 0 .and. nfree .eq. n + 1) then
c                  is a zero vector, return with the initial xcp as GCP.
         if (iprint .gt. 100) write (6,1010) (xcp(i), i = 1, n)
         return
      endif    
 
c     Initialize c = W'(xcp - x) = 0.
  
      do 60 j = 1, col2
         c(j) = zero
  60  continue 
 
c     Initialize derivative f2.
 
      f2 =  -theta*f1 
      f2_org  =  f2
      if (col .gt. 0) then
         call bmv(m,sy,wt,col,p,v,info)
         if (info .ne. 0) return
         f2 = f2 - ddot(col2,v,1,p,1)
      endif
      IF (DABS(f2).LE.1.D-16) THEN
         info=1
         RETURN
      ENDIF
      dtm = -f1/f2
      tsum = zero
      nint = 1
      if (iprint .ge. 99) 
     +   write (6,*) 'There are ',nbreak,'  breakpoints '
 
c     If there are no breakpoints, locate the GCP and return. 
 
      if (nbreak .eq. 0) goto 888
             
      nleft = nbreak
      iter = 1
 
 
      tj = zero
 
c------------------- the beginning of the loop -------------------------
 
 777  continue
 
c     Find the next smallest breakpoint;
c       compute dt = t(nleft) - t(nleft + 1).
 
      tj0 = tj
      if (iter .eq. 1) then
c         Since we already have the smallest breakpoint we need not do
c         heapsort yet. Often only one breakpoint is used and the
c         cost of heapsort is avoided.
         tj = bkmin
         ibp = iorder(ibkmin)
      else
         if (iter .eq. 2) then
c             Replace the already used smallest breakpoint with the
c             breakpoint numbered nbreak > nlast, before heapsort call.
            if (ibkmin .ne. nbreak) then
               t(ibkmin) = t(nbreak)
               iorder(ibkmin) = iorder(nbreak)
            endif 
c        Update heap structure of breakpoints
c           (if iter=2, initialize heap).
         endif
         call hpsolb(nleft,t,iorder,iter-2)
         tj = t(nleft)
         ibp = iorder(nleft)  
      endif 
         
      dt = tj - tj0
 
      if (dt .ne. zero .and. iprint .ge. 100) then
         write (6,4011) nint,f1,f2
         write (6,5010) dt
         write (6,6010) dtm
      endif          
 
c     If a minimizer is within this interval, locate the GCP and return. 
 
      if (dtm .lt. dt) goto 888
 
c     Otherwise fix one variable and
c       reset the corresponding component of d to zero.
    
      tsum = tsum + dt
      nleft = nleft - 1
      iter = iter + 1
      dibp = d(ibp)
      d(ibp) = zero
      if (dibp .gt. zero) then
         zibp = u(ibp) - x(ibp)
         xcp(ibp) = u(ibp)
         iwhere(ibp) = 2
      else
         zibp = l(ibp) - x(ibp)
         xcp(ibp) = l(ibp)
         iwhere(ibp) = 1
      endif
      if (iprint .ge. 100) write (6,*) 'Variable  ',ibp,'  is fixed.'
      if (nleft .eq. 0 .and. nbreak .eq. n) then
c                                             all n variables are fixed,
c                                                return with xcp as GCP.
         dtm = dt
         goto 999
      endif
 
c     Update the derivative information.
 
      nint = nint + 1
      dibp2 = dibp**2
 
c     Update f1 and f2.
 
c        temporarily set f1 and f2 for col=0.
      f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
      f2 = f2 - theta*dibp2

      if (col .gt. 0) then
c                          update c = c + dt*p.
         call daxpy(col2,dt,p,1,c,1)
 
c           choose wbp,
c           the row of W corresponding to the breakpoint encountered.
         pointr = head
         do 70 j = 1,col
            wbp(j) = wy(ibp,pointr)
            wbp(col + j) = theta*ws(ibp,pointr)
            pointr = mod(pointr,m) + 1
  70     continue 
 
c           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
         call bmv(m,sy,wt,col,wbp,v,info)
         if (info .ne. 0) return
         wmc = ddot(col2,c,1,v,1)
         wmp = ddot(col2,p,1,v,1) 
         wmw = ddot(col2,wbp,1,v,1)
 
c           update p = p - dibp*wbp. 
         call daxpy(col2,-dibp,wbp,1,p,1)
 
c           complete updating f1 and f2 while col > 0.
         f1 = f1 + dibp*wmc
         f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
      endif

      f2 = max(epsmch*f2_org,f2)


      if (nleft .gt. 0) then
         iRCerr=0
         IF (DABS(f2).LE.1.D-16) THEN
            info=1
            iRCerr=1
            RETURN
         ENDIF
         dtm = -f1/f2
         goto 777
c                 to repeat the loop for unsearched intervals. 
      else if(bnded) then
         f1 = zero
         f2 = zero
         dtm = zero
      else
         iRCerr=0
         IF (DABS(f2).LE.1.D-16) THEN
            info=1
            iRCerr=1
            RETURN
         ENDIF
         dtm = -f1/f2
      endif 

c------------------- the end of the loop -------------------------------
 
 888  continue
      if (iprint .ge. 99) then
         write (6,*)
         write (6,*) 'GCP found in this segment'
         write (6,4010) nint,f1,f2
         write (6,6010) dtm
      endif 
      if (dtm .le. zero) dtm = zero
      tsum = tsum + dtm
 
c     Move free variables (i.e., the ones w/o breakpoints) and 
c       the variables whose breakpoints haven't been reached.
 
      call daxpy(n,tsum,d,1,xcp,1)
 
 999  continue
 
c     Update c = c + dtm*p = W'(x^c - x) 
c       which will be used in computing r = Z'(B(x^c - x) + g).
 
      if (col .gt. 0) call daxpy(col2,dtm,p,1,c,1)
      if (iprint .gt. 100) write (6,1010) (xcp(i),i = 1,n)
      if (iprint .ge. 99) write (6,2010)

 1010 format ('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))
 2010 format (/,'---------------- exit CAUCHY----------------------',/)
 3010 format (/,'---------------- CAUCHY entered-------------------')
 4010 format ('Piece    ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))
 4011 format (/,'Piece    ',i3,' --f1, f2 at start point ',
     +        1p,2(1x,d11.4))
 5010 format ('Distance to the next break point =  ',1p,d11.4)
 6010 format ('Distance to the stationary point =  ',1p,d11.4) 
 
      return
 
      end

c====================== The end of cauchy ==============================

      subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, 
     +                 theta, col, head, nfree, cnstnd, info)
 
      logical          cnstnd
      integer          n, m, col, head, nfree, info, index(n)
      double precision theta, 
     +                 x(n), g(n), z(n), r(n), wa(4*m), 
     +                 ws(n, m), wy(n, m), sy(m, m), wt(m, m)

c     ************
c
c     Subroutine cmprlb 
c
c       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using 
c         wa(2m+1)=W'(xcp-x) from subroutine cauchy.
c
c     Subprograms called:
c
c       L-BFGS-B Library ... bmv.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer          i,j,k,pointr
      double precision a1,a2

      if (.not. cnstnd .and. col .gt. 0) then 
         do 26 i = 1, n
            r(i) = -g(i)
  26     continue
      else
         do 30 i = 1, nfree
            k = index(i)
            r(i) = -theta*(z(k) - x(k)) - g(k)
  30     continue
         call bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
         if (info .ne. 0) then
            info = -8
            return
         endif
         pointr = head 
         do 34 j = 1, col
            a1 = wa(j)
            a2 = theta*wa(col + j)
            do 32 i = 1, nfree
               k = index(i)
               r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
  32        continue
            pointr = mod(pointr,m) + 1
  34     continue
      endif

      return

      end

c======================= The end of cmprlb =============================


      subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
 
      character*60     task
      integer          n, m, info, k, nbd(n)
      double precision factr, l(n), u(n)

c     ************
c
c     Subroutine errclb
c
c     This subroutine checks the validity of the input data.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          i
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

c     Check the input arguments for errors.

      if (n .le. 0) task = 'ERROR: N .LE. 0'
      if (m .le. 0) task = 'ERROR: M .LE. 0'
      if (factr .lt. zero) task = 'ERROR: FACTR .LT. 0'

c     Check the validity of the arrays nbd(i), u(i), and l(i).

      do 10 i = 1, n
         if (nbd(i) .lt. 0 .or. nbd(i) .gt. 3) then
c                                                   return
            task = 'ERROR: INVALID NBD'
            info = -6
            k = i
         endif
         if (nbd(i) .eq. 2) then
            if (l(i) .gt. u(i)) then
c                                    return
               task = 'ERROR: NO FEASIBLE SOLUTION'
               info = -7
               k = i
            endif
         endif
  10  continue

      return

      end

c======================= The end of errclb =============================


      subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat, 
     +                 updatd, wn, wn1, m, ws, wy, sy, theta, col,
     +                 head, info)

      integer          n, nsub, m, col, head, nenter, ileave, iupdat,
     +                 info, ind(n), indx2(n)
      double precision theta, wn(2*m, 2*m), wn1(2*m, 2*m),
     +                 ws(n, m), wy(n, m), sy(m, m)
      logical          updatd

c     ************
c
c     Subroutine formk 
c
c     This subroutine forms  the LEL^T factorization of the indefinite
c
c       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c                                                    where E = [-I  0]
c                                                              [ 0  I]
c     The matrix K can be shown to be equal to the matrix M^[-1]N
c       occurring in section 5.1 of [1], as well as to the matrix
c       Mbar^[-1] Nbar in section 5.3.
c
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     nsub is an integer variable
c       On entry nsub is the number of subspace variables in free set.
c       On exit nsub is not changed.
c
c     ind is an integer array of dimension nsub.
c       On entry ind specifies the indices of subspace variables.
c       On exit ind is unchanged. 
c
c     nenter is an integer variable.
c       On entry nenter is the number of variables entering the 
c         free set.
c       On exit nenter is unchanged. 
c
c     ileave is an integer variable.
c       On entry indx2(ileave),...,indx2(n) are the variables leaving
c         the free set.
c       On exit ileave is unchanged. 
c
c     indx2 is an integer array of dimension n.
c       On entry indx2(1),...,indx2(nenter) are the variables entering
c         the free set, while indx2(ileave),...,indx2(n) are the
c         variables leaving the free set.
c       On exit indx2 is unchanged. 
c
c     iupdat is an integer variable.
c       On entry iupdat is the total number of BFGS updates made so far.
c       On exit iupdat is unchanged. 
c
c     updatd is a logical variable.
c       On entry 'updatd' is true if the L-BFGS matrix is updatd.
c       On exit 'updatd' is unchanged. 
c
c     wn is a double precision array of dimension 2m x 2m.
c       On entry wn is unspecified.
c       On exit the upper triangle of wn stores the LEL^T factorization
c         of the 2*col x 2*col indefinite matrix
c                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c
c     wn1 is a double precision array of dimension 2m x 2m.
c       On entry wn1 stores the lower triangular part of 
c                     [Y' ZZ'Y   L_a'+R_z']
c                     [L_a+R_z   S'AA'S   ]
c         in the previous iteration.
c       On exit wn1 stores the corresponding updated matrices.
c       The purpose of wn1 is just to store these inner products
c       so they can be easily updated and inserted into wn.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     ws, wy, sy, and wtyy are double precision arrays;
c     theta is a double precision variable;
c     col is an integer variable;
c     head is an integer variable.
c       On entry they store the information defining the
c                                          limited memory BFGS matrix:
c         ws(n,m) stores S, a set of s-vectors;
c         wy(n,m) stores Y, a set of y-vectors;
c         sy(m,m) stores S'Y;
c         wtyy(m,m) stores the Cholesky factorization
c                                   of (theta*S'S+LD^(-1)L')
c         theta is the scaling factor specifying B_0 = theta I;
c         col is the number of variable metric corrections stored;
c         head is the location of the 1st s- (or y-) vector in S (or Y).
c       On exit they are unchanged.
c
c     info is an integer variable.
c       On entry info is unspecified.
c       On exit info =  0 for normal return;
c                    = -1 when the 1st Cholesky factorization failed;
c                    = -2 when the 2st Cholesky factorization failed.
c
c     Subprograms called:
c
c       Linpack ... dcopy, dpofa, dtrsl.
c
c
c     References:
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
c       limited memory FORTRAN code for solving bound constrained
c       optimization problems'', Tech. Report, NAM-11, EECS Department,
c       Northwestern University, 1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k,
     +                 col2,pbegin,pend,dbegin,dend,upcl
      double precision ddot,temp1,temp2,temp3,temp4
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

c     Form the lower triangular part of
c               WN1 = [Y' ZZ'Y   L_a'+R_z'] 
c                     [L_a+R_z   S'AA'S   ]
c        where L_a is the strictly lower triangular part of S'AA'Y
c              R_z is the upper triangular part of S'ZZ'Y.
      
      if (updatd) then
         if (iupdat .gt. m) then 
c                                 shift old part of WN1.
            do 10 jy = 1, m - 1
               js = m + jy
               call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
               call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
               call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
  10        continue
         endif
 
c          put new rows in blocks (1,1), (2,1) and (2,2).
         pbegin = 1
         pend = nsub
         dbegin = nsub + 1
         dend = n
         iy = col
         is = m + col
         ipntr = head + col - 1
         if (ipntr .gt. m) ipntr = ipntr - m    
         jpntr = head
         do 20 jy = 1, col
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
c             compute element jy of row 'col' of Y'ZZ'Y
            do 15 k = pbegin, pend
               k1 = ind(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
  15        continue
c             compute elements jy of row 'col' of L_a and S'AA'S
            do 16 k = dbegin, dend
               k1 = ind(k)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  16        continue
            wn1(iy,jy) = temp1
            wn1(is,js) = temp2
            wn1(is,jy) = temp3
            jpntr = mod(jpntr,m) + 1
  20     continue
 
c          put new column in block (2,1).
         jy = col       
         jpntr = head + col - 1
         if (jpntr .gt. m) jpntr = jpntr - m
         ipntr = head
         do 30 i = 1, col
            is = m + i
            temp3 = zero
c             compute element i of column 'col' of R_z
            do 25 k = pbegin, pend
               k1 = ind(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  25        continue 
            ipntr = mod(ipntr,m) + 1
            wn1(is,jy) = temp3
  30     continue
         upcl = col - 1
      else
         upcl = col
      endif
 
c       modify the old parts in blocks (1,1) and (2,2) due to changes
c       in the set of free variables.
      ipntr = head      
      do 45 iy = 1, upcl
         is = m + iy
         jpntr = head
         do 40 jy = 1, iy
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
            temp4 = zero
            do 35 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
  35        continue
            do 36 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
               temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
  36        continue
            wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3 
            wn1(is,js) = wn1(is,js) - temp2 + temp4 
            jpntr = mod(jpntr,m) + 1
  40     continue
         ipntr = mod(ipntr,m) + 1
  45  continue
 
c       modify the old parts in block (2,1).
      ipntr = head      
      do 60 is = m + 1, m + upcl
         jpntr = head 
         do 55 jy = 1, upcl
            temp1 = zero
            temp3 = zero
            do 50 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
  50        continue
            do 51 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  51        continue
         if (is .le. jy + m) then
               wn1(is,jy) = wn1(is,jy) + temp1 - temp3  
            else
               wn1(is,jy) = wn1(is,jy) - temp1 + temp3  
            endif
            jpntr = mod(jpntr,m) + 1
  55     continue
         ipntr = mod(ipntr,m) + 1
  60  continue
 
c     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ] 
c                                     [-L_a +R_z        S'AA'S*theta]

      m2 = 2*m
      do 70 iy = 1, col
         is = col + iy
         is1 = m + iy
         do 65 jy = 1, iy
            js = col + jy
            js1 = m + jy
            wn(jy,iy) = wn1(iy,jy)/theta
            wn(js,is) = wn1(is1,js1)*theta
  65     continue
         do 66 jy = 1, iy - 1
            wn(jy,is) = -wn1(is1,jy)
  66     continue
         do 67 jy = iy, col
            wn(jy,is) = wn1(is1,jy)
  67     continue
         wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
  70  continue

c     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')] 
c                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ]

c        first Cholesky factor (1,1) block of wn to get LL'
c                          with L' stored in the upper triangle of wn.
      call dpofa(wn,m2,col,info)
      if (info .ne. 0) then
         info = -1
         return
      endif
c        then form L^-1(-L_a'+R_z') in the (1,2) block.
      col2 = 2*col
      do 71 js = col+1 ,col2
         call dtrsl(wn,m2,col,wn(1,js),11,info)
  71  continue

c     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
c        upper triangle of (2,2) block of wn.
                      

      do 72 is = col+1, col2
         do 74 js = is, col2
               wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
  74        continue
  72     continue

c     Cholesky factorization of (2,2) block of wn.

      call dpofa(wn(col+1,col+1),m2,col,info)
      if (info .ne. 0) then
         info = -2
         return
      endif

      return

      end

c======================= The end of formk ==============================

      subroutine formt(m, wt, sy, ss, col, theta, info)
 
      integer          m, col, info
      double precision theta, wt(m, m), sy(m, m), ss(m, m)

c     ************
c
c     Subroutine formt
c
c       This subroutine forms the upper half of the pos. def. and symm.
c         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
c         of the array wt, and performs the Cholesky factorization of T
c         to produce J*J', with J' stored in the upper triangle of wt.
c
c     Subprograms called:
c
c       Linpack ... dpofa.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          i,j,k,k1
      double precision ddum
      double precision zero
      parameter        (zero=0.0d0)


c     Form the upper half of  T = theta*SS + L*D^(-1)*L',
c        store T in the upper triangle of the array wt.
 
      do 52 j = 1, col
         wt(1,j) = theta*ss(1,j)
  52  continue
      do 55 i = 2, col
         do 54 j = i, col
            k1 = min(i,j) - 1
            ddum  = zero
            do 53 k = 1, k1
               ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
  53        continue
            wt(i,j) = ddum + theta*ss(i,j)
  54     continue
  55  continue
 
c     Cholesky factorize T to J*J' with 
c        J' stored in the upper triangle of wt.
 
      call dpofa(wt,m,col,info)
      if (info .ne. 0) then
         info = -3
      endif

      return

      end

c======================= The end of formt ==============================
 
      subroutine freev(n, nfree, index, nenter, ileave, indx2, 
     +                 iwhere, wrk, updatd, cnstnd, iprint, iter)

      integer n, nfree, nenter, ileave, iprint, iter, 
     +        index(n), indx2(n), iwhere(n)
      logical wrk, updatd, cnstnd

c     ************
c
c     Subroutine freev 
c
c     This subroutine counts the entering and leaving variables when
c       iter > 0, and finds the index set of free and active variables
c       at the GCP.
c
c     cnstnd is a logical variable indicating whether bounds are present
c
c     index is an integer array of dimension n
c       for i=1,...,nfree, index(i) are the indices of free variables
c       for i=nfree+1,...,n, index(i) are the indices of bound variables
c       On entry after the first iteration, index gives 
c         the free variables at the previous iteration.
c       On exit it gives the free variables based on the determination
c         in cauchy using the array iwhere.
c
c     indx2 is an integer array of dimension n
c       On entry indx2 is unspecified.
c       On exit with iter>0, indx2 indicates which variables
c          have changed status since the previous iteration.
c       For i= 1,...,nenter, indx2(i) have changed from bound to free.
c       For i= ileave+1,...,n, indx2(i) have changed from free to bound.
c 
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer iact,i,k

      nenter = 0
      ileave = n + 1
      if (iter .gt. 0 .and. cnstnd) then
c                           count the entering and leaving variables.
         do 20 i = 1, nfree
            k = index(i)
            if (iwhere(k) .gt. 0) then
               ileave = ileave - 1
               indx2(ileave) = k
               if (iprint .ge. 100) write (6,*)
     +             'Variable ',k,' leaves the set of free variables'
            endif
  20     continue
         do 22 i = 1 + nfree, n
            k = index(i)
            if (iwhere(k) .le. 0) then
               nenter = nenter + 1
               indx2(nenter) = k
               if (iprint .ge. 100) write (6,*)
     +             'Variable ',k,' enters the set of free variables'
            endif
  22     continue
         if (iprint .ge. 99) write (6,*)
     +       n+1-ileave,' variables leave; ',nenter,' variables enter'
      endif
      wrk = (ileave .lt. n+1) .or. (nenter .gt. 0) .or. updatd
 
c     Find the index set of free and active variables at the GCP.
 
      nfree = 0 
      iact = n + 1
      do 24 i = 1, n
         if (iwhere(i) .le. 0) then
            nfree = nfree + 1
            index(nfree) = i
         else
            iact = iact - 1
            index(iact) = i
         endif
  24  continue
      if (iprint .ge. 99) write (6,*)
     +      nfree,' variables are free at GCP ',iter + 1  

      return

      end

c======================= The end of freev ==============================


      subroutine hpsolb(n, t, iorder, iheap)
      integer          iheap, n, iorder(n)
      double precision t(n)
  
c     ************
c
c     Subroutine hpsolb 
c
c     This subroutine sorts out the least element of t, and puts the
c       remaining elements of t in a heap.
c 
c     n is an integer variable.
c       On entry n is the dimension of the arrays t and iorder.
c       On exit n is unchanged.
c
c     t is a double precision array of dimension n.
c       On entry t stores the elements to be sorted,
c       On exit t(n) stores the least elements of t, and t(1) to t(n-1)
c         stores the remaining elements in the form of a heap.
c
c     iorder is an integer array of dimension n.
c       On entry iorder(i) is the index of t(i).
c       On exit iorder(i) is still the index of t(i), but iorder may be
c         permuted in accordance with t.
c
c     iheap is an integer variable specifying the task.
c       On entry iheap should be set as follows:
c         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
c         iheap .ne. 0 if otherwise.
c       On exit iheap is unchanged.
c
c
c     References:
c       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c     ************
  
      integer          i,j,k,indxin,indxou
      double precision ddum,out

      if (iheap .eq. 0) then

c        Rearrange the elements t(1) to t(n) to form a heap.

         do 20 k = 2, n
            ddum  = t(k)
            indxin = iorder(k)

c           Add ddum to the heap.
            i = k
   10       continue
            if (i.gt.1) then
               j = i/2
               if (ddum .lt. t(j)) then
                  t(i) = t(j)
                  iorder(i) = iorder(j)
                  i = j
                  goto 10 
               endif  
            endif  
            t(i) = ddum
            iorder(i) = indxin
   20    continue
      endif
 
c     Assign to 'out' the value of t(1), the least member of the heap,
c        and rearrange the remaining members to form a heap as
c        elements 1 to n-1 of t.
 
      if (n .gt. 1) then
         i = 1
         out = t(1)
         indxou = iorder(1)
         ddum  = t(n)
         indxin  = iorder(n)

c        Restore the heap 
   30    continue
         j = i+i
         if (j .le. n-1) then
            if (t(j+1) .lt. t(j)) j = j+1
            if (t(j) .lt. ddum ) then
               t(i) = t(j)
               iorder(i) = iorder(j)
               i = j
               goto 30
            endif 
         endif 
         t(i) = ddum
         iorder(i) = indxin
 
c     Put the least member in t(n). 

         t(n) = out
         iorder(n) = indxou
      endif 

      return

      end

c======================= The end of hpsolb ==============================


      subroutine lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t,
     +                  z, stp, dnorm, dtd, xstep, stpmx, iter, ifun,
     +                  iback, nfgv, info, task, boxed, cnstnd, csave,
     +                  isave, dsave)

      character*60     task, csave
      logical          boxed, cnstnd
      integer          n, iter, ifun, iback, nfgv, info,
     +                 nbd(n), isave(2)
      double precision f, fold, gd, gdold, stp, dnorm, dtd, xstep,
     +                 stpmx, x(n), l(n), u(n), g(n), d(n), r(n), t(n),
     +                 z(n), dsave(13)
c     **********
c
c     Subroutine lnsrlb
c
c     This subroutine calls subroutine dcsrch from the Minpack2 library
c       to perform the line search.  Subroutine dscrch is safeguarded so
c       that all trial points lie within the feasible region.
c
c     Subprograms called:
c
c       Minpack2 Library ... dcsrch.
c
c       Linpack ... dtrsl, ddot.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     **********

      integer          i
      double           precision ddot,a1,a2
      double precision one,zero,big
      parameter        (one=1.0d0,zero=0.0d0,big=1.0d+10)
      double precision ftol,gtol,xtol
      parameter        (ftol=1.0d-3,gtol=0.9d0,xtol=0.1d0)

      if (task(1:5) .eq. 'FG_LN') goto 556

      dtd = ddot(n,d,1,d,1)
      dnorm = sqrt(dtd)

c     Determine the maximum step length.

      stpmx = big
      if (cnstnd) then
         if (iter .eq. 0) then
            stpmx = one
         else
            do 43 i = 1, n
               a1 = d(i)
               if (nbd(i) .ne. 0) then
                  if (a1 .lt. zero .and. nbd(i) .le. 2) then
                     a2 = l(i) - x(i)
                     if (a2 .ge. zero) then
                        stpmx = zero
                     else if (a1*stpmx .lt. a2) then
                        stpmx = a2/a1
                     endif
                  else if (a1 .gt. zero .and. nbd(i) .ge. 2) then
                     a2 = u(i) - x(i)
                     if (a2 .le. zero) then
                        stpmx = zero
                     else if (a1*stpmx .gt. a2) then
                        stpmx = a2/a1
                     endif
                  endif
               endif
  43        continue
         endif
      endif
 
      if (iter .eq. 0 .and. .not. boxed) then
         stp = min(one/dnorm, stpmx)
      else
         stp = one
      endif 

      call dcopy(n,x,1,t,1)
      call dcopy(n,g,1,r,1)
      fold = f
      ifun = 0
      iback = 0
      csave = 'START'
 556  continue
      gd = ddot(n,g,1,d,1)
      if (ifun .eq. 0) then
         gdold=gd
         if (gd .ge. zero) then
c                               the directional derivative >=0.
c                               Line search is impossible.
            info = -4
            return
         endif
      endif

      call dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave)

      xstep = stp*dnorm
      if (csave(1:4) .ne. 'CONV' .and. csave(1:4) .ne. 'WARN') then
         task = 'FG_LNSRCH'
         ifun = ifun + 1
         nfgv = nfgv + 1
         iback = ifun - 1 
         if (stp .eq. one) then
            call dcopy(n,z,1,x,1)
         else
            do 41 i = 1, n
               x(i) = stp*d(i) + t(i)
  41        continue
         endif
      else
         task = 'NEW_X'
      endif

      return

      end


c======================= The end of lnsrlb ==============================


      subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail, 
     +                  iupdat, col, head, theta, rr, dr, stp, dtd)
 
      integer          n, m, itail, iupdat, col, head
      double precision theta, rr, dr, stp, dtd, d(n), r(n), 
     +                 ws(n, m), wy(n, m), sy(m, m), ss(m, m)

c     ************
c
c     Subroutine matupd
c
c       This subroutine updates matrices WS and WY, and forms the
c         middle matrix in B.
c
c     Subprograms called:
c
c       Linpack ... dcopy, ddot.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************
 
      integer          j,pointr
      double precision ddot
      double precision one
      parameter        (one=1.0d0)

c     Set pointers for matrices WS and WY.
 
      if (iupdat .le. m) then
         col = iupdat
         itail = mod(head+iupdat-2,m) + 1
      else
         itail = mod(itail,m) + 1
         head = mod(head,m) + 1
      endif
 
c     Update matrices WS and WY.

      call dcopy(n,d,1,ws(1,itail),1)
      call dcopy(n,r,1,wy(1,itail),1)
 
c     Set theta=yy/ys.
 
      theta = rr/dr
 
c     Form the middle matrix in B.
 
c        update the upper triangle of SS,
c                                         and the lower triangle of SY:
      if (iupdat .gt. m) then
c                              move old information
         do 50 j = 1, col - 1
            call dcopy(j,ss(2,j+1),1,ss(1,j),1)
            call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
  50     continue
      endif
c        add new information: the last row of SY
c                                             and the last column of SS:
      pointr = head
      do 51 j = 1, col - 1
         sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
         ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
         pointr = mod(pointr,m) + 1
  51  continue
      if (stp .eq. one) then
         ss(col,col) = dtd
      else
         ss(col,col) = stp*stp*dtd
      endif
      sy(col,col) = dr
 
      return

      end

c======================= The end of matupd =============================

      subroutine prn1lb(n, m, l, u, x, iprint, itfile, epsmch)
 
      integer n, m, iprint, itfile
      double precision epsmch, x(n), l(n), u(n)

c     ************
c
c     Subroutine prn1lb
c
c     This subroutine prints the input data, initial point, upper and
c       lower bounds of each variable, machine precision, as well as 
c       the headings of the output.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i

      if (iprint .ge. 0) then
         write (6,7001) epsmch
         write (6,*) 'N = ',n,'    M = ',m
         if (iprint .ge. 1) then
            write (itfile,2001) epsmch
            write (itfile,*)'N = ',n,'    M = ',m
            write (itfile,9001)
            if (iprint .gt. 100) then
               write (6,1004) 'L =',(l(i),i = 1,n)
               write (6,1004) 'X0 =',(x(i),i = 1,n)
               write (6,1004) 'U =',(u(i),i = 1,n)
            endif 
         endif
      endif 

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format ('RUNNING THE L-BFGS-B CODE',/,/,
     + 'it    = iteration number',/,
     + 'nf    = number of function evaluations',/,
     + 'nint  = number of segments explored during the Cauchy search',/,
     + 'nact  = number of active bounds at the generalized Cauchy point'
     + ,/,
     + 'sub   = manner in which the subspace minimization terminated:'
     + ,/,'        con = converged, bnd = a bound was reached',/,
     + 'itls  = number of iterations performed in the line search',/,
     + 'stepl = step length used',/,
     + 'tstep = norm of the displacement (total step)',/,
     + 'projg = norm of the projected gradient',/,
     + 'f     = function value',/,/,
     + '           * * *',/,/,
     + 'Machine precision =',1p,d10.3)
 7001 format ('RUNNING THE L-BFGS-B CODE',/,/,
     + '           * * *',/,/,
     + 'Machine precision =',1p,d10.3)
 9001 format (/,3x,'it',3x,'nf',2x,'nint',2x,'nact',2x,'sub',2x,'itls',
     +        2x,'stepl',4x,'tstep',5x,'projg',8x,'f')

      return

      end

c======================= The end of prn1lb =============================

      subroutine prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, 
     +                  sbgnrm, nint, word, iword, iback, stp, xstep)
 
      character*3      word
      integer          n, iprint, itfile, iter, nfgv, nact, nint,
     +                 iword, iback
      double precision f, sbgnrm, stp, xstep, x(n), g(n)

c     ************
c
c     Subroutine prn2lb
c
c     This subroutine prints out new information after a successful
c       line search. 
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i,imod

c           'word' records the status of subspace solutions.
      if (iword .eq. 0) then
c                            the subspace minimization converged.
         word = 'con'
      else if (iword .eq. 1) then
c                          the subspace minimization stopped at a bound.
         word = 'bnd'
      else if (iword .eq. 5) then
c                             the truncated Newton step has been used.
         word = 'TNT'
      else
         word = '---'
      endif
      if (iprint .ge. 99) then
         write (6,*) 'LINE SEARCH',iback,' times; norm of step = ',xstep
         write (6,2001) iter,f,sbgnrm
         if (iprint .gt. 100) then      
            write (6,1004) 'X =',(x(i), i = 1, n)
            write (6,1004) 'G =',(g(i), i = 1, n)
         endif
      else if (iprint .gt. 0) then 
         imod = mod(iter,iprint)
         if (imod .eq. 0) write (6,2001) iter,f,sbgnrm
      endif
      if (iprint .ge. 1) write (itfile,3001)
     +          iter,nfgv,nint,nact,word,iback,stp,xstep,sbgnrm,f

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format
     +  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 3001 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1p,2(1x,d10.3))

      return

      end

c======================= The end of prn2lb =============================

      subroutine prn3lb(n, x, f, task, iprint, info, itfile, 
     +                  iter, nfgv, nintol, nskip, nact, sbgnrm, 
     +                  time, nint, word, iback, stp, xstep, k, 
     +                  cachyt, sbtime, lnscht)
 
      character*60     task
      character*3      word
      integer          n, iprint, info, itfile, iter, nfgv, nintol,
     +                 nskip, nact, nint, iback, k
      double precision f, sbgnrm, time, stp, xstep, cachyt, sbtime,
     +                 lnscht, x(n)

c     ************
c
c     Subroutine prn3lb
c
c     This subroutine prints out information when either a built-in
c       convergence test is satisfied or when an error message is
c       generated.
c       
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i

      if (task(1:5) .eq. 'ERROR') goto 999

      if (iprint .ge. 0) then
         write (6,3003)
         write (6,3004)
         write(6,3005) n,iter,nfgv,nintol,nskip,nact,sbgnrm,f
         if (iprint .ge. 100) then
            write (6,1004) 'X =',(x(i),i = 1,n)
         endif  
         if (iprint .ge. 1) write (6,*) ' F =',f
      endif 
 999  continue
      if (iprint .ge. 0) then
         write (6,3009) task
         if (info .ne. 0) then
            if (info .eq. -1) write (6,9011)
            if (info .eq. -2) write (6,9012)
            if (info .eq. -3) write (6,9013)
            if (info .eq. -4) write (6,9014)
            if (info .eq. -5) write (6,9015)
            if (info .eq. -6) write (6,*)' Input nbd(',k,') is invalid.'
            if (info .eq. -7) 
     +         write (6,*)' l(',k,') > u(',k,').  No feasible solution.'
            if (info .eq. -8) write (6,9018)
            if (info .eq. -9) write (6,9019)
         endif
         if (iprint .ge. 1) write (6,3007) cachyt,sbtime,lnscht
         write (6,3008) time
         if (iprint .ge. 1) then
            if (info .eq. -4 .or. info .eq. -9) then
               write (itfile,3002)
     +             iter,nfgv,nint,nact,word,iback,stp,xstep
            endif
            write (itfile,3009) task
            if (info .ne. 0) then
               if (info .eq. -1) write (itfile,9011)
               if (info .eq. -2) write (itfile,9012)
               if (info .eq. -3) write (itfile,9013)
               if (info .eq. -4) write (itfile,9014)
               if (info .eq. -5) write (itfile,9015)
               if (info .eq. -8) write (itfile,9018)
               if (info .eq. -9) write (itfile,9019)
            endif
            write (itfile,3008) time
         endif
      endif

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 3002 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6x,'-',10x,'-')
 3003 format (/,
     + '           * * *',/,/,
     + 'Tit   = total number of iterations',/,
     + 'Tnf   = total number of function evaluations',/,
     + 'Tnint = total number of segments explored during',
     +           ' Cauchy searches',/,
     + 'Skip  = number of BFGS updates skipped',/,
     + 'Nact  = number of active bounds at final generalized',
     +          ' Cauchy point',/,
     + 'Projg = norm of the final projected gradient',/,
     + 'F     = final function value',/,/,
     + '           * * *')
 3004 format (/,3x,'N',3x,'Tit',2x,'Tnf',2x,'Tnint',2x,
     +       'Skip',2x,'Nact',5x,'Projg',8x,'F')
 3005 format (i5,2(1x,i4),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
 3006 format (i5,2(1x,i4),2(1x,i6),(1x,i4),(1x,i5),7x,'-',10x,'-')
 3007 format (/,' Cauchy                time',1p,e10.3,' seconds.',/ 
     +        ' Subspace minimization time',1p,e10.3,' seconds.',/
     +        ' Line search           time',1p,e10.3,' seconds.')
 3008 format (/,' Total User time',1p,e10.3,' seconds.',/)
 3009 format (/,a60)
 9011 format (/,
     +' Matrix in 1st Cholesky factorization in formk is not Pos. Def.')
 9012 format (/,
     +' Matrix in 2st Cholesky factorization in formk is not Pos. Def.')
 9013 format (/,
     +' Matrix in the Cholesky factorization in formt is not Pos. Def.')
 9014 format (/,
     +' Derivative >= 0, backtracking line search impossible.',/,
     +'   Previous x, f and g restored.',/,
     +' Possible causes: 1 error in function or gradient evaluation;',/,
     +'                  2 rounding errors dominate computation.')
 9015 format (/,
     +' Warning:  more than 10 function and gradient',/,
     +'   evaluations in the last line search.  Termination',/,
     +'   may possibly be caused by a bad search direction.')
 9018 format (/,' The triangular system is singular.')
 9019 format (/,
     +' Line search cannot locate an adequate point after 20 function',/
     +,'  and gradient evaluations.  Previous x, f and g restored.',/,
     +' Possible causes: 1 error in function or gradient evaluation;',/,
     +'                  2 rounding error dominate computation.')

      return

      end

c======================= The end of prn3lb =============================

      subroutine projgr(n, l, u, nbd, x, g, sbgnrm)

      integer          n, nbd(n)
      double precision sbgnrm, x(n), l(n), u(n), g(n)

c     ************
c
c     Subroutine projgr
c
c     This subroutine computes the infinity norm of the projected
c       gradient.
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer i
      double precision gi
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

      sbgnrm = zero
      do 15 i = 1, n
        gi = g(i)
        if (nbd(i) .ne. 0) then
           if (gi .lt. zero) then
              if (nbd(i) .ge. 2) gi = max((x(i)-u(i)),gi)
           else
              if (nbd(i) .le. 2) gi = min((x(i)-l(i)),gi)
           endif
        endif
        sbgnrm = max(sbgnrm,abs(gi))
  15  continue

      return

      end

c======================= The end of projgr =============================

      subroutine subsm(n, m, nsub, ind, l, u, nbd, x, d, ws, wy, theta, 
     +                 col, head, iword, wv, wn, iprint, info)
 
      integer          n, m, nsub, col, head, iword, iprint, info, 
     +                 ind(nsub), nbd(n)
      double precision theta, 
     +                 l(n), u(n), x(n), d(n), 
     +                 ws(n, m), wy(n, m), 
     +                 wv(2*m), wn(2*m, 2*m)

c     ************
c
c     Subroutine subsm
c
c     Given xcp, l, u, r, an index set that specifies
c       the active set at xcp, and an l-BFGS matrix B 
c       (in terms of WY, WS, SY, WT, head, col, and theta), 
c       this subroutine computes an approximate solution
c       of the subspace problem
c
c       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
c
c             subject to l<=x<=u
c                       x_i=xcp_i for all i in A(xcp)
c                     
c       along the subspace unconstrained Newton direction 
c       
c          d = -(Z'BZ)^(-1) r.
c
c       The formula for the Newton direction, given the L-BFGS matrix
c       and the Sherman-Morrison formula, is
c
c          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
c 
c       where
c                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c
c     Note that this procedure for computing d differs 
c     from that described in [1]. One can show that the matrix K is
c     equal to the matrix M^[-1]N in that paper.
c
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     nsub is an integer variable.
c       On entry nsub is the number of free variables.
c       On exit nsub is unchanged.
c
c     ind is an integer array of dimension nsub.
c       On entry ind specifies the coordinate indices of free variables.
c       On exit ind is unchanged.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound of x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound of x.
c       On exit u is unchanged.
c
c     nbd is a integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds, and
c                3 if x(i) has only an upper bound.
c       On exit nbd is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x specifies the Cauchy point xcp. 
c       On exit x(i) is the minimizer of Q over the subspace of
c                                                        free variables. 
c
c     d is a double precision array of dimension n.
c       On entry d is the reduced gradient of Q at xcp.
c       On exit d is the Newton direction of Q. 
c
c     ws and wy are double precision arrays;
c     theta is a double precision variable;
c     col is an integer variable;
c     head is an integer variable.
c       On entry they store the information defining the
c                                          limited memory BFGS matrix:
c         ws(n,m) stores S, a set of s-vectors;
c         wy(n,m) stores Y, a set of y-vectors;
c         theta is the scaling factor specifying B_0 = theta I;
c         col is the number of variable metric corrections stored;
c         head is the location of the 1st s- (or y-) vector in S (or Y).
c       On exit they are unchanged.
c
c     iword is an integer variable.
c       On entry iword is unspecified.
c       On exit iword specifies the status of the subspace solution.
c         iword = 0 if the solution is in the box,
c                 1 if some bound is encountered.
c
c     wv is a double precision working array of dimension 2m.
c
c     wn is a double precision array of dimension 2m x 2m.
c       On entry the upper triangle of wn stores the LEL^T factorization
c         of the indefinite matrix
c
c              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                  [L_a -R_z           theta*S'AA'S ]
c                                                    where E = [-I  0]
c                                                              [ 0  I]
c       On exit wn is unchanged.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     info is an integer variable.
c       On entry info is unspecified.
c       On exit info = 0       for normal return,
c                    = nonzero for abnormal return 
c                                  when the matrix K is ill-conditioned.
c
c     Subprograms called:
c
c       Linpack dtrsl.
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          pointr,m2,col2,ibd,jy,js,i,j,k
      double precision alpha,dk,temp1,temp2
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

      if (nsub .le. 0) return
      if (iprint .ge. 99) write (6,1001)

c     Compute wv = W'Zd.

      pointr = head 
      do 20 i = 1, col
         temp1 = zero
         temp2 = zero
         do 10 j = 1, nsub
            k = ind(j)
            temp1 = temp1 + wy(k,pointr)*d(j)
            temp2 = temp2 + ws(k,pointr)*d(j)
  10     continue
         wv(i) = temp1
         wv(col + i) = theta*temp2
         pointr = mod(pointr,m) + 1
  20  continue
 
c     Compute wv:=K^(-1)wv.

      m2 = 2*m
      col2 = 2*col
      call dtrsl(wn,m2,col2,wv,11,info)
      if (info .ne. 0) return
      do 25 i = 1, col
         wv(i) = -wv(i)
  25     continue
      call dtrsl(wn,m2,col2,wv,01,info)
      if (info .ne. 0) return
 
c     Compute d = (1/theta)d + (1/theta**2)Z'W wv.
 
      pointr = head
      do 40 jy = 1, col
         js = col + jy
         do 30 i = 1, nsub
            k = ind(i)
            d(i) = d(i) + wy(k,pointr)*wv(jy)/theta     
     +                  + ws(k,pointr)*wv(js)
  30     continue
         pointr = mod(pointr,m) + 1
  40  continue
      do 50 i = 1, nsub
         d(i) = d(i)/theta
  50  continue
 
c     Backtrack to the feasible region.
 
      alpha = one
      temp1 = alpha     
      do 60 i = 1, nsub
         k = ind(i)
         dk = d(i)
         if (nbd(k) .ne. 0) then
            if (dk .lt. zero .and. nbd(k) .le. 2) then
               temp2 = l(k) - x(k)
               if (temp2 .ge. zero) then
                  temp1 = zero
               else if (dk*alpha .lt. temp2) then
                  temp1 = temp2/dk
               endif
            else if (dk .gt. zero .and. nbd(k) .ge. 2) then
               temp2 = u(k) - x(k)
               if (temp2 .le. zero) then
                  temp1 = zero
               else if (dk*alpha .gt. temp2) then
                  temp1 = temp2/dk
               endif
            endif
            if (temp1 .lt. alpha) then
               alpha = temp1
               ibd = i
            endif
         endif
  60  continue
 
      if (alpha .lt. one) then
         dk = d(ibd)
         k = ind(ibd)
         if (dk .gt. zero) then
            x(k) = u(k)
            d(ibd) = zero
         else if (dk .lt. zero) then
            x(k) = l(k)
            d(ibd) = zero
         endif
      endif
      do 70 i = 1, nsub
         k = ind(i)
         x(k) = x(k) + alpha*d(i)
  70  continue
 
      if (iprint .ge. 99) then
         if (alpha .lt. one) then
            write (6,1002) alpha
         else
            write (6,*) 'SM solution inside the box'
         end if 
         if (iprint .gt.100) write (6,1003) (x(i),i=1,n)
      endif
 
      if (alpha .lt. one) then
         iword = 1
      else
         iword = 0
      endif 
      if (iprint .ge. 99) write (6,1004)

 1001 format (/,'----------------SUBSM entered-----------------',/)
 1002 format ( 'ALPHA = ',f7.5,' backtrack to the BOX') 
 1003 format ('Subspace solution X =  ',/,(4x,1p,6(1x,d11.4)))
 1004 format (/,'----------------exit SUBSM --------------------',/)

      return

      end
      
c====================== The end of subsm ===============================

      subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
     +                  task,isave,dsave)
      character*(*) task
      integer isave(2)
      double precision f,g,stp,ftol,gtol,xtol,stpmin,stpmax
      double precision dsave(13)
c     **********
c
c     Subroutine dcsrch
c
c     This subroutine finds a step that satisfies a sufficient
c     decrease condition and a curvature condition.
c
c     Each call of the subroutine updates an interval with 
c     endpoints stx and sty. The interval is initially chosen 
c     so that it contains a minimizer of the modified function
c
c           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
c
c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c     interval is chosen so that it contains a minimizer of f. 
c
c     The algorithm is designed to find a step that satisfies 
c     the sufficient decrease condition 
c
c           f(stp) <= f(0) + ftol*stp*f'(0),
c
c     and the curvature condition
c
c           abs(f'(stp)) <= gtol*abs(f'(0)).
c
c     If ftol is less than gtol and if, for example, the function
c     is bounded below, then there is always a step which satisfies
c     both conditions. 
c
c     If no step can be found that satisfies both conditions, then 
c     the algorithm stops with a warning. In this case stp only 
c     satisfies the sufficient decrease condition.
c
c     A typical invocation of dcsrch has the following outline:
c
c     task = 'START'
c  10 continue
c        call dcsrch( ... )
c        if (task .eq. 'FG') then
c           Evaluate the function and the gradient at stp 
c           goto 10
c           end if
c
c     NOTE: The user must no alter work arrays between calls.
c
c     The subroutine statement is
c
c        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
c                          task,isave,dsave)
c     where
c
c       f is a double precision variable.
c         On initial entry f is the value of the function at 0.
c            On subsequent entries f is the value of the 
c            function at stp.
c         On exit f is the value of the function at stp.
c
c       g is a double precision variable.
c         On initial entry g is the derivative of the function at 0.
c            On subsequent entries g is the derivative of the 
c            function at stp.
c         On exit g is the derivative of the function at stp.
c
c       stp is a double precision variable. 
c         On entry stp is the current estimate of a satisfactory 
c            step. On initial entry, a positive initial estimate 
c            must be provided. 
c         On exit stp is the current estimate of a satisfactory step
c            if task = 'FG'. If task = 'CONV' then stp satisfies
c            the sufficient decrease and curvature condition.
c
c       ftol is a double precision variable.
c         On entry ftol specifies a nonnegative tolerance for the 
c            sufficient decrease condition.
c         On exit ftol is unchanged.
c
c       gtol is a double precision variable.
c         On entry gtol specifies a nonnegative tolerance for the 
c            curvature condition. 
c         On exit gtol is unchanged.
c
c       xtol is a double precision variable.
c         On entry xtol specifies a nonnegative relative tolerance
c            for an acceptable step. The subroutine exits with a
c            warning if the relative difference between sty and stx
c            is less than xtol.
c         On exit xtol is unchanged.
c
c       stpmin is a double precision variable.
c         On entry stpmin is a nonnegative lower bound for the step.
c         On exit stpmin is unchanged.
c
c       stpmax is a double precision variable.
c         On entry stpmax is a nonnegative upper bound for the step.
c         On exit stpmax is unchanged.
c
c       task is a character variable of length at least 60.
c         On initial entry task must be set to 'START'.
c         On exit task indicates the required action:
c
c            If task(1:2) = 'FG' then evaluate the function and 
c            derivative at stp and call dcsrch again.
c
c            If task(1:4) = 'CONV' then the search is successful.
c
c            If task(1:4) = 'WARN' then the subroutine is not able
c            to satisfy the convergence conditions. The exit value of
c            stp contains the best point found during the search.
c
c            If task(1:5) = 'ERROR' then there is an error in the
c            input arguments.
c
c         On exit with convergence, a warning or an error, the
c            variable task contains additional information.
c
c       isave is an integer work array of dimension 2.
c         
c       dsave is a double precision work array of dimension 13.
c
c     Subprograms called
c
c       MINPACK-2 ... dcstep
c
c     MINPACK-1 Project. June 1983.
c     Argonne National Laboratory. 
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. October 1993.
c     Argonne National Laboratory and University of Minnesota. 
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'. 
c
c     **********
      double precision zero,p5,p66
      parameter(zero=0.0d0,p5=0.5d0,p66=0.66d0)
      double precision xtrapl,xtrapu
      parameter(xtrapl=1.1d0,xtrapu=4.0d0)

      logical brackt
      integer stage
      double precision finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest,
     +       gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1

c     Initialization block.

      if (task(1:5) .eq. 'START') then

c        Check the input arguments for errors.

         if (stp .lt. stpmin) task = 'ERROR: STP .LT. STPMIN'
         if (stp .gt. stpmax) task = 'ERROR: STP .GT. STPMAX'
         if (g .ge. zero) task = 'ERROR: INITIAL G .GE. ZERO'
         if (ftol .lt. zero) task = 'ERROR: FTOL .LT. ZERO'
         if (gtol .lt. zero) task = 'ERROR: GTOL .LT. ZERO'
         if (xtol .lt. zero) task = 'ERROR: XTOL .LT. ZERO'
         if (stpmin .lt. zero) task = 'ERROR: STPMIN .LT. ZERO'
         if (stpmax .lt. stpmin) task = 'ERROR: STPMAX .LT. STPMIN'

c        Exit if there are errors on input.

         if (task(1:5) .eq. 'ERROR') return

c        Initialize local variables.

         brackt = .false.
         stage = 1
         finit = f
         ginit = g
         gtest = ftol*ginit
         width = stpmax - stpmin
         width1 = width/p5

c        The variables stx, fx, gx contain the values of the step, 
c        function, and derivative at the best step. 
c        The variables sty, fy, gy contain the value of the step, 
c        function, and derivative at sty.
c        The variables stp, f, g contain the values of the step, 
c        function, and derivative at stp.

         stx = zero
         fx = finit
         gx = ginit
         sty = zero
         fy = finit
         gy = ginit
         stmin = zero
         stmax = stp + xtrapu*stp
         task = 'FG'

         goto 1000

      else

c        Restore local variables.

         if (isave(1) .eq. 1) then
            brackt = .true.
         else
            brackt = .false.
         endif
         stage = isave(2) 
         ginit = dsave(1) 
         gtest = dsave(2) 
         gx = dsave(3) 
         gy = dsave(4) 
         finit = dsave(5) 
         fx = dsave(6) 
         fy = dsave(7) 
         stx = dsave(8) 
         sty = dsave(9) 
         stmin = dsave(10) 
         stmax = dsave(11) 
         width = dsave(12) 
         width1 = dsave(13) 

      endif

c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c     algorithm enters the second stage.

      ftest = finit + stp*gtest
      if (stage .eq. 1 .and. f .le. ftest .and. g .ge. zero) 
     +   stage = 2

c     Test for warnings.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     +   task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      if (brackt .and. stmax - stmin .le. xtol*stmax) 
     +   task = 'WARNING: XTOL TEST SATISFIED'
      if (stp .eq. stpmax .and. f .le. ftest .and. g .le. gtest) 
     +   task = 'WARNING: STP = STPMAX'
      if (stp .eq. stpmin .and. (f .gt. ftest .or. g .ge. gtest)) 
     +   task = 'WARNING: STP = STPMIN'

c     Test for convergence.

      if (f .le. ftest .and. abs(g) .le. gtol*(-ginit)) 
     +   task = 'CONVERGENCE'

c     Test for termination.

      if (task(1:4) .eq. 'WARN' .or. task(1:4) .eq. 'CONV') goto 1000

c     A modified function is used to predict the step during the
c     first stage if a lower function value has been obtained but 
c     the decrease is not sufficient.

      if (stage .eq. 1 .and. f .le. fx .and. f .gt. ftest) then

c        Define the modified function and derivative values.

         fm = f - stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest

c        Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,
     +               brackt,stmin,stmax)

c        Reset the function and derivative values for f.

         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      else

c       Call dcstep to update stx, sty, and to compute the new step.

        call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,
     +              brackt,stmin,stmax)

      endif

c     Decide if a bisection step is needed.

      if (brackt) then
         if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty - stx)
         width1 = width
         width = abs(sty-stx)
      endif

c     Set the minimum and maximum steps allowed for stp.

      if (brackt) then
         stmin = min(stx,sty)
         stmax = max(stx,sty)
      else
         stmin = stp + xtrapl*(stp - stx)
         stmax = stp + xtrapu*(stp - stx)
      endif
 
c     Force the step to be within the bounds stpmax and stpmin.
 
      stp = max(stp,stpmin)
      stp = min(stp,stpmax)

c     If further progress is not possible, let stp be the best
c     point obtained during the search.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax)
     +   .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx

c     Obtain another function and derivative.

      task = 'FG'

 1000 continue

c     Save local variables.

      if (brackt) then
         isave(1) = 1
      else
         isave(1) = 0
      endif
      isave(2) = stage
      dsave(1) =  ginit
      dsave(2) =  gtest
      dsave(3) =  gx
      dsave(4) =  gy
      dsave(5) =  finit
      dsave(6) =  fx
      dsave(7) =  fy
      dsave(8) =  stx
      dsave(9) =  sty
      dsave(10) = stmin
      dsave(11) = stmax
      dsave(12) = width
      dsave(13) = width1

      end
      
c====================== The end of dcsrch ==============================

      subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
     +                  stpmin,stpmax)
      logical brackt
      double precision stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
c     **********
c
c     Subroutine dcstep
c
c     This subroutine computes a safeguarded step for a search
c     procedure and updates an interval that contains a step that
c     satisfies a sufficient decrease and a curvature condition.
c
c     The parameter stx contains the step with the least function
c     value. If brackt is set to .true. then a minimizer has
c     been bracketed in an interval with endpoints stx and sty.
c     The parameter stp contains the current step. 
c     The subroutine assumes that if brackt is set to .true. then
c
c           min(stx,sty) < stp < max(stx,sty),
c
c     and that the derivative at stx is negative in the direction 
c     of the step.
c
c     The subroutine statement is
c
c       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
c                         stpmin,stpmax)
c
c     where
c
c       stx is a double precision variable.
c         On entry stx is the best step obtained so far and is an
c            endpoint of the interval that contains the minimizer. 
c         On exit stx is the updated best step.
c
c       fx is a double precision variable.
c         On entry fx is the function at stx.
c         On exit fx is the function at stx.
c
c       dx is a double precision variable.
c         On entry dx is the derivative of the function at 
c            stx. The derivative must be negative in the direction of 
c            the step, that is, dx and stp - stx must have opposite 
c            signs.
c         On exit dx is the derivative of the function at stx.
c
c       sty is a double precision variable.
c         On entry sty is the second endpoint of the interval that 
c            contains the minimizer.
c         On exit sty is the updated endpoint of the interval that 
c            contains the minimizer.
c
c       fy is a double precision variable.
c         On entry fy is the function at sty.
c         On exit fy is the function at sty.
c
c       dy is a double precision variable.
c         On entry dy is the derivative of the function at sty.
c         On exit dy is the derivative of the function at the exit sty.
c
c       stp is a double precision variable.
c         On entry stp is the current step. If brackt is set to .true.
c            then on input stp must be between stx and sty. 
c         On exit stp is a new trial step.
c
c       fp is a double precision variable.
c         On entry fp is the function at stp
c         On exit fp is unchanged.
c
c       dp is a double precision variable.
c         On entry dp is the the derivative of the function at stp.
c         On exit dp is unchanged.
c
c       brackt is an logical variable.
c         On entry brackt specifies if a minimizer has been bracketed.
c            Initially brackt must be set to .false.
c         On exit brackt specifies if a minimizer has been bracketed.
c            When a minimizer is bracketed brackt is set to .true.
c
c       stpmin is a double precision variable.
c         On entry stpmin is a lower bound for the step.
c         On exit stpmin is unchanged.
c
c       stpmax is a double precision variable.
c         On entry stpmax is an upper bound for the step.
c         On exit stpmax is unchanged.
c
c     MINPACK-1 Project. June 1983
c     Argonne National Laboratory. 
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. October 1993.
c     Argonne National Laboratory and University of Minnesota. 
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision zero,p66,two,three
      parameter(zero=0.0d0,p66=0.66d0,two=2.0d0,three=3.0d0)
   
      double precision gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta

      sgnd = dp*(dx/abs(dx))

c     First case: A higher function value. The minimum is bracketed. 
c     If the cubic step is closer to stx than the quadratic step, the 
c     cubic step is taken, otherwise the average of the cubic and 
c     quadratic steps is taken.

      if (fp .gt. fx) then
         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)*
     +                                                       (stp - stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
            stpf = stpc + (stpq - stpc)/two
         endif
         brackt = .true.

c     Second case: A lower function value and derivatives of opposite 
c     sign. The minimum is bracketed. If the cubic step is farther from 
c     stp than the secant step, the cubic step is taken, otherwise the 
c     secant step is taken.

      else if (sgnd .lt. zero) then
         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = ((gamma - dp) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp/(dp - dx))*(stx - stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
         endif
         brackt = .true.

c     Third case: A lower function value, derivatives of the same sign,
c     and the magnitude of the derivative decreases.

      else if (abs(dp) .lt. abs(dx)) then

c        The cubic step is computed only if the cubic tends to infinity 
c        in the direction of the step or if the minimum of the cubic
c        is beyond stp. Otherwise the cubic step is defined to be the 
c        secant step.

         theta = three*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))

c        The case gamma = 0 only arises if the cubic does not tend
c        to infinity in the direction of the step.

         gamma = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. zero .and. gamma .ne. zero) then
            stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
         endif
         stpq = stp + (dp/(dp - dx))*(stx - stp)

         if (brackt) then

c           A minimizer has been bracketed. If the cubic step is 
c           closer to stp than the secant step, the cubic step is 
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .lt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            endif
            if (stp .gt. stx) then
               stpf = min(stp+p66*(sty-stp),stpf)
            else
               stpf = max(stp+p66*(sty-stp),stpf)
            endif
         else

c           A minimizer has not been bracketed. If the cubic step is 
c           farther from stp than the secant step, the cubic step is 
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .gt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            endif
            stpf = min(stpmax,stpf)
            stpf = max(stpmin,stpf)
         endif

c     Fourth case: A lower function value, derivatives of the same sign, 
c     and the magnitude of the derivative does not decrease. If the 
c     minimum is not bracketed, the step is either stpmin or stpmax, 
c     otherwise the cubic step is taken.

      else
         if (brackt) then
            theta = three*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
         endif
      endif

c     Update the interval which contains a minimizer.

      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. zero) then
            sty = stx
            fy = fx
            dy = dx
         endif
         stx = stp
         fx = fp
         dx = dp
      endif

c     Compute the new step.

      stp = stpf

      end
      
c====================== The end of dcstep ==============================

      subroutine timer(ttime)
      double precision ttime
c     *********
c
c     Subroutine timer
c
c     This subroutine is used to determine user time. In a typical 
c     application, the user time for a code segment requires calls 
c     to subroutine timer to determine the initial and final time.
c
c     The subroutine statement is
c
c       subroutine timer(ttime)
c
c     where
c
c       ttime is an output variable which specifies the user time.
c
c     Argonne National Laboratory and University of Minnesota.
c     MINPACK-2 Project.
c
c     Modified October 1990 by Brett M. Averick.
c
c     **********
      real temp
      real tarray(2)
      real etime

c     The first element of the array tarray specifies user time

      temp = etime(tarray) 

      ttime = dble(tarray(1))
 
      return

      end
      
c====================== The end of timer ===============================


      double precision function dpmeps()
c     **********
c
c     Subroutine dpeps
c
c     This subroutine computes the machine precision parameter
c     dpmeps as the smallest floating point number such that
c     1 + dpmeps differs from 1.
c
c     This subroutine is based on the subroutine machar described in
c
c     W. J. Cody,
c     MACHAR: A subroutine to dynamically determine machine parameters,
c     ACM Transactions on Mathematical Software, 14, 1988, pages 303-311.
c
c     The subroutine statement is:
c
c       subroutine dpeps(dpmeps)
c
c     where
c
c       dpmeps is a double precision variable.
c         On entry dpmeps need not be specified.
c         On exit dpmeps is the machine precision.
c
c     MINPACK-2 Project. February 1991.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     *******
      integer i,ibeta,irnd,it,itemp,negep
      double precision a,b,beta,betain,betah,temp,tempa,temp1,
     +       zero,one,two
      data zero,one,two /0.0d0,1.0d0,2.0d0/
 
c     determine ibeta, beta ala malcolm.

      a = one
      b = one
   10 continue
         a = a + a
         temp = a + one
         temp1 = temp - a
      if (temp1 - one .eq. zero) go to 10
   20 continue
         b = b + b
         temp = a + b
         itemp = int(temp - a)
      if (itemp .eq. 0) go to 20
      ibeta = itemp
      beta = dble(ibeta)

c     determine it, irnd.

      it = 0
      b = one
   30 continue
         it = it + 1
         b = b * beta
         temp = b + one
         temp1 = temp - b
      if (temp1 - one .eq. zero) go to 30
      irnd = 0
      betah = beta/two
      temp = a + betah
      if (temp - a .ne. zero) irnd = 1
      tempa = a + beta
      temp = tempa + betah
      if ((irnd .eq. 0) .and. (temp - tempa .ne. zero)) irnd = 2

c     determine dpmeps.

      negep = it + 3
      betain = one/beta
      a = one
      do 40 i = 1, negep
         a = a*betain
   40 continue
   50 continue
        temp = one + a
        if (temp - one .ne. zero) go to 60
        a = a*beta
        go to  50
   60 continue
      dpmeps = a
      if ((ibeta .eq. 2) .or. (irnd .eq. 0)) go to 70
      a = (a*(one + a))/two
      temp = one + a
      if (temp - one .ne. zero) dpmeps = a

   70 return

      end

c====================== The end of dpmeps ===============================


      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      
c====================== The end of daxpy ===============================

      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      
c====================== The end of dcopy ===============================


      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      
c====================== The end of dcopy ===============================


      subroutine dpofa(a,lda,n,info)
      integer lda,n,info
      double precision a(lda,1)
c
c     dpofa factors a double precision symmetric positive definite
c     matrix.
c
c     dpofa is usually called by dpoco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the symmetric matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a = trans(r)*r
c                where  trans(r)  is the transpose.
c                the strict lower triangle is unaltered.
c                if  info .ne. 0 , the factorization is not complete.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas ddot
c     fortran sqrt
c
c     internal variables
c
      double precision ddot,t
      double precision s
      integer j,jm1,k
c     begin block with ...exits to 40
c
c
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
c     ......exit
            if (s .le. 0.0d0) go to 40
            a(j,j) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
      
c====================== The end of dpofa ===============================

      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      
c====================== The end of dscal ===============================

      subroutine dtrsl(t,ldt,n,b,job,info)
      integer ldt,n,job,info
      double precision t(ldt,1),b(1)
c
c
c     dtrsl solves systems of the form
c
c                   t * x = b
c     or
c                   trans(t) * x = b
c
c     where t is a triangular matrix of order n. here trans(t)
c     denotes the transpose of the matrix t.
c
c     on entry
c
c         t         double precision(ldt,n)
c                   t contains the matrix of the system. the zero
c                   elements of the matrix are not referenced, and
c                   the corresponding elements of the array can be
c                   used to store other information.
c
c         ldt       integer
c                   ldt is the leading dimension of the array t.
c
c         n         integer
c                   n is the order of the system.
c
c         b         double precision(n).
c                   b contains the right hand side of the system.
c
c         job       integer
c                   job specifies what kind of system is to be solved.
c                   if job is
c
c                        00   solve t*x=b, t lower triangular,
c                        01   solve t*x=b, t upper triangular,
c                        10   solve trans(t)*x=b, t lower triangular,
c                        11   solve trans(t)*x=b, t upper triangular.
c
c     on return
c
c         b         b contains the solution, if info .eq. 0.
c                   otherwise b is unaltered.
c
c         info      integer
c                   info contains zero if the system is nonsingular.
c                   otherwise info contains the index of
c                   the first zero diagonal element of t.
c
c     linpack. this version dated 08/14/78 .
c     g. w. stewart, university of maryland, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran mod
c
c     internal variables
c
      double precision ddot,temp
      integer case,j,jj
c
c     begin block permitting ...exits to 150
c
c        check for zero diagonal elements.
c
         do 10 info = 1, n
c     ......exit
            if (t(info,info) .eq. 0.0d0) go to 150
   10    continue
         info = 0
c
c        determine the task and go to it.
c
         case = 1
         if (mod(job,10) .ne. 0) case = 2
         if (mod(job,100)/10 .ne. 0) case = case + 2
         go to (20,50,80,110), case
c
c        solve t*x=b for t lower triangular
c
   20    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 40
            do 30 j = 2, n
               temp = -b(j-1)
               call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
               b(j) = b(j)/t(j,j)
   30       continue
   40       continue
         go to 140
c
c        solve t*x=b for t upper triangular.
c
   50    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 70
            do 60 jj = 2, n
               j = n - jj + 1
               temp = -b(j+1)
               call daxpy(j,temp,t(1,j+1),1,b(1),1)
               b(j) = b(j)/t(j,j)
   60       continue
   70       continue
         go to 140
c
c        solve trans(t)*x=b for t lower triangular.
c
   80    continue
            b(n) = b(n)/t(n,n)
            if (n .lt. 2) go to 100
            do 90 jj = 2, n
               j = n - jj + 1
               b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
               b(j) = b(j)/t(j,j)
   90       continue
  100       continue
         go to 140
c
c        solve trans(t)*x=b for t upper triangular.
c
  110    continue
            b(1) = b(1)/t(1,1)
            if (n .lt. 2) go to 130
            do 120 j = 2, n
               b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
               b(j) = b(j)/t(j,j)
  120       continue
  130       continue
  140    continue
  150 continue
      return
      end
      
c====================== The end of dtrsl ===============================





c
c $Id: gordon_o2.f,v 1.1 1995/02/10 23:10:56 jim Exp $
c
c $Log: gordon_o2.f,v $
c Revision 1.1  1995/02/10  23:10:56  jim
c Add O_2 correction for SeaWIFS band 7.
c Corrections/space savings to version 2 atmospheric algorithms.
c
C
C Copyright 1988-1995 by Rosenstiel School of Marine and Atmospheric Science,
C University of Miami, Miami, Florida.
C
C                       All Rights Reserved
C
C Permission to use, copy, modify, and distribute this software and its
C documentation for non-commercial purposes and without fee is hereby granted,
C provided that the above copyright notice appear in all copies and that both
C that copyright notice and this permission notice appear in supporting
C documentation, and that the names of University of Miami and/or RSMAS not be
C used in advertising or publicity pertaining to distribution of the software
C without specific, written prior permission. 
C
C UNIVERSITY OF MIAMI DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
C INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT
C SHALL UNIVERSITY OF MIAMI BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
C DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
C WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING
C OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 
C
c Subject: Oxygen correction for band 7 on SeaWiFS.
c 
c Jim:  Here is how the O2 correction works. First, recall that only 
c SeaWiFS Band 7 need be corrected, so the notes below refer to Band 7.  
c The other bands are not corrected.
c
c       (1) compute the rayleigh as usual to yield "ray_rad"
c
c       (2) compute the two-way air mass "air_mass"
c
c       air_mass=(1.0/cos(theta))+(1.0/cos(theta_0))
c
c       (3) use function "funct_oxygen_ray(air_mass)" to compute the O2 
c       correction to the rayliegh component yielding "ray_rad_with"
c
c       ray_rad_with = ray_rad * funct_oxygen_ray(air_mass)  ! with O2
c
c       (4) subtract the "ray_rad_with" from the total to get the 
c       aerosol, Note if we are working with radiance, we need to switch to 
c       reflectance to operate the algorithm. (If we are already using 
c       reflectance ignore the "pi/cos(theta_0)" below
c
c       refl(7) = pi*(tal_rad - ray_rad_with)/(cos(theta_0)) 
c
c       (5) correct aerosol reflectance in band 7 for O2
c
c       refl(7) = refl(7) * funct_oxygen_aer(air_mass)          ! with O2
c
c       (6) use "refl(7)" in the algorithm as usual. 
c
c  *********************************************************************

        function funct_oxygen_aer(air_mass)

c       converts aerosol at 765 with Oxygen absorption to 
c       aerosol at 765 without Oxygen absorption 
c
c       Base case: m80_t50_strato  visibility (550nm ,0-2km):25km 

        real funct_oxygen_aer
        real air_mass,  a(3)
        real power

        data a / -1.0796,9.0481e-2,-6.8452e-3 /

        power = a(1) + air_mass * a(2) + (air_mass ** 2) * a(3)
        funct_oxygen_aer = 1. + 10.**power              

        return 
        end

c  *********************************************************************

        function funct_oxygen_ray(air_mass)

c       converts rayleigh at 765 without Oxygen absorption to 
c       rayleigh at 765 with Oxygen absorption 
c
c       Base case here is the 1976 Standard atmosphere without aerosols
        
        real funct_oxygen_ray
        real air_mass,  a(3)
        real power, x

        data a / -1.3491, 0.1155, -7.0218e-3 /

        power = a(1) + air_mass * a(2) + (air_mass ** 2) * a(3)
        x = 1. + 10.**power             
        funct_oxygen_ray = 1. / x               

        return 
        end
