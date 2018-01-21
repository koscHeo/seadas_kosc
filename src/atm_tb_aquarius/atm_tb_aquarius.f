      PROGRAM MAIN
      USE HDF5                  ! This module contains all necessary modules 
        
      IMPLICIT NONE

      REAL(4)    :: XLAT, XLON, XTHT, XTRAN, XTBUP, XTBDW, XVAP, XCLD, XAO, XAV, XAL
      INTEGER(4) :: ILAT, ILON, ITHT
      INTEGER(4) :: IERR

      CHARACTER(LEN=100) :: atmfile
      REAL(4), TARGET    :: TBUP(0:359,-90:90,0:90), TBDW(0:359,-90:90,0:90), TRAN(0:359,-90:90,0:90)
      REAL(4), TARGET    ::   AO(0:359,-90:90),   AV(0:359,-90:90),   AL(0:359,-90:90)

      REAL(4), PARAMETER                :: FREQ=1.413
      INTEGER(4), PARAMETER     :: ICLD=1

      REAL(4), PARAMETER :: INVALID=9999.0
      CHARACTER(len=128) :: metfile

      INTEGER(HID_T) :: file_id   
      INTEGER(HID_T), DIMENSION(6) :: dset_id
      INTEGER :: rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dims
      INTEGER(HSIZE_T), DIMENSION(3) :: offset=(/0,0,0/)
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace
      integer :: nargs

      INTEGER(HID_T) :: attr_id ! Attribute identifier 
      INTEGER(HID_T) :: aspace_id ! Attribute Dataspace identifier 
      INTEGER(HID_T) :: atype_id      ! Attribute Dataspace identifier 
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string

      CHARACTER(LEN=8), DIMENSION(1) ::  attr_data  ! Attribute data

      INTEGER     ::   error ! Error flag
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
     
      real(4), POINTER :: pt
      CHARACTER(LEN=8) version

      include "aquarius.fin"

      version = '1.000'
      write(*,*) 'atm_tb_aquarius '//version//'('//__DATE__//' '//__TIME__//')'
      write(*,*) ''

      nargs = IArgC()
      IF (nargs.eq.0) THEN
         write(*,*) 'atm_tb_aquarius input_qmet_file output_qatm_file'
         write(*,*) ' '
         write(*,*) 'where (typically) input_qmet_file:  inpath/Nyyyydoyhh_QMET_NCEP_6h'
         write(*,*) '                  output_qmet_file: outpath/Nyyyydoyhh_QATM_NCEP_6h.h5'
         write(*,*) ' '
         write(*,*) '      yyyy is 4-digit year'
         write(*,*) '      doy  is 3-digit day of year'
         write(*,*) ' '
         stop
      ENDIF

      CALL getarg(1, metfile)
      CALL getarg(2, atmfile)

c     Get basename of metfile and extract hour 
c      slash_pos = INDEX(metfile,'/',BACK=.TRUE.)
c      write(*,*) slash_pos
c      write(*,*) metfile(slash_pos+9:slash_pos+10)
c      READ(UNIT=metfile(slash_pos+9:slash_pos+10), FMT='(I4)') HH
c      write(*,*) hh

c      call getarg(2, str)
c      READ(UNIT=str, FMT='(I4)') LYEAR

      call create_hdf5(atmfile, file_id)

C     Write version attribute
      attrlen = 8
      attr_data(1) = version
      write(*,*) attr_data(1)
      CALL h5screate_simple_f( arank, adims, aspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attrlen, error)

      CALL h5acreate_f(file_id, 'Version', atype_id, aspace_id, attr_id, error)
      data_dims(1) = 8
      CALL h5awrite_f(attr_id, atype_id, attr_data, data_dims, error)

      rank = 2
      dims(1) = 360
      dims(2) = 181
      dims(3) = 91
      call create_hdf5_real(file_id, 'Absorption_O2', rank, dims, dset_id(1))
      call create_hdf5_real(file_id, 'Absorption_H2O', rank, dims, dset_id(2))
      call create_hdf5_real(file_id, 'Absorption_CLD', rank, dims, dset_id(3))

      rank = 3
      call create_hdf5_real(file_id, 'Transmittance', rank, dims, dset_id(4))
      call create_hdf5_real(file_id, 'Upwelling_TB', rank, dims, dset_id(5))
      call create_hdf5_real(file_id, 'Downwelling_TB', rank, dims, dset_id(6))

      DO ITHT=0,90,1
         WRITE(*,*) ITHT
         call flush(6)

         DO ILAT=-90,90,1
            DO ILON=0,359,1     
               XLAT=ILAT
               XLON=ILON
               XTHT=ITHT 
               CALL  GET_ATM(icld,xlat,xlon,xtht,freq,
     2              xvap,xcld,xtbup,xtbdw,xtran,xao,xav,xal,metfile,ierr) 

               IF (IERR /=0) THEN
                  WRITE(*,*) ' ERROR',IERR
                  WRITE(*,*) ILAT,ILON
                  TRAN(ILON,ILAT,ITHT)=INVALID
                  TBUP(ILON,ILAT,ITHT)=INVALID
                  TBDW(ILON,ILAT,ITHT)=INVALID
                  IF (ITHT==0) AO(ILON,ILAT)=INVALID
                  IF (ITHT==0) AV(ILON,ILAT)=INVALID
                  IF (ITHT==0) AL(ILON,ILAT)=INVALID
               ELSE
                  TRAN(ILON,ILAT,ITHT)=XTRAN
                  TBUP(ILON,ILAT,ITHT)=XTBUP
                  TBDW(ILON,ILAT,ITHT)=XTBDW
                  IF (ITHT==0) AO(ILON,ILAT)=XAO ! Absorption_O2 
                  IF (ITHT==0) AV(ILON,ILAT)=XAV ! Absorption_H20
                  IF (ITHT==0) AL(ILON,ILAT)=XAL ! Absorption_CLD
               ENDIF
c     WRITE(*,*) itht,ilat,ilon,xtran,xtbup,xtbdw
            ENDDO               !LON
         ENDDO                  !LAT
      ENDDO                     !THT 

           
      rank = 2
      call set_space_hdf5(dset_id(1), rank, dims, filespace, dataspace)
      pt => AO(0,-90)
      call write_hdf5_real(dset_id(1), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(1), dataspace, filespace)

      call set_space_hdf5(dset_id(2), rank, dims, filespace, dataspace)
      pt => AV(0,-90)
      call write_hdf5_real(dset_id(2), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(2), dataspace, filespace)

      call set_space_hdf5(dset_id(3), rank, dims, filespace, dataspace)
      pt => AL(0,-90)
      call write_hdf5_real(dset_id(3), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(3), dataspace, filespace)


!               WRITE(1) AO
!               WRITE(1) AV
!               WRITE(1) AL

      rank = 3
      call set_space_hdf5(dset_id(4), rank, dims, filespace, dataspace)
      pt => TRAN(0,-90,0)
      call write_hdf5_real(dset_id(4), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(4), dataspace, filespace)

      call set_space_hdf5(dset_id(5), rank, dims, filespace, dataspace)
      pt => TBUP(0,-90,0)
      call write_hdf5_real(dset_id(5), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(5), dataspace, filespace)

      call set_space_hdf5(dset_id(6), rank, dims, filespace, dataspace)
      pt => TBDW(0,-90,0)
      call write_hdf5_real(dset_id(6), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(6), dataspace, filespace)

!            WRITE(1) TRAN
!            WRITE(1) TBUP
!            WRITE(1) TBDW

      call close_hdf5_df(file_id)

      STOP ' NORM END'
      END   

      SUBROUTINE GET_ATM(icld,xlat,xlon,tht,freq,
     2     vap,cld,tbup,tbdw,tran,AO,AV,AL,metfile,ierr)

c       USE PORTLIB 

c    adapted for aquarius
c    EIA BETWEEN 0 and 90 DEG
c    FEB 06 2006


c    NCEP profiles 
c     
c     RSS vapor absorption 
c     Rosenkranz O2 absorption
c     Meissner dielectric model for cloud water
c    
c     icld = 0: no clouds
c     icld = 1: NCEP cloud density, only liquid part 
c
c     ierr: = 0 normal return
c           =-1 no NCEP maps
c           =+1 other error

         implicit none

         INTEGER, PARAMETER :: NMAX = 26
         REAL(4) ::   xlat,xlon,tht,freq,vap,cld,tbup,tbdw,tran,cld0
         INTEGER:: ierr,ipr,ibegin,icld
         real(4)  :: AV,AO,AL
         REAL, DIMENSION(0:NMAX) :: T,P,PV,Z,RHOV,RHOL,RHOL0
         REAL(4) :: PWAT,CWAT
         REAL :: ABH2O(0:nmax),ABO2(0:nmax),ABCLD(0:nmax), 
     1        TABS(0:nmax) 
         INTEGER :: IERR1,IERR2
         CHARACTER(len=128) :: metfile
         REAL(4), PARAMETER :: TC = 2.7


c      NCEP PARAMETERS   
       CALL FINDNCEP(NMAX,xlat,xlon,VAP,CLD0,
     1        PWAT,CWAT,P,T,PV,RHOV,RHOL0,z,ibegin,metfile,ierr1)

         ierr = ierr1
         if (ierr1 .ne. 0) return

         if (icld == 0) then ! no cloud
                 RHOL = 0.0
               CLD = 0.0
         else if (icld == 1)    then ! NCEP, only liquid part
c              liquid / ice densities
               RHOL = RHOL0
                 CLD = CLD0
       endif

c      
c     Atmospheric Absorption
c     O2 from Rosenkranz 
c     H2O vapor from AMSR ATBD
c     CLOUD with Meissner dielectric constant 


         DO IPR  = IBEGIN,NMAX
         CALL FDABSCOEFF(freq,P(ipr),T(ipr),PV(ipr), 
     1           ABH2O(ipr),ABO2(ipr))
       ! [ neper/km]

         if (RHOL(IPR) > 1.0E-7) then
         CALL FDCLDABS(freq,T(ipr),RHOL(ipr),   ABCLD(ipr))
         else
         ABCLD(ipr) = 0.0
         endif


         ABH2O(ipr) = 
     1ABH2O(ipr)/1000.0  ! [neper/m]

         ABO2(ipr) = 
     1ABO2(ipr)/1000.0           ! [neper/m]

         ABCLD(ipr) = 
     1ABCLD(ipr)/1000.0          ! [neper/m]


         ENDDO ! IPR

c      vertical integrals

         CALL COLUMN (nmax-ibegin,z(ibegin:nmax),
     1ABH2O(ibegin:nmax),2,   AV,ierr2)

         CALL COLUMN (nmax-ibegin,z(ibegin:nmax),
     1ABO2(ibegin:nmax),2,    AO,ierr2)

         CALL COLUMN (nmax-ibegin,z(ibegin:nmax),
     1ABCLD(ibegin:nmax),1,   AL,ierr2)


       ierr = ierr2
         if (ierr2 .ne. 0) then
         return
         endif

c       total absorption 
         TABS(ibegin:nmax) = ABH2O(ibegin:nmax) + ABO2 (ibegin:nmax) + 
     2        ABCLD(ibegin:nmax)

         CALL ATM_TRAN(NMAX-IBEGIN,THT,T(ibegin:nmax),Z(ibegin:nmax),
     1         TABS(ibegin:nmax),
     3         TRAN,TBDW,TBUP)


         RETURN

         END


       SUBROUTINE ATM_TRAN(NLEV,THT,T,Z,TABS,  TRAN,TBDW,TBUP)
       IMPLICIT NONE

c     SUBSTITUTE SECTHT BY DSDH IN ORDER TO INCLUDE GRAZING ANGLES

C       COMPUTER ATMOSPHERIC DOWNWELLING AND UPWELLING BRIGHTNESS TEMPERATURES
C       AND UPWARD TRANSMITTANCE AT EACH PRESSURE LEVEL (ALTITUDE) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       Input:
c      NLEV           NUMBER OF ATMOSPHERE LEVELS
c      THT            Earth Incidence Angle [in deg]
c      TABS(0:NLEV)   ATMOSPHRIC ABSORPTRION COEFFICIENTS [NEPERS/M]
c      T(0:NLEV)      TEMPERATURE PROFILE[in K]
C       Z(0:NLEV)      ELEVATION (M) 


c      Output:
c      TRAN             total atmospheric transmission
c      TBDW                     downwelling brightness temperature T_BD [in K]
c      TBUP                     upwelling   brightness temperature T_BU [in K]  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


         INTEGER(4) :: NLEV,I
         REAL(4) :: T(0:NLEV),Z(0:NLEV),TABS(0:NLEV)
         REAL(4) :: OPACTY(NLEV),TAVG(NLEV),EMS(NLEV)
         REAL(4), PARAMETER :: RE=6378.135, DELTA=0.00035
         REAL(4) :: SUMOP, SUMDW, SUMUP, TRAN ,TBAVG, TBDW, TBUP, SECTHT, DSDH, THT

        real(4) cosd

       SECTHT=1./COSD(THT)
         DSDH = (1.0+DELTA)/SQRT(COSD(THT)*COSD(THT) + DELTA*(2+DELTA))


       DO I=1,NLEV
         OPACTY(I)=-DSDH*0.5*(TABS(I-1)+TABS(I))*(Z(I)-Z(I-1))
       TAVG(I)  =0.5*(T(I-1)+T(I))
       EMS(I)   =1.-EXP(OPACTY(I))
         ENDDO

       SUMOP=0 
         SUMDW=0
       DO I=1,NLEV
       SUMDW=SUMDW+(TAVG(I)-T(1))*EMS(I)*EXP(SUMOP)
       SUMOP=SUMOP+OPACTY(I)
         ENDDO

       SUMOP=0 
         SUMUP=0.
       DO I=NLEV,1,-1
       SUMUP=SUMUP+(TAVG(I)-T(1))*EMS(I)*EXP(SUMOP)
       SUMOP=SUMOP+OPACTY(I)
         ENDDO


       TRAN=EXP(SUMOP)
       TBAVG=(1.-TRAN)*T(1)
       TBDW=TBAVG+SUMDW
       TBUP=TBAVG+SUMUP


         RETURN
         END

c    DIELECTRIC CONSTANT
        SUBROUTINE      meissner(freq,t,s,   eps)
C    COMPLEX DIELECTRIC CONSTANT: EPS
c    T. Meissner, February 2002
c    UPDATED OCT 2004 

C   INPUT:
c   NAME   PARAMETER  UNIT  RANGE
c   FREQ   FREQUENCY  [GHz] 1 to 400
c   T      SST        [K]   248.16 K (-25 C) to 313.16 K (40 C) for pure water
c                           271.16 K (-2  C) to 307.16 K (34 C) for saline water
c   S      SALINITY   [ppt]  0 to 40
c
c   OUTPUT:
c   EPS    COMPLEX DIELECTRIC CONSTANT 
c          NEGATIVE IMAGINARY PART TO BE CONSISTENT WITH WENTZ1 CONVENTIONc
c
c
        USE HDF5                ! This module contains all necessary modules 
        implicit none

        REAL(4) :: FREQ,T,SST,S
        REAL(4) :: E0,E1,E2,N1,N2
        REAL(4) :: A0,A1,A2,B1,B2
        REAL(4) :: E0S,E1S,E2S,N1S,N2S
        REAL(4) :: SIG35,R15,RTR15,ALPHA0,ALPHA1,SIG,F0=17.97510
        COMPLEX(4) :: J = (0.0,1.0), EPS
        INTEGER :: ISTART = 0, len, lenstr

        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
        REAL(4), POINTER :: pt

        CHARACTER datadir*255

        include "aquarius.fin"

c     EXTERNAL ASCII FILE 
c     PATH + FILENAME NEED TO BE SPECIFIED BY USER      
c     UPDATE 2004 

        REAL, TARGET, DIMENSION(11), SAVE ::    X
        REAL, TARGET, DIMENSION(13), SAVE ::    Z
        If (ISTART == 0) then

           call getenv('OCDATAROOT', datadir)
           len = lenstr(datadir)

           pt => x(1)
           data_dims = 11
           call dump_hdf5(datadir(1:len)//'/aquarius/radiometer/rss_dielectric_2004.h5', 'X', data_dims, pt) 

           pt => z(1)
           data_dims = 13
           call dump_hdf5(datadir(1:len)//'/aquarius/radiometer/rss_dielectric_2004.h5', 'Z', data_dims, pt) 

!           open(unit=3,file=inputfile,status='old',form='formatted')
!           READ(3,*) X
!           READ(3,*) Z(1:11)
!           READ(3,*) Z(12:13)
!           CLOSE(3)
        endif

        SST = T - 273.15 ! [Celsius]
        if(sst.lt.-30.16) sst=-30.16 ! Added 10/29/10

c     PURE WATER

      E0    = (3.70886E4 - 8.2168E1*SST)/(4.21854E2 + SST) ! Stogryn et al.
      E1    = X(1) + X(2)*SST + X(3)*SST**2
      N1    = (45.00 + SST)/(X(4) + X(5)*SST + X(6)*SST**2)
      E2    = X(7) + X(8)*SST
      N2    = (45.00 + SST)/(X(9) + X(10)*SST + X(11)*SST**2)



c     Saline Water
c     Conductivity [S/m] taken from Stogryn et al. 
      SIG35 = 2.903602 + 8.60700E-2*SST + 4.738817E-4*SST**2 - 2.9910E-6*SST**3 + 
     1          4.3047E-9*SST**4
      R15   = S*(37.5109+5.45216*S+1.4409E-2*S**2)/(1004.75+182.283*S+S**2)

      alpha0 = (6.9431+3.2841*S-9.9486E-2*S**2)/(84.850+69.024*S+S**2)
      alpha1 = 49.843 - 0.2276*S + 0.198E-2*S**2
      RTR15 = 1.0 + (SST-15.0)*ALPHA0/(ALPHA1+SST)

      SIG = SIG35*R15*RTR15



c     Permittivity
      A0  = exp(Z(1)*S + Z(2)*S**2 + Z(3)*S*SST)
      E0S = A0*E0
      B1  = 1.0 + S*(Z(4) + Z(5)*SST + Z(6)*SST**2)
      N1S = N1*B1
      A1  = exp(Z(7)*S + Z(8)*S**2 + Z(9)*S*SST)
      E1S = E1*A1
      B2 = 1.0 + S*(Z(10) + Z(11)*SST)
      N2S = N2*B2
      A2 = 1.0  + S*(Z(12) + Z(13)*SST)
      E2S = E2*A2


c     Debye Law (2 relaxation wavelengths)
      EPS = (E0S - E1S)/(1.0 - J*(FREQ/N1S)) + 
     1      (E1S - E2S)/(1.0 - J*(FREQ/N2S)) + E2S + 
     2      J*SIG*F0/FREQ


        EPS = CONJG(EPS)

        ISTART = 1


        RETURN 
        END


!      Subroutine FDABSCOEFF(freq,P,T,PV, AV,AO)

c
c       Input:
c     Oxygen Absorption from Rosenkranz
c     freq  frequency [in GHz]
c     P      Pressure [in h Pa]
c     T      Temperature [in K]
c     PV     water vapor pressure  [in hPa]
c
c     Output:   
c     AV          water vapor absorption coefficients [neper/km]
c     AO          oxygen absortption coefficient                [neper/km]

!       implicit none

!      REAL(4) :: FREQ
!      REAL(4) :: P,T,PV,VAP
!       REAL(4) :: AV,AO,ALPHA(3)
!      REAL(4), parameter :: AVOGADRO = 6.02214199E+23,Rd=287.05,
!     1 c=29.97924580,EPSILON=0.622,RV=RD/EPSILON
!       VAP = (1.0E5/RV) * PV / T ! water vapor density [in g /m**3]
!      CALL ABH2O_WENTZ_OLD(P,T,PV,FREQ, ALPHA)
!       AV = alpha(1) + alpha(2) + alpha(3)
!       CALL ABO2_RK(T,P,VAP,FREQ, AO) ! RK O2 ABS 
!      return
!      END


!    Added 10/31/10
c       input:
c     oxygen absorption from rosenkranz
c     freq  frequency [in ghz]
c     p      pressure [in h pa]
c     t      temperature [in k]
c     pv     water vapor pressure  [in hpa]
c
c     output:   
c     av          water vapor absorption coefficients [neper/km]
c     ao          oxygen absortption coefficient                [neper/km]

      subroutine fdabscoeff(freq,p,t,pv, av,ao)
        implicit none

      real(4), parameter :: xnaper=0.2302585094 !convert db/km to naper/km

      real(4) freq,p,t,pv
        real(4) av,ao
        real(4) gamoxy,gamh2o

        call fdabsoxy_1992_modified(p,t,pv,freq, gamoxy)  !gamoxy is db/km
      call abh2o_rk_modified(     p,t,pv,freq, gamh2o)  !gamh2o is db/km

        ao=xnaper*gamoxy
        av=xnaper*gamh2o 
      return
      end

c     following two routines come from 'o:\skytemp2\atms_abs_routines.f' dated June 22, 2009


c     this module contains the oxygen absoprtion routine, the vapor absorption routine and fdaray
c     there were used for the april-june 2009 skytemp update.  see 'memo8.txt'


c     ================================================================================================================
c     ========================== modified version of Liebe 1992 oxygen model =========================================
c     ================================================================================================================

c     This is from Atmospheric 60-GHz Oxygen Spectrum:.. Liebe, Rosenkranz, Hufford, 1992
c     coded: June 2009 1992 by f.wentz
c           inputs: t, temperature (k)
c                   p, total pressure (mb)
c                   pv, water vapor pressure (mb)
c                   freq, frequency (ghz)
c           output: gamoxy, oxygen absorption coefficient (db/km)


c     It is the same as fdabsoxy_1989 except for the a5 and a6 coefs for finding delta have different values.
c     Also this 1992 versions says delta is proprotional to total pressure p rather than pdry. 
c     Also note in this routine, the 1.e-3 scaling is done in the startup block.

c     compared to abo2_rk (the code Rosenkranz sent us), this routine gives very similar results if you set the apterm to 0.
c     for my freqs 6-85 ghz, the largest dif was 0.0003 at the coldest vapor at 85.5 GHz.
c     the apterm adds 0.003 at 85 ghz
c     Apart from the apterm, you can essentially say fdabsoxy_1992 and abo2_rk are the same for my purposes.

c     this routine has been modified in the following ways (see 'memo8.txt')
c     1.  non-resonance continuum temperature coef changed from 0.8 to 1.5
c     2.  a p*p continuum was added
c     these modifications were done June 22 2009 

      subroutine fdabsoxy_1992_modified(p,t,pv,freq, gamoxy)
        implicit none

      integer(4), parameter :: nlines=44
      integer(4) istart,i
      real(4) p,t,pv,freq,gamoxy
        real(4) tht,pwet,pdry,ga,gasq,delta,rnuneg,rnupos,ff,zterm,apterm,sftot,xterm
      real(4) h(6,nlines),f0(nlines),a1(nlines),a2(nlines),a3(nlines),a4(nlines),a5(nlines),a6(nlines)

      real(8) sum

      data istart/1/
        data a4/38*0., 6*0.6/

c          freq          a1      a2       a3       a5          a6
      data h/
     1  50.474238,    0.94e-6,  9.694,  8.60e-3,  0.210,  0.685,
     1  50.987749,    2.46e-6,  8.694,  8.70e-3,  0.190,  0.680,
     1  51.503350,    6.08e-6,  7.744,  8.90e-3,  0.171,  0.673,
     1  52.021410,   14.14e-6,  6.844,  9.20e-3,  0.144,  0.664,
     1  52.542394,   31.02e-6,  6.004,  9.40e-3,  0.118,  0.653,
     1  53.066907,   64.10e-6,  5.224,  9.70e-3,  0.114,  0.621,
     1  53.595749,  124.70e-6,  4.484, 10.00e-3,  0.200,  0.508,
     1  54.130000,  228.00e-6,  3.814, 10.20e-3,  0.291,  0.375,
     1  54.671159,  391.80e-6,  3.194, 10.50e-3,  0.325,  0.265,
     1  55.221367,  631.60e-6,  2.624, 10.79e-3,  0.224,  0.295,
     1  55.783802,  953.50e-6,  2.119, 11.10e-3, -0.144,  0.613,
     1  56.264775,  548.90e-6,  0.015, 16.46e-3,  0.339, -0.098,
     1  56.363389, 1344.00e-6,  1.660, 11.44e-3, -0.258,  0.655,
     1  56.968206, 1763.00e-6,  1.260, 11.81e-3, -0.362,  0.645,
     1  57.612484, 2141.00e-6,  0.915, 12.21e-3, -0.533,  0.606,
     1  58.323877, 2386.00e-6,  0.626, 12.66e-3, -0.178,  0.044,
     1  58.446590, 1457.00e-6,  0.084, 14.49e-3,  0.650, -0.127,
     1  59.164207, 2404.00e-6,  0.391, 13.19e-3, -0.628,  0.231,
     1  59.590983, 2112.00e-6,  0.212, 13.60e-3,  0.665, -0.078,
     1  60.306061, 2124.00e-6,  0.212, 13.82e-3, -0.613,  0.070,
     1  60.434776, 2461.00e-6,  0.391, 12.97e-3,  0.606, -0.282,
     1  61.150560, 2504.00e-6,  0.626, 12.48e-3,  0.090, -0.058,
     1  61.800154, 2298.00e-6,  0.915, 12.07e-3,  0.496, -0.662,
     1  62.411215, 1933.00e-6,  1.260, 11.71e-3,  0.313, -0.676,
     1  62.486260, 1517.00e-6,  0.083, 14.68e-3, -0.433,  0.084,
     1  62.997977, 1503.00e-6,  1.665, 11.39e-3,  0.208, -0.668,
     1  63.568518, 1087.00e-6,  2.115, 11.08e-3,  0.094, -0.614,
     1  64.127767,  733.50e-6,  2.620, 10.78e-3, -0.270, -0.289,
     1  64.678903,  463.50e-6,  3.195, 10.50e-3, -0.366, -0.259,
     1  65.224071,  274.80e-6,  3.815, 10.20e-3, -0.326, -0.368,
     1  65.764772,  153.00e-6,  4.485, 10.00e-3, -0.232, -0.500,
     1  66.302091,   80.09e-6,  5.225,  9.70e-3, -0.146, -0.609,
     1  66.836830,   39.46e-6,  6.005,  9.40e-3, -0.147, -0.639,
     1  67.369598,   18.32e-6,  6.845,  9.20e-3, -0.174, -0.647,
     1  67.900867,    8.01e-6,  7.745,  8.90e-3, -0.198, -0.655,
     1  68.431005,    3.30e-6,  8.695,  8.70e-3, -0.210, -0.660,
     1  68.960311,    1.28e-6,  9.695,  8.60e-3, -0.220, -0.665,
     1 118.750343,  945.00e-6,  0.009, 16.30e-3, -0.031,  0.008,
     1 368.498350,   67.90e-6,  0.049, 19.20e-3,  0.0,    0.0,
     1 424.763124,  638.00e-6,  0.044, 19.16e-3,  0.0,    0.0,
     1 487.249370,  235.00e-6,  0.049, 19.20e-3,  0.0,    0.0,
     1 715.393150,   99.60e-6,  0.145, 18.10e-3,  0.0,    0.0,
     1 773.839675,  671.00e-6,  0.130, 18.10e-3,  0.0,    0.0,
     1 834.145330,  180.00e-6,  0.147, 18.10e-3,  0.0,    0.0/

      if(istart.eq.1) then
      istart=0
      f0(:)=h(1,:)
      a1(:)=h(2,:)/h(1,:)
      a2(:)=h(3,:)
      a3(:)=h(4,:)
      a5(:)=0.001*h(5,:)
      a6(:)=0.001*h(6,:)
      endif

      tht = 300/t
      pwet=0.1*pv
      pdry=0.1*p-pwet
      xterm=1-tht

      sum = 0.
      do i=1,nlines
      ga = a3(i)*(pdry*tht**(0.8-a4(i)) + 1.1*tht*pwet)
      gasq=ga*ga
        delta=(a5(i) + a6(i)*tht)*p*tht**0.8
      rnuneg = f0(i)-freq
      rnupos = f0(i)+freq
      ff = (ga-rnuneg*delta)/(gasq+rnuneg**2) +  (ga-rnupos*delta)/(gasq+rnupos**2)
      sum = sum + ff*a1(i)*exp(a2(i)*xterm)
        enddo
      if(sum.lt.0) sum=0

c     add nonresonant contribution

c     ga=5.6e-3*(pdry+1.1*pwet)*tht**0.8  
      ga=5.6e-3*(pdry+1.1*pwet)*tht**1.5  !modification 1
      zterm=ga*(1.+(freq/ga)**2)
      apterm=1.4e-10*(1-1.2e-5*freq**1.5)*pdry*tht**1.5
      if(apterm.lt.0) apterm=0
cx    sftot=pdry*freq*tht**2 * (tht*sum + 6.14e-4/zterm + apterm)
      sftot=pdry*freq*tht**2 * (tht*sum + 6.14e-4/zterm + apterm + 3.e-10*pdry*tht) !modification 2
      gamoxy=0.1820*freq*sftot
      return
      end


c     ================================================================================================================
c     ========================== modified version of Rosenkranz water vapor model ====================================
c     ================================================================================================================


c purpose- compute absorption coef in atmosphere due to water vapor
c
c  calling sequence parameters-
c    specifications
c      name    units    i/o  descripton            valid range
c      t       kelvin    i   temperature
c      p       millibar  i   pressure              .1 to 1000
c      f       ghz       i   frequency             0 to 800
c      gamh2o  db/km     o   absorption coefficient
c
c   references-
c    p.w. rosenkranz, radio science v.33, pp.919-928 (1998).
c
c   line intensities selection threshold=
c     half of continuum absorption at 1000 mb.
c   widths measured at 22,183,380 ghz, others calculated.
c     a.bauer et al.asa workshop (sept. 1989) (380ghz).
c
c   revision history-
c    date- oct.6, 1988  p.w.rosenkranz - eqs as publ. in 1993.
c          oct.4, 1995  pwr- use clough's definition of local line
c                   contribution,  hitran intensities, add 7 lines.
c          oct. 24, 95  pwr -add 1 line.
c          july 7, 97   pwr -separate coeff. for self-broadening,
c                       revised continuum.
c          dec. 11, 98  pwr - added comments

c     the routine is a modified version of abh2o_rk_reformat. 
c     this routine has been modified in the following three ways (see 'memo8.txt')
c     1.  b1(1)=1.01*b1(1)  :22 ghz line strength increase slightly
c     2.  22 ghz line shape below 22 ghz has been modified
c     3.  foreign and self broadening continuum has been adjusted 
c     these modification were done June 22 2009 
 




      subroutine abh2o_rk_modified(p,t,pv,freq,  gamh2o)
      implicit none
  
      integer(4), parameter :: nlines=15

      integer(4) istart,i
      real(4) t,p,freq
      real(4) b1(nlines),b2(nlines),b3(nlines),f0(nlines),b4(nlines),b5(nlines),b6(nlines)
      real(4) pv,s,base,gamh2o
        real(4) tht,pwet,pdry,ga,gasq,sftot,xterm,rnuneg,rnupos
        real(8) sum

        real(4) chi,chisq,freqsq,f0sq,u

      data istart/1/

c     line frequencies:
      data f0/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,
     & 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,
     & 620.7008, 752.0332, 916.1712/
c     line intensities at 300k:
      data b1/ .1310e-13, .2273e-11, .8036e-13, .2694e-11, .2438e-10,
     &         .2179e-11, .4624e-12, .2562e-10, .8369e-12, .3263e-11, .6659e-12,
     &         .1531e-08, .1707e-10, .1011e-08, .4227e-10/
c     t coeff. of intensities:
      data b2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
c     air-broadened width parameters at 300k:
      data b3/.0281, .0281, .023, .0278, .0287, .021, .0186, .0263, .0215, .0236, .026, .0321, .0244, .0306, .0267/
c     self-broadened width parameters at 300k:
      data b5/.1349, .1491, .108, .135, .1541, .090, .0788, .1275, .0983, .1095, .1313, .1320, .1140, .1253, .1275/
c     t-exponent of air-broadening:
      data b4 /.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, .71, .68, .70/
c     t-exponent of self-broadening:
      data b6/.61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72, 1.0, .68, .84, .78/
c
      if(istart.eq.1) then
        istart=0
      b1=1.8281089E+14*b1/f0**2
        b5=b5/b3  !convert b5 to Leibe notation
      b1(1)=1.01*b1(1)  !modification 1
        endif
 
      if(pv.le.0.) then
      gamh2o=0
        return
      endif
 

      pwet=0.1*pv
      pdry=0.1*p-pwet
      tht = 300./t
      xterm=1-tht
        freqsq=freq*freq

      sum = 0.
      do i=1,nlines
        f0sq=f0(i)*f0(i)
      ga=b3(i)*(pdry*tht**b4(i) + b5(i)*pwet*tht**b6(i))
      gasq = ga*ga
      s = b1(i)*exp(b2(i)*xterm)
      rnuneg = f0(i)-freq
      rnupos = f0(i)+freq
      base = ga/(562500. + gasq)  !use clough's definition of local line contribution

      if(i.ne.1) then
      if(abs(rnuneg).lt.750) sum = sum + s*(ga/(gasq + rnuneg**2) - base)
      if(abs(rnupos).lt.750) sum = sum + s*(ga/(gasq + rnupos**2) - base)

      else
        chi=0
        if(freq.lt.f0(i)) then
        u=abs(freq-f0(i))/f0(i)
        if(u.gt.1) u=1
        chi=ga*u*u*(3-2*u)  !modification 2
      endif
        chisq=chi*chi
      sum=sum +     s*2*((ga-chi)*freqsq + (ga+chi)*(f0sq+gasq-chisq))/((freqsq-f0sq-gasq+chisq)**2 + 4*freqsq*gasq)
        endif

        enddo
      if(sum.lt.0) sum=0

cx    sftot=pwet*freq*tht**3.5*(sum +     1.2957246e-6*pdry/tht**0.5 +                    4.2952193e-5*pwet*tht**4)
      sftot=pwet*freq*tht**3.5*(sum + 1.1*1.2957246e-6*pdry/tht**0.5 + 0.375*(freq**0.15)*4.2952193e-5*pwet*tht**4) !modification 3

      gamh2o=0.1820*freq*sftot
        return
      end

      SUBROUTINE ABO2_RK(TEMP,PRES,VAPDEN,FREQ,  O2ABS)
C
C     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
C              IN NEPERS/KM
C
C      5/1/95  P. Rosenkranz 
C      11/5/97  P. Rosenkranz - 1- line modification.
c      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
C
C     ARGUMENTS:
      REAL TEMP,PRES,VAPDEN,FREQ, O2ABS
C
C     NAME    UNITS    DESCRIPTION        VALID RANGE
C
C     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
c                                          valid for atmosphere
C     PRES   MILLIBARS PRESSURE           3 TO 1000
C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
C     FREQ    GHZ      FREQUENCY          0 TO 900
C
C     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
C     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
C      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
C     H.J. Liebe et al, JQSRT V.48, PP.629-643 (1992).
c     M.J. Schwartz, Ph.D. thesis, M.I.T. (1997).
C     SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
c     This version differs from Liebe's MPM92 in two significant respects:
c     1. It uses the modification of the 1- line width temperature dependence
c     recommended by Schwartz: (1/T).
c     2. It uses the same temperature dependence (X) for submillimeter 
c     line widths as in the 60 GHz band: (1/T)**0.8 

C
      REAL(4), parameter :: Rd=287.05,EPSILON=0.622,RV=RD/EPSILON

      real(8) sum

      COMMON /O2COM/ X,WB300,W300(40),F(40),Y300(40),S300(40),
     & V(40),BE(40)
C      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      DATA F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,
     2  59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
     3  56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,
     4  55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,
     5  53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,
     6  52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632,
     7  487.2494, 715.3931, 773.8397, 834.1458/
        DATA S300/.2936E-14,.8079E-15, .2480E-14,.2228E-14,
     &  .3351E-14,.3292E-14, .3721E-14,.3891E-14,
     &  .3640E-14,.4005E-14, .3227E-14,.3715E-14,
     &  .2627E-14,.3156E-14, .1982E-14,.2477E-14,
     &  .1391E-14,.1808E-14, .9124E-15,.1230E-14,
     &  .5603E-15,.7842E-15, .3228E-15,.4689E-15,
     &  .1748E-15,.2632E-15, .8898E-16,.1389E-15,
     &  .4264E-16,.6899E-16, .1924E-16,.3229E-16,
     &  .8191E-17,.1423E-16, .6494E-15, .7083E-14, .3025E-14,
     &  .1835E-14, .1158E-13, .3993E-14/
      DATA BE/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626,
     & 2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625,
     & 2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844,
     & 2*7.744, .048, .044, .049, .145, .141, .145/
C      WIDTHS IN MHZ/MB
      DATA WB300/.56/, X/.8/
      DATA W300/1.63, 1.646, 1.468, 1.449, 1.382, 1.360,
     & 1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,
     & 1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05,
     & 2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.92, 3*1.81/
      DATA Y300/  -0.0233,  0.2408, -0.3486,  0.5227,
     & -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,
     &  0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,
     &  0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,
     &  0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,
     &  0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      DATA V/  0.0079, -0.0978,  0.0844, -0.1273,
     &  0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,
     &  0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,
     &  0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,
     &  0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,
     &  0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
C
      TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      PRESWV = VAPDEN*TEMP*RV*1.0E-5

      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DENS = .001*(PRESDA + 1.1*PRESWV)*TH
      DFNR = WB300*DEN
      SUM = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      DO 32 K=1,40
      IF(K.EQ.1) THEN !exception for 1- line
        DF = W300(1)*DENS
      ELSE
        DF = W300(K)*DEN
      ENDIF
      Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
      STR = S300(K)*EXP(-BE(K)*TH1)
      SF1 = (DF + (FREQ-F(K))*Y)/((FREQ-F(K))**2 + DF*DF)
      SF2 = (DF - (FREQ+F(K))*Y)/((FREQ+F(K))**2 + DF*DF)
32    SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))**2
      O2ABS = .5034E12*SUM*PRESDA*TH**3/3.14159
      RETURN
      END


c     Wentz Water Vapor (BEFORE March 2002 Version) 
      SUBROUTINE ABH2O_WENTZ_OLD(P,T,PV,FREQ, ALPHA)
c     same as ABH2O_NEW, but FRQTERM = 1 (no frequency dependence)
c     identical to sky_routines
C
C  FROM LIEBE RADIO SCI., VOL 20, NO 5, PP 1069-1089 SEPT-OCT 1985
C      CODED: DEC 1992 BY F.WENTZ
c
c  corrections:
C  self broadening term of the water vapor continuum( CTERM1b) term is multiplied by 0.52 
c  (ATBD: page 15) and the width (i.e. B3 term) of the 22 GHz line is multiplied by 1.015
c
C           INPUTS: T, TEMPERATURE (K)
C                   P, TOTAL PRESSURE (MB)
C                   PV, WATER VAPOR PRESSURE (MB)
C                   FREQ, FREQUENCY (GHZ)
C           OUTPUT: ALPHA: WATER VAPOR ABSORPTION COEFFICIENTS (NEPER/KM)
c               1 = LINE,  2= FB,  3 = SB
C
      REAL*4 F0(3),B1(3),B2(3),B3(3),ALPHA(3)
      DATA F0/22.235080,183.310117,325.152919/
      DATA B1/0.1090,2.3000,1.5400/
      DATA B2/2.143,0.653,1.515/
CCC   DATA B3/27.84E-3,28.35E-3,27.00E-3/
      DATA B3/28.26E-3,28.35E-3,27.00E-3/
      DATA ISTART/1/
C
      IF(ISTART.EQ.1) THEN
      ISTART=0
      B1(1)=B1(1)/F0(1)
      B1(2)=B1(2)/F0(2)
      B1(3)=B1(3)/F0(3)
      ENDIF
C
      THT=300/T
      PWET=0.1*PV
      PDRY=0.1*P-PWET
      TTERM=THT**.8
      PTERM=PDRY*TTERM + 4.8*THT*PWET
      XTERM=1-THT
C
      SUM=0.
      DO 100 I=1,3
      GA=B3(I)*PTERM
      GASQ=GA*GA
      FF=GA/(GASQ+(F0(I)-FREQ)**2) + GA/(GASQ+(F0(I)+FREQ)**2)
      SUM=SUM+FF*B1(I)*EXP(B2(I)*XTERM)
  100 CONTINUE
C
      CTERM1A=1.4E-6*PDRY/THT
      CTERM1B= 0.52*5.41E-5*PWET*THT**2
      CTERM2= 2.9E-8*PDRY*PWET**0.1*SQRT(FREQ)/THT**1.5
        FRQTERM=1.
        CFB = CTERM1A
        CSB = CTERM1B*FRQTERM
      SFTOT= PWET*FREQ*THT**3.5*(SUM + CTERM1A + CTERM1B*FRQTERM + CTERM2)

      GAMH2O=0.1820*FREQ*SFTOT
      XNAPER =0.2302585094 !CONVERT DB/KM TO NEPER/KM

        alpha(1) = XNAPER*0.1820*FREQ*PWET*FREQ*(THT**3.5)*(SUM+CTERM2)
        alpha(2) = XNAPER*0.1820*FREQ*PWET*FREQ*(THT**3.5)*CFB
        alpha(3) = XNAPER*0.1820*FREQ*PWET*FREQ*(THT**3.5)*CSB

      RETURN
      END
        SUBROUTINE FDCLDABS(freq,T,RHOL,   AL)
c       USE RSS_SURFACE_TB ! for dielectric constant

c    liquid cloud water absorption
c    Rayleigh 
c    freq:      Frequency [GHz]
c    T:         Temperature [K]
c    RHOL:      liquid cloud water density [g/m**3]

c    Output:
c    AL:        cloud water absorption coefficient [neper/km]

      implicit none 
        real(4) :: freq,t,rhol,rhol0,AL,WAVLEN,C=29.979,PI=3.14159
!       c changed from 29.79 to 29.979 10/31/10

        complex(4) :: permit

        RHOL0 = 1.0E-6*RHOL ![g/cm**3]
c       CALL MEISSNER_WENTZ_2004_DIELECTRIC(FREQ,T,0.0, PERMIT)
        CALL meissner(freq,t,0.0,   permit)     
        WAVLEN = C/FREQ
        AL = (6.0*PI*RHOL0/WAVLEN)*AIMAG((1.0-PERMIT)/(2.0+PERMIT))  ! nepers/cm
        AL = 1.0E5*AL ! nepers/km
        return
        end

      subroutine FINDNCEP(NMAX,xlat,xlon,COLVAP,COLWAT,
     1  PWAT,CWAT,P,T,PV,RHOV,RHOL,z,ibegin,metfile,ierr)
c
c
c     01/05/2006
c     WAS: IF T(ipr) < 100 K : ierr=-10, bail out
c     CHANGED: IF T(ipr) < 100 K, set RHOCWAT=0.0
c
c     01/28/2003
c     RHOL = liquid cloud water density 
c
c     ordered NCEP profiles, SURTEP, WIND, PWAT ,COLVAP
c     
c     Input:
c     nmax: maximum number of levels
c     yy:   long year (e.g. 1999)
c     mm:   month
c     isec: UCT second
c     xlat: Latitude  [between -90 and 90]
c     xlon: Longitude   [between 0 and 360]


c    Output:
c     COLVAP: NCEP water vapor integrated [mm]
c     COLWAT: NCEP columnar liquid cloud water integrated [mm]
c     PWAT:   NCEP value for precipitable water [mm]
c     CWAT:   NCEP value for columnar cloud water       (liquid + ice) [mm]
c     P: air pressure profile      [mb] (0:nmax)   ordered
c     T: air temperature profile [K]  (0:nmax)
c     Z: elevation profile [m]
c     PV: water vapor pressure profile [mb]  (0:nmax)
c     RHOV: water vapor density [g /m**3]
c     RHOL: liquid water density [g /m**3] 
c     ibegin: index of surface level
c     ierr   =  0: valid
c            /= 0: invalid  

      implicit none

      integer(4):: nmax,nrh
      integer(4) :: ibegin,ierr,ipr
      real(4)    :: xlat,xlon,colvap,colwat,pwat,cwat,
     $     p_sfc,P0,t2M,RHOAIR
      real(4), dimension(0:nmax) ::  T,RH,HGT,z,P,CLWMR
      real(4), dimension(0:nmax) ::  PV,RHOV,RHOCWAT,RHOL,RHOI
      real(4), dimension(0:nmax) ::    wksp
      integer(4), dimension(0:nmax) :: iwksp
      real(4), parameter  :: R_E = 6371000, RD=287.05, EPSILON = 0.622 
      CHARACTER(len=128) :: metfile

      NRH = 21
      RH(nrh+1:nmax) = 0.0
      CLWMR(nrh+1:nmax) = 0.0

c     Only read a single qmet ancillary file
      call FND_NCEP_XX(xlat,xlon,
     $        t,hgt,rh,clwmr,p_sfc,pwat,cwat,metfile,ierr)

        if (ierr .ne. 0) return

        P0 = 0.01*p_sfc    ! hPa
        T2M = T(0)

c    sort according to prs in descending order
c    if surface pressure >= 1000:        level (0,1,2,...) = (sfc,1000,975,...)
c    if 975 <= surface pressure <= 1000: level (0,1,2,...) = (1000,sfc,975,...)
c    etc.
c
c
      P(0:nmax) = (/   P0,
     1         1000.,975.,950.,925.,900.,850.,800.,750.,700.,650.,600.,
     2       550.,500.,450.,400.,350.,300.,250.,200.,150.,100.,
     3           70., 50., 30., 20. ,10.  /)
c


      P = -P
        call sort5(nmax+1,P,t,rh,clwmr,hgt,wksp,iwksp)
        P=-P    
        iwksp = iwksp-1

c
      z = HGT * R_E / (R_E - HGT)

c     find index for surface
        do ipr =0,nmax
        if (iwksp(ipr).eq.0) ibegin = ipr
      enddo
        if (z(ibegin) >= z(ibegin+1)) z(ibegin) = z(ibegin+1) - 0.1 
        if(ABS(P0-P(IBEGIN)) > 1.0E-3 *P_SFC) THEN
        WRITE(*,*) 'INCOSISTENCY IN PRESSURE'
        WRITE(*,*) P0,IBEGIN,P(IBEGIN)
        DO IPR = 0,NMAX
        WRITE(*,*) IPR,P(IPR)
        ENDDO
      STOP
        ENDIF


c       transform RH -> water vapor pressure and density
      call goff_gratch_vap(nmax,T,RH,P,  PV,rhoV)

c     convert NCEP clwmr to RHOCWAT
      clwmr(ibegin) = clwmr(ibegin+1) ! clwmr at surface  
        do ipr = 0,nmax
        IF (T(IPR) >= 100) THEN
        RHOAIR = P(ipr)/(RD*T(IPR)) ! density of dry air
        RHOAIR = RHOAIR*(1.0 - (1.0-EPSILON)*PV(IPR)/P(IPR)) 
        ! density of moist air (without cloud), unusual NCEP definition for mixing ratio

      RHOCWAT(ipr) = 1.0E5*CLWMR(ipr)*RHOAIR    ! [g/m**3]
        CALL CLDWATICE(T(IPR),RHOCWAT(IPR),  RHOL(IPR),RHOI(IPR))  
        ELSE ! unphysical temperature
      RHOCWAT(IPR)=0.0
        RHOL(IPR) = 0.0
      RHOI(IPR) = 0.0
        ENDIF

c

        enddo

      where (rhov < 0) RHOV=0.0
      where (rhol < 0) RHOL=0.0

      call column(nmax-ibegin,z(ibegin:nmax),rhoV(ibegin:nmax),2,
     1                                                       COLVAP,ierr)
      call column(nmax-ibegin,z(ibegin:nmax),RHOL(ibegin:nmax),1,
     1                                                       COLWAT,ierr)


        COLVAP = COLVAP *1.E-3   ! g/m**2 -> mm = kg/m**2
        COLWAT = COLWAT *1.E-3   ! g/m**2 -> mm = kg/m**2
        end




        SUBROUTINE CLDWATICE(T,RHOC,  RHOL,RHOI)
c     T: TEMPERATURE [K]
c     RHOC:   cloud water density [arbitrary unit]
c 
c     RHOL : liquid part [same unit as RHOC]
c     RHOI : ice    part [same unit as RHOC]
      implicit none
        real, parameter :: TWAT = 273.15, TICE = 253.15
        real :: T,RHOC,RHOL,RHOI

c     water or ice 
        if (T >= TWAT) then
                 RHOL = RHOC ! all liquid water
               RHOI = 0.0 
        else if (T <= TICE) then 
                 RHOL = 0.0       ! all ice
               RHOI = RHOC
        else
                 RHOL = RHOC*(T-TICE)/(TWAT-TICE)
               RHOI = RHOC - RHOL  
                         ! linear temperature interpolation in between
        endif
        RETURN
        END


      SUBROUTINE FND_NCEP_XX(xlat,xlon,
     &                     t,hgt,rh,clwmr,p_sfc,pwat,cwat,metfile,ierr)

c
c    NCEP variables
c    1 deg resolution
c    4 times daily: 00Z 06Z 12Z 18Z
c 
c    Thomas Meissner : APR 1999 
c
c    needs to be linked with time_routines.f 
c
c
c         Varaible    Level (subfolder)
c           SST               sfc    (0) 
c           T                 sfc    (0) 2 m above ground
c           RH                sfc    (0)
c           HGT               sfc    (0)
c           P                 sfc    (0)
c           PWAT              col
c           CWAT              col
c           T             1:nmax  (prf)
c           RH            1:nrh   (prf)
c           HGT           1:nmax  (prf)
c           CLWMR         1:nmax  (prf) 
c
c    named as year_month_day_xxZ.dat :
c    binary file:
c    year month day hour (integer(4))
c    real*4 array 
c    GRID: LAT from +90 to -90 by -1 deg
c          LON from 0 to 360 by 1 deg
c    total size 4*4 + 360*181*4 = 260656  
c
c           xlat   (latitude)                                     REAL(4)
c           xlon   (longitude)                                    REAL(4)
c    
c    Output: var (interpolated variable) REAL(4) 
c              IERR error parameter =0,  regular return
c                                 =-1, no data (var set to -999)
c                                   
c
      integer, parameter :: nmax =26, nrh =21
      real, save :: T1(360,181,0:nmax), RH1(360,181,0:nrh),
     &              HGT1(360,181,0:nmax),CLWMR1(360,181,0:nrh),
     &              p_sfc1(360,181),pwat1(360,181),cwat1(360,181)
      CHARACTER(len=128) :: metfile
c
      integer, save :: ifirst=1
c
      REAL t(0:nmax),rh(0:nrh),hgt(0:nmax),clwmr(0:nmax),
     &p_sfc,pwat,cwat

      if (ifirst.eq.1) then
         CALL yread(T1,RH1,CLWMR1,HGT1,P_sfc1,pwat1,cwat1,metfile)
         ifirst = 0
      endif

      A1=1
C                                                                       
      IF(XLON.LT.0..OR.XLON.GT.360.)    then
         ierr = -4
         return
      endif
                                  
      IF(ABS(XLAT).GT.90.) then
         ierr = -5
         return
      endif


      BRIEF= -(XLAT-90.)
      J1=INT(1+BRIEF)                                                        
      J2=J1+1 
      if(j2.ge.181) j2=181                                                           
      B1=J1-BRIEF                                                       
      B2=1.-B1
c 
      BRIEF=XLON
      IF(BRIEF.GT.359.999) BRIEF=0.
      K1=INT(1+BRIEF)                                                       
      K2=K1+1                                                           
      IF(K2.EQ.361) K2=1                                                
      C1=K1-BRIEF                                                       
      C2=1.-C1     
C
      T(0:nmax)=                                                              
     1     A1*B1*(C1*T1(K1,J1,0:nmax)+C2*T1(K2,J1,0:nmax))+                 
     2     A1*B2*(C1*T1(K1,J2,0:nmax)+C2*T1(K2,J2,0:nmax))
C
      RH(0:nrh)=                                                              
     1     A1*B1*(C1*RH1(K1,J1,0:nrh)+C2*RH1(K2,J1,0:nrh))+                 
     2     A1*B2*(C1*RH1(K1,J2,0:nrh)+C2*RH1(K2,J2,0:nrh))
C     
      CLWMR(0:nrh)=                                                              
     1     A1*B1*(C1*CLWMR1(K1,J1,0:nrh)+C2*CLWMR1(K2,J1,0:nrh))+                 
     2     A1*B2*(C1*CLWMR1(K1,J2,0:nrh)+C2*CLWMR1(K2,J2,0:nrh))

C
      HGT(0:nmax)=                                                              
     1     A1*B1*(C1*HGT1(K1,J1,0:nmax)+C2*HGT1(K2,J1,0:nmax))+                 
     2     A1*B2*(C1*HGT1(K1,J2,0:nmax)+C2*HGT1(K2,J2,0:nmax))
c  
      p_sfc =                                                              
     1     A1*B1*(C1*p_sfc1(K1,J1)+C2*p_sfc1(K2,J1))+                 
     2     A1*B2*(C1*p_sfc1(K1,J2)+C2*p_sfc1(K2,J2))
c  
      pwat =                                                              
     1     A1*B1*(C1*pwat1(K1,J1)+C2*pwat1(K2,J1))+                 
     2     A1*B2*(C1*pwat1(K1,J2)+C2*pwat1(K2,J2))
c  
      cwat =                                                              
     1     A1*B1*(C1*cwat1(K1,J1)+C2*cwat1(K2,J1))+                 
     2     A1*B2*(C1*cwat1(K1,J2)+C2*cwat1(K2,J2))
      
c 
      IERR=0
      RETURN
      END


c    01/23/2006
c    KYLE CHANGED SUBROUTINE YREAD
c    read fileheader between each profile level
c    



      SUBROUTINE yread(T,RH,CLWMR,HGT,P_sfc,pwat,cwat,metfile)
c
      integer , parameter :: nmax = 26, nrh=21
      real :: T(360,181,0:nmax), RH(360,181,0:nrh),CLWMR(360,181,0:nrh),
     &              HGT(360,181,0:nmax),
     &              p_sfc(360,181),pwat(360,181),CWAT(360,181)

      CHARACTER*128 metfile
      INTEGER*4 ILEVEL, irow
      integer, parameter :: DFACC_READ = 1
      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(3), count(3), stride(3)
      integer retn 
      real flt(360)

      stride(1) = 1
      stride(2) = 1
      stride(3) = 1

      sd_id_anc = sfstart( metfile, DFACC_READ)
      if (sd_id_anc.eq.-1) then
         write(*,*) 'Unable to open: '//metfile
         stop
      endif

      start(1) = 0
      start(2) = 0
      count(1) = 360
      count(2) = 181
      count(3) = 1

c     Temp
      sds_index = sfn2index(sd_id_anc, 'tmp_prfl')
      if (sds_index.eq.-1) then
         write(*,*) '"tmp_prfl" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)

      do ilevel=1,nmax
         start(3) = ilevel-1
         retn = sfrdata(sds_id, start, stride, count, T(1,1,ILEVEL))
         do irow=1,90
            flt = T(:,irow,ilevel)
            T(:,irow,ilevel) = T(:,181-irow+1,ilevel)
            T(:,181-irow+1,ilevel) = flt
         enddo
      enddo
      retn = sfendacc(sds_id)

      sds_index = sfn2index(sd_id_anc, 'tmp_sfc')
      if (sds_index.eq.-1) then
         write(*,*) '"tmp_sfc" field not found in: '//metfile
         stop
      endif
      
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, T(1,1,0))
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = T(:,irow,0)
         T(:,irow,0) = T(:,181-irow+1,0)
         T(:,181-irow+1,0) = flt
      enddo
c      write(*,*) t(100,85,0),t(100,120,0)
c      write(*,*) t(100,85,1),t(100,120,1)

c     RH
      sds_index = sfn2index(sd_id_anc, 'rh_prfl')
      if (sds_index.eq.-1) then
         write(*,*) '"rh_prfl" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)

      do ilevel=1,nrh
         start(3) = ilevel-1
         retn = sfrdata(sds_id, start, stride, count, RH(1,1,ILEVEL))
         do irow=1,90
            flt = RH(:,irow,ilevel)
            RH(:,irow,ilevel) = RH(:,181-irow+1,ilevel)
            RH(:,181-irow+1,ilevel) = flt
         enddo
      enddo
      retn = sfendacc(sds_id)

      sds_index = sfn2index(sd_id_anc, 'rh')
      if (sds_index.eq.-1) then
         write(*,*) '"rh" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, RH(1,1,0))
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = RH(:,irow,0)
         RH(:,irow,0) = RH(:,181-irow+1,0)
         RH(:,181-irow+1,0) = flt
      enddo


c     HGT
      sds_index = sfn2index(sd_id_anc, 'geohgt_prfl')
      if (sds_index.eq.-1) then
         write(*,*) '"geohgt_prfl" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)

      do ilevel=1,nmax
         start(3) = ilevel-1
         retn = sfrdata(sds_id, start, stride, count, HGT(1,1,ILEVEL))
         do irow=1,90
            flt = HGT(:,irow,ilevel)
            HGT(:,irow,ilevel) = HGT(:,181-irow+1,ilevel)
            HGT(:,181-irow+1,ilevel) = flt
         enddo
      enddo
      retn = sfendacc(sds_id)

      sds_index = sfn2index(sd_id_anc, 'geohgt')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, HGT(1,1,0))
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = HGT(:,irow,0)
         HGT(:,irow,0) = HGT(:,181-irow+1,0)
         HGT(:,181-irow+1,0) = flt
      enddo

c     CLWMR
      sds_index = sfn2index(sd_id_anc, 'clwmr_prfl')
      if (sds_index.eq.-1) then
         write(*,*) '"clwmr_prfl" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)

      do ilevel=1,nrh
         start(3) = ilevel-1
         retn = sfrdata(sds_id, start, stride, count, CLWMR(1,1,ILEVEL))
         do irow=1,90
            flt = CLWMR(:,irow,ilevel)
            CLWMR(:,irow,ilevel) = CLWMR(:,181-irow+1,ilevel)
            CLWMR(:,181-irow+1,ilevel) = flt
         enddo
      enddo
      retn = sfendacc(sds_id)
      CLWMR(1:360,1:181,0) = 0.0 !!!!! Set to 0


c     Surface Pressure
      sds_index = sfn2index(sd_id_anc, 'press')
      if (sds_index.eq.-1) then
         write(*,*) '"press" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, P_SFC)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = P_SFC(:,irow)
         P_SFC(:,irow) = P_SFC(:,181-irow+1)
         P_SFC(:,181-irow+1) = flt
      enddo
c      write(*,*) metfile
c      write(*,*) p_sfc(1,:)

c     P water
      sds_index = sfn2index(sd_id_anc, 'p_water')
      if (sds_index.eq.-1) then
         write(*,*) '"p_water" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, PWAT)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = PWAT(:,irow)
         PWAT(:,irow) = PWAT(:,181-irow+1)
         PWAT(:,181-irow+1) = flt
      enddo


c     C water
      sds_index = sfn2index(sd_id_anc, 'c_water')
      if (sds_index.eq.-1) then
         write(*,*) '"c_water" field not found in: '//metfile
         stop
      endif

      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, CWAT)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = CWAT(:,irow)
         CWAT(:,irow) = CWAT(:,181-irow+1)
         CWAT(:,181-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 
c
      RETURN
      END 


      subroutine goff_gratch_vap(nprofiles,T,RH,P,  P_V,rho_V)
c
c    
c     calculates water vapor pressure P_V and density rho_V
c     from relative humidity RH
c     S. Cruz Pol, C. Ruf and S. Keihm, Radio Science 33 (5),1319 (1998)
c 
c    water vapor pressure: P_v = P_s * RH
c    water vapor density   rho_v = F_w *(P_v * eps) / (R_d * T)
c
      real, dimension(0:nprofiles) :: T,RH,P,P_v,rho_v,P_s
        real, dimension(0:nprofiles) :: F_w,xi1,xi2,xi3  
        real  :: p_standard = 1013.246  ! hPa
        real  :: TS = 373.14     ! K  (water boiling)
        real  :: T0 = 273.14  ! K 
        real  :: eps = 1./1.607795  ! M(H2O)/M(air)
        real  :: R_d = 287.05 ! J /(K*kg)  gas constant for dry air
c
      F_w = 1.0 + 1.E-4 * (5.92854 + 3.740346e-2 * P +
     &              1.971198E-4 * (T-T0) * (800-P)  +
     &                    6.045511E-6 * P * (T-T0)**2 ) ! deviation from ideal gas
c
      xi1 = -7.90298*(Ts/T - 1) + 5.02808*Log10(TS/T)
        xi2 = -1.3816E-7 * 10**(11.344*(1.-T/TS) -1 )
        xi3 = 8.1328E-3 * (10**(-3.49149*(TS/T - 1)) - 1)
c
      P_s = p_standard * 10**(xi1 + xi2 + xi3)
        P_v = P_s * RH* 0.01                                                                               !  mbar
        rho_v = ((F_w * P_v * eps) / (R_d * T)) * 1.E5  !  [g/m**3]
c
      return
        end  


      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
!        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END



      subroutine vapor(nlev,T,RH,P,  PV,vap)
c
c     calculates water vapor pressure PV and density vap
c     from relative humidity RH
c     S. Cruz Pol, C. Ruf and S. Keihm, Radio Science 33 (5),1319 (1998)
c 
c    water vapor pressure: PV = P_s * RH
c    water vapor density  vap = F_w *(PV * eps) / (R_d * T)
c
      real, dimension(0:nlev) :: T,RH,P,PV,vap,P_s
        real, dimension(0:nlev) :: F_w,xi1,xi2,xi3  
        real  :: p_standard = 1013.246  ! hPa
        real  :: TS = 373.14     ! K  (boiling point)
        real  :: T0 = 273.14     ! K  (freezing point)
        real  :: eps = 1./1.607795  ! M(H2O)/M(air)
        real  :: R_d = 287.05 ! J /(K*kg)  gas constant for dry air
c
      F_w = 1.0 + 1.E-4 * (5.92854 + 3.740346e-2 * P +
     1                1.971198E-4 * (T-T0) * (800-P)  +
     2                    6.045511E-6 * P * (T-T0)**2 ) ! deviation from ideal gas
c
      xi1 = -7.90298*(Ts/T - 1) + 5.02808*Log10(TS/T)
        xi2 = -1.3816E-7 * 10**(11.344*(1.-T/TS) -1 )
        xi3 = 8.1328E-3 * (10**(-3.49149*(TS/T - 1)) - 1)
c
      P_s = p_standard * 10**(xi1 + xi2 + xi3)
        PV = P_s * RH* 0.01                                                                                !  mbar
        vap = ((F_w * PV * eps) / (R_d * T)) * 1.E5  !  [g/m**3]
c
      return
        end  


      subroutine column(nlevel,z,rho,ip,  COL,ierr)
c     columnar integral
c     input:  
c     nlevel: number of profiles (nlevel = 0 -> sfc)
c     z:      altitude levels [in m]
c     rho:    density [in unit/m]         (unit arbitrary)
c     ip:     =1 : linear varying profile -> arithmetic mean for integration
c             =2 : approximately exponentially varying profile ->
c                  use average of arithmetic and geometric mean for integration
c             THIS HAS SMALLER ERROR IN INTEGRATION THAN EITHER 
c             ARITHMETIC (A) OR GEOMETRIC (G) MEAN
c             (if profile is varying exactly exponentially with z, then the integration error
c              is minimal for the combination: 2/3 G + 1/3 A)
c
c     output:
c     col [in unit/m]       
c
      implicit none
      integer ip, nlevel, i ,ierr
      real  col, dz , avg
      real, dimension(0:nlevel) :: z,rho
c     
      ierr = 0
      col = 0.
      profiles: do i=0,nlevel,1

      if (i.eq.0) cycle profiles

      if (z(i).le.z(i-1)) then
         write(*,*) 'error 1 in column: order of profiles wrong !'
         ierr=12
         return
      endif

      dz = z(i) - z(i-1)

      if (ip.eq.1) then
         avg = 0.5* (rho(i) + rho(i-1) )
      else if (ip.eq.2) then
         avg = 0.25* ( rho(i) + rho(i-1) + 2.*sqrt(rho(i-1)*rho(i)) )
      endif

      col = col + avg*dz

      enddo profiles

      return
      end

      subroutine fdpsat(T,P,  PSAT)
c
c     calculates water vapor saturation pressure PSAT
c     S. Cruz Pol, C. Ruf and S. Keihm, Radio Science 33 (5),1319 (1998)
c
      real :: T,P,Psat
        real :: F_w,xi1,xi2,xi3  
        real  :: p_standard = 1013.246  ! hPa
        real  :: TS = 373.14     ! K  (boiling point)
        real  :: T0 = 273.14     ! K  (freezing point)
c       real  :: eps = 1./1.607795  ! M(H2O)/M(air)
c       real  :: R_d = 287.05 ! J /(K*kg)  gas constant for dry air
c
      F_w = 1.0 + 1.E-4 * (5.92854 + 3.740346e-2 * P +
     1                1.971198E-4 * (T-T0) * (800-P)  +
     2                    6.045511E-6 * P * (T-T0)**2 ) ! deviation from ideal gas
c
      xi1 = -7.90298*(Ts/T - 1) + 5.02808*Log10(TS/T)
        xi2 = -1.3816E-7 * 10**(11.344*(1.-T/TS) -1 )
        xi3 = 8.1328E-3 * (10**(-3.49149*(TS/T - 1)) - 1)
c
      Psat = p_standard * 10**(xi1 + xi2 + xi3)

      return
        end  



      SUBROUTINE sort5(n,ra,rb,rc,rd,re,wksp,iwksp)
      INTEGER n,iwksp(n)
      REAL ra(n),rb(n),rc(n),rd(n),re(n),wksp(n)
CU    USES indexx
      INTEGER j
      call indexx(n,ra,iwksp)
      do 11 j=1,n
        wksp(j)=ra(j)
11    continue
      do 12 j=1,n
        ra(j)=wksp(iwksp(j))
12    continue
      do 13 j=1,n
        wksp(j)=rb(j)
13    continue
      do 14 j=1,n
        rb(j)=wksp(iwksp(j))
14    continue
      do 15 j=1,n
        wksp(j)=rc(j)
15    continue
      do 16 j=1,n
        rc(j)=wksp(iwksp(j))
16    continue
      do 17 j=1,n
        wksp(j)=rd(j)
17    continue
      do 18 j=1,n
        rd(j)=wksp(iwksp(j))
18    continue
      do 19 j=1,n
        wksp(j)=re(j)
19    continue
      do 20 j=1,n
        re(j)=wksp(iwksp(j))
20    continue

      return
      END



      SUBROUTINE YFDMONDAY(LYEAR,IDAYJL,IMON,IDAY)

        INTEGER*4 IDAYFX(12,0:1)
      DATA IDAYFX/1,32,60,91,121,152,182,213,244,274,305,335,           
     1        1,32,61,92,122,153,183,214,245,275,306,336/           
                                                  
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


      SUBROUTINE YFDSEC75(LYEAR,IDAYJL,ISECDY, ITIME)
      ITIME=31536000*(LYEAR-1975) + 86400*(IDAYJL-1) + ISECDY
     & + 86400*int((lyear - 1973)/4)
        RETURN
        END



