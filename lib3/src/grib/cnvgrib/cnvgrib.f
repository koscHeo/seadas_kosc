      program cnvgrib
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C  
C MAIN PROGRAM:  cnvgrib
C   PRGMMR: Gilbert        ORG: NP11        DATE: 2003-06-06
C
C ABSTRACT: This program converts every GRIB field in a file from
C   (1) GRIB1 to GRIB2   (2) GRIB2 to GRIB1  or (3) GRIB2 to GRIB2.
C
C PROGRAM HISTORY LOG:
C 2003-06-06  Gilbert
C
C USAGE:
C   INPUT FILES:
C     UNIT 10  - Input GRIB file
C
C   OUTPUT FILES: 
C     UNIT 50  - Output GRIB file
C
C   SUBPROGRAMS CALLED: (LIST ALL CALLED FROM ANYWHERE IN CODES)
C     UNIQUE:    - cnv12, cnv21, cnv22, usage
C     LIBRARY:
C       W3LIB    - errexit
C       BACIO    - baopenr, baopenw, baclose
C
C   EXIT STATES:
C     COND =   0 - SUCCESSFUL RUN
C          =   2 - Problem processing command line arguments
C          =   3 - Problem opening input GRIB file
C          =   4 - Problem opening output GRIB file
C          =   5 - Unknown conversion option
C
C REMARKS: LIST CAVEATS, OTHER HELPFUL HINTS OR INFORMATION
C
C ATTRIBUTES:
C   LANGUAGE: Fortran 90
C   MACHINE:  IBM SP
C
C$$$

      integer :: inver=0,outver=0,ipack=-1
      character(len=500) :: gfilein,gfileout,copt
      INTEGER(4) NARG,IARGC
      logical :: usemiss=.false., uvvect=.true.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET ARGUMENTS
      NARG=IARGC()
      IF(NARG.lt.3) THEN       ! may be a problem with args
         IF(NARG.eq.0) THEN
            !CALL ERRMSG('cnvgrib:  Incorrect usage')
            call usage(0)
            CALL ERREXIT(2)
         ELSE                !  look for -h "help" option
            do j=1,NARG
               call getarg(j,copt)
               if (copt.eq.'-h') then
                  call usage(1)
                  CALL ERREXIT(0)
               endif
            ENDDO
            call usage(0)
            CALL ERREXIT(2)
         ENDIF
      ELSE           
         j=1
         dowhile (j.le.NARG-2)        ! parse first narg-2 args
            call getarg(j,copt)
            j=j+1
            selectcase(copt)
             case('-g12')
                inver=1
                outver=2
             case('-g21')
                inver=2
                outver=1
             case('-g22')
                inver=2
                outver=2
             case('-p0')
                ipack=0
             case('-p2')
                ipack=2
             case('-p31')
                ipack=31
             case('-p32')
                ipack=32
             case('-p40')
                ipack=40
             case('-p41')
                ipack=41
             case('-p40000')       ! Obsolete 
                ipack=40000
             case('-p40010')       ! Obsolete 
                ipack=40010
             case('-m')
                usemiss=.true.
             case('-nv')
                uvvect=.false.
             case default
                call usage(0)
                CALL ERREXIT(2)
            end select
         ENDDO
         !
         !   get filenames from last two arguments
         !
         CALL GETARG(NARG-1,gfilein)
         CALL GETARG(NARG,gfileout)
         !
         !   If -p option specified, must be writing out grib2
         !
         if ( (ipack.ne.-1).and.(outver.eq.1) ) then
            CALL ERRMSG('cnvgrib: -pxx option ignored when using -g21')
         endif
         !
         !   Must have -g option
         !
         if ( (inver.eq.0).or.(outver.eq.0) ) then
            CALL ERRMSG('cnvgrib: must use one -gxx option')
            call usage(0)
            CALL ERREXIT(2)
         endif
         !
         !   If -m option specified, must be writing out grib2
         !   and using DRT 5.2 or 5.3
         !
         if ( (usemiss).and.(ipack.ne.2 .AND. ipack.ne.31 .AND.
     &                       ipack.ne.32) ) then
            CALL ERRMSG('cnvgrib: -m option ignored when not using '//
     &                 '-p2, -p31 or -p32.')
            usemiss=.false.
         endif
      ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Open input and output grib files
!
      IFL1=10
      IFL2=50
      NCGB=LEN_TRIM(gfilein)
      CALL BAOPENR(ifl1,gfilein(1:NCGB),IOS)
      if (IOS.NE.0) then
         call errmsg('cnvgrib: cannot open input GRIB file '//
     &               gfilein(1:NCGB))
         call errexit(3)
      endif
      NCGB=LEN_TRIM(gfileout)
      CALL BAOPENW(ifl2,gfileout(1:NCGB),IOS)
      if (IOS.NE.0) then
         call errmsg('cnvgrib: cannot open output GRIB file '//
     &               gfileout(1:NCGB))
         call errexit(4)
      endif
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  convert grib file
!
      if ((inver.eq.1).AND.(outver.eq.2)) then
         call cnv12(ifl1,ifl2,ipack,usemiss,uvvect)
      elseif ((inver.eq.2).AND.(outver.eq.1)) then
         call cnv21(ifl1,ifl2)
      elseif ((inver.eq.2).AND.(outver.eq.2)) then
         call cnv22(ifl1,ifl2,ipack,usemiss)
      else
         print *,' Unknown conversion option.'
         call errexit(5)
      endif
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  close grib files
!
      CALL BACLOSE(ifl1,IOS)
      CALL BACLOSE(ifl2,IOS)

      stop
      end

      subroutine usage(iopt)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    usage
C   PRGMMR: Gilbert     ORG: W/NP11    DATE: 2003-06-06
C
C ABSTRACT: This routine prints out the command "usage" 
C   or a brief description of the command line options.
C
C PROGRAM HISTORY LOG:
C 2003-06-06  Gilbert
C 2007-04-25  Vuong   -  Changed the cnvgrib_ver
C
C USAGE:    CALL usage(iopt)
C   INPUT ARGUMENT LIST:
C     iopt   - ouput option:
C                   1 = print description of arguments
C                   otherwise, print command usage summary
C
C REMARKS: None
C
C ATTRIBUTES:
C   LANGUAGE: Fortran 90
C   MACHINE:  IBM SP
C
C$$$
         character(len=15) :: cnvgrib_ver="cnvgrib-1.1.6"
         integer,intent(in) :: iopt 
         call errmsg('Usage: cnvgrib [-h] {-g12|-g21|-g22} [-m] [-nv]')
         call errmsg('               [{-p0|-p2|-p31|-p32|-p40'//
     &               '|-p41|-p40000|-p40010}] ingribfile outgribfile')
         if ( iopt.eq.1 ) then
            call errmsg('cnvgrib:  version '//cnvgrib_ver)
            call errmsg('Must use one of the following options:')
            call errmsg('   -g12     converts GRIB1 to GRIB2')
            call errmsg('   -g21     converts GRIB2 to GRIB1')
            call errmsg('   -g22     converts GRIB2 to GRIB2 '//
     &                  ' (used to change packing option)')
            call errmsg('Optional packing options: (for use with '//
     &                   ' -g12 and -g22 only)')
            call errmsg('   -p0      simple packing')
            call errmsg('   -p2      complex packing')
            call errmsg('   -p31     complex pack with 1st order diffs')
            call errmsg('   -p32     complex pack with 2nd order diffs')
            call errmsg('   -p40     JPEG2000 encoding')
            call errmsg('   -p41     PNG encoding')
            call errmsg('   -p40000  JPEG2000 encoding (Obsolete)')
            call errmsg('   -p40010  PNG encoding (Obsolete)')
            call errmsg('Other Optional options: ')
         call errmsg('   -m       Use missing values instead of bitmap')
            call errmsg('            (valid with -p2, -p31 or -p32 '//
     &                               'options only)')
          call errmsg('   -nv      Do not combine U, V wind components')
         endif
      return
      end
