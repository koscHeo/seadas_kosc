      PROGRAM MK_AQUARIUS_ANCILLARY_DATA
      USE HDF5                  ! This module contains all necessary modules 
      IMPLICIT NONE

!     Generates the "y" 6 hour ancillary files
      CHARACTER(100) FILENAME

!      REAL(4), TARGET :: SST(360,181),UU(360,181),VV(360,181),surtemp(360,181)
      REAL(4), TARGET :: UU(360,181),VV(360,181),surtemp(360,181)
      REAL(4), TARGET :: surp(360,181),cwat(360,181),soilw(360,181)
      REAL(4), TARGET :: bgrtemp(360,181)
      REAL(4), TARGET :: swe(720,361), swh(288,157)
      REAL(4), TARGET ::  TRAN(360,181,91),TBUP(360,181,91),TBDW(360,181,91)
      real(4), target :: fice(720,360), solarflux(1)
      real(4), target ::  surtemp_out(360*4,180*4), sssref(360*4,180*4)

      REAL(4) XLAT,XLON, sst1(360*4,180*4), sst2(360*4,180*4)
      real(4) solarflux1, solarflux2
      real(4) sss1(360*4,180*4), sss2(360*4,180*4)

      INTEGER(4) LYEAR,IMON,IDAYMO,IDAYJL,IHOUR,ILAT,ILON
      INTEGER(4) LYEAR_qmet,IDAYJL_qmet,IHOUR_qmet
      INTEGER(4) lsfyear1, lsfmon1, lsfdaymo1
      INTEGER(4) lsfyear2, lsfmon2, lsfdaymo2
      INTEGER(4) lsstyear1, lsstdayjl1
      INTEGER(4) lsstyear2, lsstdayjl2
      INTEGER(4) lsssyear1, lsssdayjl1
      INTEGER(4) lsssyear2, lsssdayjl2
      INTEGER(4) liceyear1, licedayjl1
      INTEGER(4) liceyear2, licedayjl2
      INTEGER(4) imon1, imon2, idaymo1, idaymo2
      real*8 julsec, julsec1, julsec2
      
      CHARACTER(128) OISSTFILE1                            
      CHARACTER(128) OISSTFILE2                            
      CHARACTER(128) atmos_file
      CHARACTER(128) met_file, swh_file
      CHARACTER(128) seaice_file1, seaice_file2
      CHARACTER(128) sssref_file1, sssref_file2
      CHARACTER(128) solarflux_date1, solarflux_date2, outfile_date

      INTEGER(HID_T) :: file_id   
      INTEGER(HID_T), DIMENSION(16) :: dset_id
      INTEGER :: rank
      INTEGER(HSIZE_T), DIMENSION(3) :: dims, dims_gfs, dims_sst, dims_swh
      INTEGER(HSIZE_T), DIMENSION(3) :: offset=(/0,0,0/)
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace
      real(4) a1, a2

      real(4), POINTER :: pt

      integer :: slash_pos, nargs, len, lenstr

      CHARACTER solarflux_file*255

      CHARACTER*80, DIMENSION(9) ::  attr_data ! Attribute data
      INTEGER(HID_T) :: attr_id ! Attribute identifier 
      INTEGER(HID_T) :: aspace_id ! Attribute Dataspace identifier 
      INTEGER(HID_T) :: atype_id ! Attribute Dataspace identifier 
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/8/) ! Attribute dimension
      INTEGER     ::   arank = 1 ! Attribure rank
      INTEGER(SIZE_T) :: attrlen = 80 ! Length of the attribute string
      INTEGER     ::   error    ! Error flag

      INTERFACE
         SUBROUTINE READ_SSSREF(sssref_file, sss)
         USE HDF5               ! This module contains all necessary modules
         IMPLICIT NONE

         CHARACTER(100) sssref_file
         REAL(4), TARGET ::  SSS(1440,720)
         END SUBROUTINE READ_SSSREF

         SUBROUTINE READ_SWH_NCEP(swh_file, swh)
         USE HDF5               ! This module contains all necessary modules 
         IMPLICIT NONE

         CHARACTER(100) swh_file
         REAL(4), TARGET ::  swh(288,157)
         END SUBROUTINE READ_SWH_NCEP
      END INTERFACE

      include "aquarius.fin"

      write(*,*) 'mk_aquarius_ancillary_data 1.08 ('//__DATE__//' '//__TIME__//')'
      write(*,*) ''
      nargs = IArgC()
      IF (nargs.eq.0) THEN
         write(*,*) 'mk_aquarius_ancillary beg_sstfile end_sstfile qatm_file qmet_file swh_file'
         write(*,*) '                      beg_seaicefile end_seaicefile'
         write(*,*) '                      beg_sssreffile end_sssreffile'
         write(*,*) '                      solarflux_date1 solarflux_date2'
         write(*,*) ' '
         write(*,*) 'where beg_sstfile is the beginning SST ancillary file'
         write(*,*) '      end_sstfile is the ending SST ancillary file'
         write(*,*) '      qatm_file is the QATM ancillary file'
         write(*,*) '      qmet_file is the QMET ancillary file'
         write(*,*) '      swh_file is the SWH ancillary file'
         write(*,*) '      beg_seaicefile is the beginning SEAICE ancillary file'
         write(*,*) '      end_seaicefile is the ending SEAICE ancillary file'
         write(*,*) '      beg_sssreffile is the beginning SSS reference ancillary file'
         write(*,*) '      end_sssreffile is the ending SSS reference ancillary file'
         write(*,*) '      solarflux_date1 is the beginning date for the solar flux'
         write(*,*) '      solarflux_date2 is the ending date for the solar flux'
         write(*,*) ' '
         stop
      ENDIF

      call getarg(1, oisstfile1)
      call getarg(2, oisstfile2)
      slash_pos = INDEX(OISSTFILE1,'/',BACK=.TRUE.)
      attr_data(1) = OISSTFILE1(slash_pos+1:)
      slash_pos = INDEX(OISSTFILE2,'/',BACK=.TRUE.)
      attr_data(2) = OISSTFILE2(slash_pos+1:)

      call getarg(3, atmos_file)
      call getarg(4, met_file)
      slash_pos = INDEX(atmos_file,'/',BACK=.TRUE.)
      attr_data(3) = atmos_file(slash_pos+1:)
      slash_pos = INDEX(met_file,'/',BACK=.TRUE.)
      attr_data(4) = met_file(slash_pos+1:)

      call getarg(5, swh_file)
      slash_pos = INDEX(swh_file,'/',BACK=.TRUE.)
      attr_data(5) = swh_file(slash_pos+1:)

      call getarg(6, seaice_file1)
      call getarg(7, seaice_file2)
      slash_pos = INDEX(seaice_file1,'/',BACK=.TRUE.)
      attr_data(6) = seaice_file1(slash_pos+1:)
      slash_pos = INDEX(seaice_file2,'/',BACK=.TRUE.)
      attr_data(7) = seaice_file2(slash_pos+1:)

      call getarg(8, sssref_file1)
      call getarg(9, sssref_file2)
      slash_pos = INDEX(sssref_file1,'/',BACK=.TRUE.)
      attr_data(8) = sssref_file1(slash_pos+1:)
      slash_pos = INDEX(sssref_file2,'/',BACK=.TRUE.)
      attr_data(9) = sssref_file2(slash_pos+1:)

      call getarg(10, solarflux_date1)
      call getarg(11, solarflux_date2)

      call getarg(12, outfile_date)


c     Get basename of metfile and extract hour 
      if (outfile_date(1:1).eq.' ') then
         slash_pos = INDEX(met_file,'/',BACK=.TRUE.)

         READ(UNIT=met_file(slash_pos+2:slash_pos+5), FMT='(I4)') lyear
!     write(*,*) lyear

         READ(UNIT=met_file(slash_pos+6:slash_pos+8), FMT='(I4)') idayjl
!     write(*,*) idayjl

         READ(UNIT=met_file(slash_pos+9:slash_pos+10), FMT='(I2)') ihour
!     write(*,*) ihour
      else
         READ(UNIT=outfile_date(1:4), FMT='(I4)') lyear
         READ(UNIT=outfile_date(5:7), FMT='(I3)') idayjl
         READ(UNIT=outfile_date(8:9), FMT='(I2)') ihour

         READ(UNIT=met_file(slash_pos+2:slash_pos+5), FMT='(I4)') lyear_qmet
         READ(UNIT=met_file(slash_pos+6:slash_pos+8), FMT='(I4)') idayjl_qmet
         READ(UNIT=met_file(slash_pos+9:slash_pos+10), FMT='(I2)') ihour_qmet

         if ((lyear.ne.lyear_qmet).or.(idayjl.ne.idayjl_qmet).or.(ihour.ne.ihour_qmet)) then
            write(*,200) 'Interpolation time ',lyear,idayjl,ihour,
     1           ' not equal to QMET time: ',lyear_qmet,idayjl_qmet,ihour_qmet
 200     format(a,i4,i3.3,i2.2,a,i4,i3.3,i2.2)
         endif
      endif

      call FIND_MONTH_DAY(LYEAR, IDAYJL, IMON, IDAYMO)
!      write(*,*) imon, idaymo

c     open oisst1
      CALL READ_OISST(OISSTFILE1, SST1)

c     open oisst2
      CALL READ_OISST(OISSTFILE2, SST2)

      CALL READ_SURFACE_TEMP_NCEP(met_file, surtemp)
      CALL READ_SURFACE_PRESSURE_NCEP(met_file, surp)
      CALL READ_CWAT_NCEP(met_file, cwat)
      CALL READ_SOILW_NCEP(met_file, soilw)
      CALL READ_bgrtemp_NCEP(met_file, bgrtemp)
      CALL READ_SWE_NCEP_GFS(met_file, swe)

      CALL READ_SWH_NCEP(swh_file, swh)

c     sst processing
!      READ(UNIT=oisstfile1(2:5), FMT='(I4)') lsstyear1
!      READ(UNIT=oisstfile2(2:5), FMT='(I4)') lsstyear2
!      READ(UNIT=oisstfile1(6:8), FMT='(I4)') lsstdayjl1
!      READ(UNIT=oisstfile2(6:8), FMT='(I4)') lsstdayjl2

      READ(UNIT=attr_data(1)(2:5), FMT='(I4)') lsstyear1
      READ(UNIT=attr_data(2)(2:5), FMT='(I4)') lsstyear2
      READ(UNIT=attr_data(1)(6:8), FMT='(I4)') lsstdayjl1
      READ(UNIT=attr_data(2)(6:8), FMT='(I4)') lsstdayjl2

      call FIND_MONTH_DAY(lsstyear1, lsstdayjl1, IMON1, IDAYMO1)
      call FIND_MONTH_DAY(lsstyear2, lsstdayjl2, IMON2, IDAYMO2)

      ! For Aquarius, SST files start at 12 noon
      call  ymdhms2jul(lsstyear1,imon1,idaymo1,12,0,0,julsec1)
      call  ymdhms2jul(lsstyear2,imon2,idaymo2,12,0,0,julsec2)
      call  ymdhms2jul(lyear,imon,idaymo,ihour,0,0,julsec)

      a2 = 0
      if ( julsec1.ne.julsec2) a2 = (julsec - julsec1) / (julsec2 - julsec1)
      a1 = 1 - a2
      if ( (a1.lt.0).or.(a2.lt.0)) then
         write(*,100) 'Time ',lyear,imon,idaymo,ihour,
     1   ' out of SST interpolation range: ',
     2   lsstyear1,imon1,idaymo1,12,'-',lsstyear2,imon2,idaymo2,12
 100     format(a,i4,i2.2,i2.2,i2.2,a,i4,i2.2,i2.2,i2.2,a,i4,i2.2,i2.2,i2.2)
!         write(*,*) 'SST files start at 12 noon'
!         call exit(1)
      endif

       
!     Note: Time interpolates between sst{1|2}(360*4,180*4) and surtemp_out (360*4,180*4)
       DO ILAT=1,180*4
          XLAT=89.75 - 0.25 * (ilat-1)
          DO ILON=1,360*4
             XLON=0.125 + 0.25 * (ilon-1)
             surtemp_out(ilon,ilat) = surtemp(ilon/4,ilat*181/(180*4))
             if (sst1(ilon,721-ilat).ne.0) then
                surtemp_out(ilon,ilat) = a1*sst1(ilon,721-ilat) + a2*sst2(ilon,721-ilat) + 273.15
             endif
          ENDDO
       ENDDO


c     SSS reference data processing
      READ(UNIT=attr_data(8)(2:5), FMT='(I4)') lsssyear1
      READ(UNIT=attr_data(9)(2:5), FMT='(I4)') lsssyear2
      READ(UNIT=attr_data(8)(6:8), FMT='(I4)') lsssdayjl1
      READ(UNIT=attr_data(9)(6:8), FMT='(I4)') lsssdayjl2
      call FIND_MONTH_DAY(lsssyear1, lsssdayjl1, IMON1, IDAYMO1)
      call FIND_MONTH_DAY(lsssyear2, lsssdayjl2, IMON2, IDAYMO2)

      ! For Aquarius, SSS files start at 12 noon
      call  ymdhms2jul(lsstyear1,imon1,idaymo1,12,0,0,julsec1)
      call  ymdhms2jul(lsstyear2,imon2,idaymo2,12,0,0,julsec2)
      call  ymdhms2jul(lyear,imon,idaymo,ihour,0,0,julsec)

      a2 = 0
      if ( julsec1.ne.julsec2) a2 = (julsec - julsec1) / (julsec2 - julsec1)
      a1 = 1 - a2
      if ( (a1.lt.0).or.(a2.lt.0)) then
         write(*,100) 'Time ',lyear,imon,idaymo,ihour,
     1   ' out of SSS interpolation range: ',
     2   lsssyear1,imon1,idaymo1,12,'-',lsssyear2,imon2,idaymo2,12
!         write(*,*) 'SSS files start at 12 noon'
!         call exit(1)
      endif

      CALL READ_SSSREF(sssref_FILE1, sss1)
      CALL READ_SSSREF(sssref_FILE2, sss2)

       DO ILAT=1,180*4
          DO ILON=1,360*4
!    Don't relect about equator 11/13/10 JMG
!             sssref(ilon,ilat) = a1*sss1(ilon,721-ilat)! + a2*sst2(ilon,721-ilat)
             sssref(ilon,ilat) = a1*sss1(ilon,ilat) + a2*sss2(ilon,ilat)
             if (sssref(ilon,ilat).gt.1e30) sssref(ilon,ilat) = -999
          ENDDO
       ENDDO


c     seaice processing
      READ(UNIT=attr_data(6)(2:5), FMT='(I4)') liceyear1
      READ(UNIT=attr_data(7)(2:5), FMT='(I4)') liceyear2
      READ(UNIT=attr_data(6)(6:8), FMT='(I4)') licedayjl1
      READ(UNIT=attr_data(7)(6:8), FMT='(I4)') licedayjl2
      call FIND_MONTH_DAY(liceyear1, licedayjl1, IMON1, IDAYMO1)
      call FIND_MONTH_DAY(liceyear2, licedayjl2, IMON2, IDAYMO2)

      call  ymdhms2jul(liceyear1,imon1,idaymo1,0,0,0,julsec1)
      call  ymdhms2jul(liceyear2,imon2,idaymo2,0,0,0,julsec2)
      call  ymdhms2jul(lyear,imon,idaymo,0,0,0,julsec)

      a2 = 0
      if ( julsec1.ne.julsec2) a2 = (julsec - julsec1) / (julsec2 - julsec1)
      a1 = 1 - a2
      if ( (a1.lt.0).or.(a2.lt.0)) then
         write(*,100) 'Time ',lyear,imon,idaymo,ihour,
     1   ' out of Sea Ice interpolation range: ',
     2   liceyear1,imon1,idaymo1,0,'-',liceyear2,imon2,idaymo2,0
      endif


!      write(*,*) 'fd_fice commented out'
      call fd_fice(seaice_file1, seaice_file2, a1, a2, fice)


c     Solar flux processing
      call getenv('OCVARROOT', solarflux_file)
      len = lenstr(solarflux_file)
      solarflux_file = solarflux_file(1:len)//'/aquarius/solar_flux_noon.txt'
!      write(*,*) 'Reading ', solarflux_file

      call read_solarflux(solarflux_file, solarflux_date1, solarflux1)
      call read_solarflux(solarflux_file, solarflux_date2, solarflux2)

      READ(UNIT=solarflux_date1(1:4), FMT='(I4)') lsfyear1
      READ(UNIT=solarflux_date2(1:4), FMT='(I4)') lsfyear2
      READ(UNIT=solarflux_date1(5:6), FMT='(I4)') lsfmon1
      READ(UNIT=solarflux_date2(5:6), FMT='(I4)') lsfmon2
      READ(UNIT=solarflux_date1(7:8), FMT='(I4)') lsfdaymo1
      READ(UNIT=solarflux_date2(7:8), FMT='(I4)') lsfdaymo2

      call  ymdhms2jul(lsfyear1,lsfmon1,lsfdaymo1,0,0,0,julsec1)
      call  ymdhms2jul(lsfyear2,lsfmon2,lsfdaymo2,0,0,0,julsec2)
      call  ymdhms2jul(lyear,imon,idaymo,ihour,0,0,julsec)

      a2 = 0
      if ( julsec1.ne.julsec2) a2 = (julsec - julsec1) / (julsec2 - julsec1)
      a1 = 1 - a2
      if ( (a1.lt.0).or.(a2.lt.0)) then
         write(*,100) 'Time ',lyear,imon,idaymo,ihour,
     1   ' out of Solar Flux interpolation range: ',
     2   lsfyear1,lsfmon1,lsfdaymo1,0,'-',lsfyear2,lsfmon2,lsfdaymo2,0
      endif

      solarflux = a1 * solarflux1 + a2 * solarflux2

  
      CALL READ_WINDS_NCEP_10M(met_file, UU,VV)

      CALL READ_ATMOS(atmos_file, TRAN,TBUP,TBDW)

      WRITE(FILENAME,9001) LYEAR,IMON,IDAYMO,IHOUR
 9001 FORMAT('y',I4.4,I2.2,I2.2,I2.2,'.h5')
      call create_hdf5(filename, file_id)

      rank = 2
      dims(1) = 360
      dims(2) = 181

      dims_gfs(1) = 720
      dims_gfs(2) = 361

      dims_sst(1) = 360*4
      dims_sst(2) = 180*4

      dims_swh(1) = 288
      dims_swh(2) = 157

      call create_hdf5_real(file_id, 'UU Wind', rank, dims, dset_id(2))
      call create_hdf5_real(file_id, 'VV Wind', rank, dims, dset_id(3))
      call create_hdf5_real(file_id, 'Surface Temp', rank, dims_sst, dset_id(7))
      call create_hdf5_real(file_id, 'Surface Pressure', rank, dims, dset_id(8))
      call create_hdf5_real(file_id, 'Cloud Water', rank, dims, dset_id(9))
      call create_hdf5_real(file_id, 'Soil Moisture', rank, dims, dset_id(10))
      call create_hdf5_real(file_id, 'Sub-Surface Temp', rank, dims, dset_id(15))
      call create_hdf5_real(file_id, 'Significant Wave Height', rank, dims_swh, dset_id(16))

      if ( swe(1,1) .ne. -9999.) then
         call create_hdf5_real(file_id, 'Snow Water', rank, dims_gfs, dset_id(11))
      endif

      call set_space_hdf5(dset_id(2), rank, dims, filespace, dataspace)
      pt => UU(1,1)
      call write_hdf5_real(dset_id(2), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(2), dataspace, filespace)

      call set_space_hdf5(dset_id(3), rank, dims, filespace, dataspace)
      pt => VV(1,1)
      call write_hdf5_real(dset_id(3), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(3), dataspace, filespace)

      call set_space_hdf5(dset_id(7), rank, dims_sst, filespace, dataspace)
      pt => surtemp_out(1,1)
      call write_hdf5_real(dset_id(7), filespace, dataspace, offset, dims_sst, pt)
      call close_hdf5_ds(dset_id(7), dataspace, filespace)

      call set_space_hdf5(dset_id(8), rank, dims, filespace, dataspace)
      pt => surp(1,1)
      call write_hdf5_real(dset_id(8), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(8), dataspace, filespace)

      call set_space_hdf5(dset_id(9), rank, dims, filespace, dataspace)
      pt => cwat(1,1)
      call write_hdf5_real(dset_id(9), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(9), dataspace, filespace)

      call set_space_hdf5(dset_id(10), rank, dims, filespace, dataspace)
      pt => soilw(1,1)
      call write_hdf5_real(dset_id(10), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(10), dataspace, filespace)

      call set_space_hdf5(dset_id(15), rank, dims, filespace, dataspace)
      pt => bgrtemp(1,1)
      call write_hdf5_real(dset_id(15), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(15), dataspace, filespace)

      call set_space_hdf5(dset_id(16), rank, dims_swh, filespace, dataspace)
      pt => swh(1,1)
      call write_hdf5_real(dset_id(16), filespace, dataspace, offset, dims_swh, pt)
      call close_hdf5_ds(dset_id(16), dataspace, filespace)

      if ( swe(1,1) .ne. -9999.) then
         call set_space_hdf5(dset_id(11), rank, dims_gfs, filespace, dataspace)
         pt => swe(1,1)
         call write_hdf5_real(dset_id(11), filespace, dataspace, offset, dims_gfs, pt)
         call close_hdf5_ds(dset_id(11), dataspace, filespace)
      endif

      dims(1) = 360*2
      dims(2) = 180*2
      call create_hdf5_real(file_id, 'FracIce', rank, dims, dset_id(12))
      call set_space_hdf5(dset_id(12), rank, dims, filespace, dataspace)
      pt => fice(1,1)
      call write_hdf5_real(dset_id(12), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(12), dataspace, filespace)


      dims(1) = 360*4
      dims(2) = 180*4
      call create_hdf5_real(file_id, 'SSS_Ref', rank, dims, dset_id(14))
      call set_space_hdf5(dset_id(14), rank, dims, filespace, dataspace)
      pt => sssref(1,1)
      call write_hdf5_real(dset_id(14), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(14), dataspace, filespace)


      rank = 1
      dims(1) = 1
      call create_hdf5_real(file_id, 'Solar Flux', rank, dims, dset_id(13))
      call set_space_hdf5(dset_id(13), rank, dims, filespace, dataspace)
      pt => solarflux(1)
      call write_hdf5_real(dset_id(13), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(13), dataspace, filespace)


      rank = 3
      dims(1) = 360
      dims(2) = 181
      dims(3) = 91

      call create_hdf5_real(file_id, 'Transmittance', rank, dims, dset_id(4))
      call create_hdf5_real(file_id, 'Upwelling_TB', rank, dims, dset_id(5))
      call create_hdf5_real(file_id, 'Downwelling_TB', rank, dims, dset_id(6))

      call set_space_hdf5(dset_id(4), rank, dims, filespace, dataspace)
      pt => TRAN(1,1,1)
      call write_hdf5_real(dset_id(4), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(4), dataspace, filespace)

      call set_space_hdf5(dset_id(5), rank, dims, filespace, dataspace)
      pt => TBUP(1,1,1)
      call write_hdf5_real(dset_id(5), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(5), dataspace, filespace)

      call set_space_hdf5(dset_id(6), rank, dims, filespace, dataspace)
      pt => TBDW(1,1,1)
      call write_hdf5_real(dset_id(6), filespace, dataspace, offset, dims, pt)
      call close_hdf5_ds(dset_id(6), dataspace, filespace)


      ! Write Input Ancillary File array attribute
      CALL h5screate_simple_f(arank, adims, aspace_id, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
      CALL h5tset_size_f(atype_id, attrlen, error)
      
      CALL h5acreate_f(file_id, 'Input Ancillary Files', 
     1                 atype_id, aspace_id, attr_id, error)
      CALL h5awrite_f(attr_id, atype_id, attr_data, adims, error)

      call close_hdf5_df(file_id)
      
      STOP 'NORM END'
      END


      SUBROUTINE READ_WINDS_NCEP_10M(metfile, UU,VV)
      IMPLICIT NONE
      REAL(4) UU(360,181),VV(360,181)

      INTEGER(4) irow
      integer, parameter :: DFACC_READ = 1

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 
      real flt(360)
      CHARACTER*128 metfile

      sd_id_anc = sfstart( metfile, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 360
      count(2) = 181

      sds_index = sfn2index(sd_id_anc, 'u_wind')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, uu)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = uu(:,irow)
         uu(:,irow) = uu(:,181-irow+1)
         uu(:,181-irow+1) = flt
      enddo

      sds_index = sfn2index(sd_id_anc, 'v_wind')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, vv)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = vv(:,irow)
         vv(:,irow) = vv(:,181-irow+1)
         vv(:,181-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 

      RETURN
      END



      SUBROUTINE READ_SURFACE_TEMP_NCEP(metfile, surtemp)
      IMPLICIT NONE
      REAL(4) surtemp(360,181)

      INTEGER(4) irow
      integer, parameter :: DFACC_READ = 1

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 
      real flt(360)
      CHARACTER*128 metfile

      sd_id_anc = sfstart( metfile, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 360
      count(2) = 181

      sds_index = sfn2index(sd_id_anc, 'tmp_sfc')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, surtemp)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = surtemp(:,irow)
         surtemp(:,irow) = surtemp(:,181-irow+1)
         surtemp(:,181-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 

      RETURN
      END


      SUBROUTINE READ_SURFACE_PRESSURE_NCEP(metfile, surp)
      IMPLICIT NONE
      REAL(4) surp(360,181)

      INTEGER(4) irow
      integer, parameter :: DFACC_READ = 1

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 
      real flt(360)
      CHARACTER*128 metfile

      sd_id_anc = sfstart( metfile, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 360
      count(2) = 181

      sds_index = sfn2index(sd_id_anc, 'press')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, surp)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = surp(:,irow)
         surp(:,irow) = surp(:,181-irow+1)
         surp(:,181-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 

      RETURN
      END


      SUBROUTINE READ_CWAT_NCEP(metfile, cwat)
      IMPLICIT NONE
      REAL(4) cwat(360,181)

      INTEGER(4) irow
      integer, parameter :: DFACC_READ = 1

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 
      real flt(360)
      CHARACTER*128 metfile

      sd_id_anc = sfstart( metfile, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 360
      count(2) = 181

      sds_index = sfn2index(sd_id_anc, 'c_water')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, cwat)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = cwat(:,irow)
         cwat(:,irow) = cwat(:,181-irow+1)
         cwat(:,181-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 

      RETURN
      END


      SUBROUTINE READ_SOILW_NCEP(metfile, soilw)
      IMPLICIT NONE
      REAL(4) soilw(360,181)

      INTEGER(4) irow
      integer, parameter :: DFACC_READ = 1

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 
      real flt(360)
      CHARACTER*128 metfile

      sd_id_anc = sfstart( metfile, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 360
      count(2) = 181

      sds_index = sfn2index(sd_id_anc, 'soil_w')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, soilw)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = soilw(:,irow)
         soilw(:,irow) = soilw(:,181-irow+1)
         soilw(:,181-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 

      RETURN
      END

      SUBROUTINE READ_bgrtemp_NCEP(metfile, bgrtemp)
      IMPLICIT NONE
      REAL(4) bgrtemp(360,181)

      INTEGER(4) irow
      integer, parameter :: DFACC_READ = 1

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 
      real flt(360)
      CHARACTER*128 metfile

      sd_id_anc = sfstart( metfile, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 360
      count(2) = 181

      sds_index = sfn2index(sd_id_anc, 'tmp_bgr')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, bgrtemp)
      retn = sfendacc(sds_id)
      do irow=1,90
         flt = bgrtemp(:,irow)
         bgrtemp(:,irow) = bgrtemp(:,181-irow+1)
         bgrtemp(:,181-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 

      RETURN
      END


      SUBROUTINE READ_SWE_NCEP_GFS(metfile, swe)
      IMPLICIT NONE
      REAL(4) swe(720,361)

      INTEGER(4) irow
      integer, parameter :: DFACC_READ = 1

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 
      real flt(720)
      CHARACTER*128 metfile

      sd_id_anc = sfstart( metfile, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 720
      count(2) = 361

      sds_index = sfn2index(sd_id_anc, 'snow_gfs')
      if ( sds_index.eq.-1) then
         write(*,*) '"snow_gfs" field not found in QMET ancillary file'
         write(*,*) ''

         retn = sfend(sd_id_anc) 
         swe(1,1) = -9999.
         return
      endif
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, swe)
      retn = sfendacc(sds_id)
      do irow=1,180
         flt = swe(:,irow)
         swe(:,irow) = swe(:,361-irow+1)
         swe(:,361-irow+1) = flt
      enddo

      retn = sfend(sd_id_anc) 

      RETURN
      END



      SUBROUTINE READ_OISST(oisstfile, sst)
      IMPLICIT NONE

      real(4) SST(360*4,180*4)
      integer, parameter :: DFACC_READ = 1

      integer(2) isst(360*4,180*4)
      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sffattr, sfrnatt, sfginfo
      integer sfendacc, sfend
      integer , DIMENSION(4) :: start=(/0,0,0,0/), stride=(/1,1,1,1/)
      integer retn 
      CHARACTER*128 oisstfile
      CHARACTER*64 sds_name
      integer*4 attr_id, ilat, ilon, rank, dimsizes(8), data_type, num_attrs
      real slope, intercept

      sd_id_anc = sfstart( oisstfile, DFACC_READ)
      if ( sd_id_anc.eq.-1) then
         write(*,*) oisstfile//' cannot be opened.'
         call exit(1)
      endif

      sds_index = sfn2index(sd_id_anc, 'sst')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfginfo(sds_id, sds_name, rank, dimsizes, data_type, num_attrs)

      retn = sfrdata(sds_id, start, stride, dimsizes, isst)

      attr_id = sffattr(sds_id, 'scale_factor') 
      retn = sfrnatt(sds_id, attr_id, slope)
      attr_id = sffattr(sds_id, 'add_offset') 
      retn = sfrnatt(sds_id, attr_id, intercept)

      retn = sfendacc(sds_id)
      retn = sfend(sd_id_anc) 

      DO ILAT=1,180*4
         DO ILON=1,360*4
            if ( isst(ilon,ilat).eq.-999) then
               sst(ilon,ilat) = 0
            else
               sst(ilon,ilat) = isst(ilon,ilat)*slope + intercept
            endif
         enddo
      enddo

      RETURN
      END



      SUBROUTINE READ_SSSREF(sssref_file, sss)
      USE HDF5                  ! This module contains all necessary modules 
      IMPLICIT NONE

      CHARACTER(100) sssref_file
      REAL(4), TARGET ::  SSS(1440,720)
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      REAL(4), POINTER :: pt

      include "aquarius.fin"

      data_dims(1) = 360*4
      data_dims(2) = 180*4

      pt => sss(1,1)
      call dump_hdf5(sssref_file, 'salinity', data_dims, pt) 

      RETURN
      END



      SUBROUTINE READ_ATMOS(atmos_file, TRAN,TBUP,TBDW)
      USE HDF5                  ! This module contains all necessary modules 
      IMPLICIT NONE

      CHARACTER(100) atmos_file
      REAL(4), TARGET ::  TRAN1(360,181,91),TBUP1(360,181,91),TBDW1(360,181,91)
      REAL(4) TRAN( 360,181,91),TBUP( 360,181,91),TBDW( 360,181,91)

      INTEGER(4) ILAT,JLAT
 
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      REAL(4), POINTER :: pt

      include "aquarius.fin"

      data_dims(1) = 360
      data_dims(2) = 181
      data_dims(3) = 91

      pt => tran1(1,1,1)
      call dump_hdf5(atmos_file, 'Transmittance', data_dims, pt) 

      pt => tbup1(1,1,1)
      call dump_hdf5(atmos_file, 'Upwelling_TB', data_dims, pt) 

      pt => tbdw1(1,1,1)
      call dump_hdf5(atmos_file, 'Downwelling_TB', data_dims, pt) 

C     NEED TO FLIP THOMAS FILES, WHICH GO FROM -90S TO -90N TO NCEP CONVENTION OF 90N TO 90S

      DO ILAT=1,181
         JLAT=182-ILAT
         TRAN(:,JLAT,:)=TRAN1(:,ILAT,:)
         TBUP(:,JLAT,:)=TBUP1(:,ILAT,:)
         TBDW(:,JLAT,:)=TBDW1(:,ILAT,:)
      ENDDO


      RETURN

      END


      SUBROUTINE READ_SWH_NCEP(swh_file, swh)
      USE HDF5                  ! This module contains all necessary modules 
      IMPLICIT NONE

      CHARACTER(100) swh_file
      REAL(4), TARGET ::  swh(288,157)

      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      REAL(4), POINTER :: pt

      include "aquarius.fin"

      data_dims(1) = 288
      data_dims(2) = 157

      pt => swh(1,1)
      call dump_hdf5(swh_file, 'significant wave height', data_dims, pt) 

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
 10   CONTINUE
      IMON=12
 20   CONTINUE
      
      IDAY=1+IDAYJL-IDAYFX(IMON,ILEAP)
      RETURN
      END



C     11/17/00 VERSION CHANGE ON JULY 2, 2002 TO MAKE PROGRAM MORE ROBUST TO
C     CASES WHERE THE REYNOLDS FILE IS BEING UPDATED.   The semiphore file is used
C     as well as OPENBIG_TRY.  Also the reset function has been added


      SUBROUTINE FDREYNOLDV2(sst1, sst2, a1, a2, XLAT, XLON, SST)
      
      real SST1(360*4,180*4)
      real SST2(360*4,180*4)
      real(4) a1, a2
      real xlat, xlon, sst

C     DO TIME,LAT,LON INTERPOLATION 

      BRIEF=XLAT+89.5
      J1=1+BRIEF                                                        
      J2=J1+1                                                           
      B1=J1-BRIEF                                                       
      B2=1-B1
      IF(J1.EQ.  0) J1=  1
      IF(J2.EQ.181) J2=180
C
      BRIEF=XLON-0.5 
      K1=1+BRIEF                                                       
      K2=K1+1                                                           
      C1=K1-BRIEF                                                       
      C2=1-C1                                                          
      IF(K1.EQ.  0) K1=360
      IF(K2.EQ.361) K2=  1 
C                                                                       
      k1 = 4 * k1
      k2 = 4 * k2
      j1 = 4 * j1
      j2 = 4 * j2

       SST=             
     1 A1*B1*(C1*SST1(K1,J1)+C2*SST1(K2,J1))+                 
     2 A1*B2*(C1*SST1(K1,J2)+C2*SST1(K2,J2))+                 
     3 A2*B1*(C1*SST2(K1,J1)+C2*SST2(K2,J1))+                 
     4 A2*B2*(C1*SST2(K1,J2)+C2*SST2(K2,J2)) 
                      
      RETURN                                                            
      END                                                               


      subroutine fd_fice(seaice_file1, seaice_file2, a1, a2, fice)
      implicit none
 
      character(100), intent(in) ::  seaice_file1, seaice_file2
      real(4), intent(in) :: a1, a2
      real(4)   , intent(out) :: fice(720, 360)
 
      real(4) fice1(720,360),fice2(720,360)
      integer(4) ilon, ilat

!      write(*,*) seaice_file1
!      call mk_smooth_ice_map(seaice_file1, fice1)
      call read_seaice(seaice_file1, fice1)
!      write(*,*) seaice_file2
!      call mk_smooth_ice_map(seaice_file2, fice2)
      call read_seaice(seaice_file2, fice2)
 
      DO ILAT=1,180*2
         DO ILON=1,360*2
!     do not do xlat,xlon iterpolation for seaice.  Just do drop in the bucket
!            if(ilat.gt.-50 .and. ilat.lt.40) then 
!               fice(ilon,ilat)=0
!            else
               if(fice1(ilon,ilat).lt.0 .or. fice2(ilon,ilat).lt.0) then !this is an all-land cell, set fice to zero
                  fice(ilon,ilat)=0
               else
                  fice(ilon,ilat)=a1*fice1(ilon,ilat) + a2*fice2(ilon,ilat)
               endif
!            endif
         enddo
      enddo

      return
      end


      subroutine read_seaice(seaice_file, fice)
      implicit none

      integer, parameter :: DFACC_READ = 1

      character(100),  intent(in)  :: seaice_file
      real(4),     intent(out) :: fice(720,360)

      real(4) dummy_lat(360*2)
       integer*4 ilon, ilat

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 

      integer i

      sd_id_anc = sfstart( seaice_file, DFACC_READ)

      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = 360*2
      count(2) = 180*2

      sds_index = sfn2index(sd_id_anc, 'seaice_smoothed')
      sds_id = sfselect (sd_id_anc, sds_index)

      retn = sfrdata(sds_id, start, stride, count, fice)

!     Rotate by 180 degrees longitude
      do ilat=1,180*2
         dummy_lat=fice(:,ilat)
         do ilon=1,360*2
            i = mod((ilon-1+360),720)+1
            fice(ilon,ilat) = dummy_lat(i)
         enddo
      enddo

      retn = sfendacc(sds_id)
      retn = sfend(sd_id_anc) 

      return
      end

#if 0
!  this routine averages the ncep daily ice fields onto a 0.5 deg lat/lon grid
!  the averaging window is a circle with a radius of 50 km.

!  input:  filename1 is the file name for the ncep ice map.
!  the file contains a character (1-byte) array that is 4320 by 2160

      subroutine mk_smooth_ice_map(seaice_file, frac_ice_smoothed)
      implicit none

      real(8), parameter :: re=6378.137d3 !earth equatorial radius (meters)
      real(8), parameter :: rp=6356.752d3 !earth      polar radius (meters)
      real(8), parameter :: ffac=(rp/re)*(rp/re)

      integer(4), parameter :: nlat=2160, nlon=4320
      integer, parameter :: DFACC_READ = 1

      character(100),  intent(in)  :: seaice_file
      real(4),     intent(out) :: frac_ice_smoothed(720,360)
 
      character(1) char_ice(nlon,nlat)

      integer(4) istart
      integer(4) ilat,ilon,jlat,jlon,klat,klon,kklon

      real(4) xlat,xlatg,xlon
      real(4) sinlat0,coslat0
      real(4) cel0(3),dif(3),rsq,rsq_max,wt

      real(4) coslat(nlat),cel(3,nlon,nlat)
      integer(4) nstep(nlat)

      real(4) frac_ice(nlon,nlat), dummy_lat(nlon)
      real(8) tsum(0:1)

      real(4) cosd, sind, tand
      real(8) datand

      data istart/1/

      integer sfstart, sfrdata, sfn2index, sfselect
      integer sd_id_anc, sds_id, sds_index
      integer sfendacc, sfend
      integer start(2), count(2), stride(2)
      integer retn 

      if(istart.eq.1) then
         istart=0

         rsq_max=(2*sind(0.5*0.45d0))**2 !0.45 degree = 0.45 * 111.2 km/deg = 50 km radius

         do ilat=1,nlat
            xlat=(ilat-1)/12.-89.95833
            xlatg=datand(ffac*tand(xlat)) !convert from geodetic to geocentric
            sinlat0=sind(xlatg)
            coslat0=cosd(xlatg)

            coslat(ilat)=coslat0

            nstep(ilat)=nint(6./coslat0) + 1
            if(nstep(ilat).ge.nlon/2) nstep(ilat)=nlon/2 -1
            if(ilat.le.12 .or.ilat.ge.nlat-11)  nstep(ilat)=nlon/2 -1

            do ilon=1,nlon
               xlon=(ilon-1)/12.+ 0.04167
               cel(1,ilon,ilat)=cosd(xlon)*coslat0
               cel(2,ilon,ilat)=sind(xlon)*coslat0
               cel(3,ilon,ilat)=sinlat0
            enddo               !ilon
         enddo                  !ilat

      endif

      ! Read HDF4 file (byte array)
!      write(*,*) 'Reading ice file'
      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      count(1) = nlon
      count(2) = nlat

      sd_id_anc = sfstart( seaice_file, DFACC_READ)
      sds_index = sfn2index(sd_id_anc, 'seaice')
      sds_id = sfselect (sd_id_anc, sds_index)
      retn = sfrdata(sds_id, start, stride, count, char_ice)
      retn = sfendacc(sds_id)
      retn = sfend(sd_id_anc) 

      frac_ice=ichar(char_ice)*0.01


!     Reflect about the equator
!      do ilat=1,nlat/2
!         dummy_lat=frac_ice(:,ilat)
!         frac_ice(:,ilat) = frac_ice(:,nlat+1-ilat)
!         frac_ice(:,nlat+1-ilat) = dummy_lat
!      enddo


!     Rotate by 180 degrees longitude
      do ilat=1,nlat
         dummy_lat=frac_ice(:,ilat)
         do ilon=1,nlon
            frac_ice(ilon,ilat) = dummy_lat(mod((ilon+2160),4320))
         enddo
      enddo

!      write(*,*) 'Smoothing ice file'
      frac_ice_smoothed = 0
      do ilat=1,nlat

         xlat=(ilat-1)/12.-89.95833
         xlatg=datand(ffac*tand(xlat)) !convert from geodetic to geocentric
         sinlat0=sind(xlatg)
         coslat0=cosd(xlatg)

!         write(*,*) ilat

         do ilon=1,nlon

            xlon=(ilon-1)/12.+ 0.04167

            cel0(1)=cosd(xlon)*coslat0
            cel0(2)=sind(xlon)*coslat0
            cel0(3)=sinlat0

            tsum=0

            do klat=ilat-6,ilat+6
               if(klat.lt.1 .or. klat.gt.nlat) cycle
               wt=coslat(klat)

               do klon=ilon-nstep(klat),ilon+nstep(klat)+1
                  kklon=klon
                  if(kklon.lt.   1) kklon=kklon+nlon
                  if(kklon.gt.nlon) kklon=kklon-nlon

                  dif=cel(:,kklon,klat)-cel0

                  rsq=dif(1)*dif(1) + dif(2)*dif(2) + dif(3)*dif(3)

                  if(rsq.gt.rsq_max) cycle

                  tsum(0)=tsum(0) + wt
                  tsum(1)=tsum(1) + wt*frac_ice(kklon,klat)

               enddo            !klon
            enddo               !klat

            jlat=1 + int(2*(xlat+90))
            jlon=1 + int(2*xlon)
            if(jlat.lt.1 .or. jlat.gt.360) stop 'error1'
            if(jlon.lt.1 .or. jlon.gt.720) stop 'error2'

            frac_ice_smoothed(jlon,jlat)=tsum(1)/tsum(0)

         enddo                  !ilon
      enddo                     !ilat

      if(minval(frac_ice_smoothed).lt.0 .or. minval(frac_ice_smoothed).gt.1)  stop 'frac_ice_smooth oob, pgm stopped'
      
      return
      end
#endif

      subroutine read_solarflux(solarflux_file, solarflux_date, solarflux)
      implicit none

      integer, parameter :: DFACC_READ = 1

      character(128),  intent(in)  :: solarflux_file
      character(100),  intent(in)  :: solarflux_date
      real(4),     intent(out) :: solarflux
      character(9) date

!      integer sfstart, sfrdata, sfn2index, sfselect
!      integer sd_id_anc, sds_id, sds_index
!      integer sfendacc, sfend
!      integer start(1), count(1), stride(1)
!      integer retn 

      open(unit=30, file=solarflux_file, status='old', ERR=100)
      do
         read(30, 200) date, solarflux
         if ( date(1:8).eq.solarflux_date(1:8)) exit
      end do
 200  format(a9, f8.0)
      close(unit=30)
      return

 100  write (*,*) solarflux_file//' not found.'
      call EXIT(1)

      stop

!      sd_id_anc = sfstart( solarflux_file, DFACC_READ)
!      if (sd_id_anc.eq.-1) then
!         write (*,*) solarflux_file//' not found.'
!         stop
!      endif

!      start(1) = 0
!      stride(1) = 1
!      count(1) = 1

!      sds_index = sfn2index(sd_id_anc, 'solar_flux')
!      sds_id = sfselect (sd_id_anc, sds_index)
!      retn = sfrdata(sds_id, start, stride, count, solarflux)
!      retn = sfendacc(sds_id)
!      retn = sfend(sd_id_anc) 

      return
      end




C***************************  SUBROUTINE SWAP_BYTES  ***************************
C
C       this subroutine swaps the order of bytes in an input array.
C       (same as BSWAP in IIS intrinsics)
C
C       Parameters:
C         BUF   >< i1a  Array to have its bytes swapped.
C         NBYTE >  i4   number of bytes to be swapped
C         COUNT >  i4   Number of iteration for the swapping
C
C       Created by Gary Fu, GSC, 12/98
C*************************************************************************
C
        subroutine SWAP_BYTES(BUF)
        byte BUF(4),tmp

        tmp = buf(4)
        buf(4) = buf(1)
        buf(1) = tmp

        tmp = buf(3)
        buf(3) = buf(2)
        buf(2) = tmp

        return
        end

