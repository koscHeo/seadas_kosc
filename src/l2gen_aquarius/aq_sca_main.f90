Integer Function jday (Year,Month,Day)

        ! computes the Julian date given a gregorian calendar date (year,month,day)
        Integer Year,Month,Day,Yi,Mj,Dk
        Yi= Year
        Mj= Month
        Dk= Day
        jday= Dk-32075+1461*(Yi+4800+(Mj-14)/12)/4+367*(Mj-2-(Mj-14)/12*12)/12-3*((Yi+4900+(Mj-14)/12)/100)/4
        Return

End Function jday

Subroutine gdate (Julianday, Year,Month,Day)

        ! computes the gregorian calendar date (year,month,day) given the julian date
        Integer Julianday,Year,Month,Day,Yi,Mj,Dk
        L= Julianday+68569
        N= 4*L/146097
        L= L-(146097*N+3)/4
        Yi= 4000*(L+1)/1461001
        L= L-1461*Yi/4+31
        Mj= 80*L/2447
        Dk= L-2447*Mj/80
        L= Mj/11
        Mj= Mj+2-12*L
        Yi= 100*(N-49)+Yi+L
        Year= Yi
        Month= Mj
        Day= Dk
        Return

End Subroutine gdate

program sca_main

        ! Implements single-channel algorithm (SCA) for soil moisture retrieval
        ! using Aquarius observations (HDF5 from PO.DAAC)

        ! use hdf5 in fortran program, to add the hdf5.mod module file to
        ! compiler includes.
        
        ! Developed by Rajat Bindlish, Tom Jackson and Tianjie Zhao at USDA-ARS Hydrology and Remote Sensing Lab
        ! sm_flag order No SM retrieval, TB, Orbit Maneuver, RFI, Tsurf, Frozen Ground, Snow, Ice, NDVI quality, Dense Veg, Urban, Soil Texture, Water

        ! Change sm_flag+bflags(n)*(2**n) to sm_flag+bflags(n)*(2**(n-1)) (RB 08/26/13)
        ! Note: must adjust Dense Vegetation masking in l2gen_aquarius from 1024 to 512
        ! Use L2 NAV flag to set SM NAV flag

      Use HDF5

      Implicit none
      Integer                                   :: jday

      ! Retrieving
      Byte                      :: flg, lct
      Byte, Dimension(16)       :: bflags
      Integer*2                 :: sm_flag
      Integer*4                 :: i, j, k, n
      Integer*4                                 :: ii, jj, kk
      Integer*4                                 :: yy, mm, dd, jd, nyr, doy, ddy
      Real                                      :: sandp, clayp, poro, xporo, bukd
      Real                                      :: ang, te, tbh, tbv, hh, hv, vh, vv
      Real                                      :: ssa, mdv, ndv, nndv, vwc, bfac, ft, ntau, minndvi
      Real                                      :: mv, hr
      Real, Dimension(3)        :: angs
      data angs /29.358, 38.487, 46.286/

      ! Reading
      Real, Dimension(7200)                     :: slt1d
      Real, Dimension(7200,3600)                :: snd2d, cly2d
      Real, Dimension(7200)                     :: bkd1d
      Real, Dimension(7200,3600)                :: bkd2d
      Integer*2, Dimension(7200)                :: ndv1d
      Integer*2, Dimension(7200,3600)           :: ndv2d, nndv2d, mdv2d, minndvi2d
      Byte, Dimension(7200)                     :: lct1d, flg1d
      Byte, Dimension(7200,3600)                :: lct2d, flg2d

      ! Aquarius
      Integer(HID_T)                            :: file_id, nfile_id ! File identifier.
      Integer(HID_T)                            :: attr_id, atype_id ! Attribute identifier.
      INTEGER(HSIZE_T)                          :: attrlen    ! Length of the attribute string
      INTEGER(SIZE_T)                          :: dattrlen    ! Length of the attribute string
      INTEGER(HSIZE_T), DIMENSION(1)            :: attrdims
      Integer(HID_T)                            :: sec_id
      Integer(HID_T)                            :: lat_id, lon_id, zng_id
      Integer(HID_T)                            :: tbh_id, tbv_id ! Dataset identifier.
      Integer(HID_T)                            :: hh_id, hv_id, vh_id, vv_id ! Dataset identifier.
      Integer(HID_T)                            :: vsm_id, lst_id, lkt_id
      Integer(HID_T)                            :: lfr_id, ifr_id
      Integer(HID_T)                            :: phi_id, rpy_id, acs_id, rflag_id
      Integer                                   :: hdferr ! Error flag.

      Integer(HID_T)                            :: space, secspace ! Space handle.
      Integer(HSIZE_T), Dimension(1:2)          :: dims, secdims
      Integer(HSIZE_T), Dimension(1:2)          :: ndims, secndims
      Integer(HSIZE_T), Dimension(3)            :: rflagdims
      Real*8, Dimension(5000)                   :: sec1d
      Integer*1                                 :: acs1d

      Real*8                                    :: secGPS
      Real, Dimension(3,5000)                   :: lat2d, lon2d, tbh2d, tbv2d, hh2d, hv2d, vh2d, vv2d
      Real, Dimension(3,5000)                   :: vsm2d, lst2d, lfr2d, ifr2d, lkt2d  ! Read buffer.
      Real, Dimension(3,5000)                   :: phi2d, rpy2d, swe2d
      integer*4, Dimension(4,3,5000)            :: rflag3d
      Real*8, Dimension(5000)                   :: zang1d
      Real                            :: lat, lon, lfr, ifr, phi, roll, pitch, yaw, swe, sec, lkt, lst, vsm
      Integer*1                       :: acs
      Real*8                          :: zang
      Integer(HID_T)                  :: swe_id
      integer*4                       :: nav

      character ver_n*5

      character outfile*320
      character inyr*4,indt*3,inff*26,ser0*4,ser*4,mmstr*2,tend*1,nser0*4,nser*4
      character anomaly_status*128, nominalNav*8
      LOGICAL attr_exists

      character*320 read_ancillary_dir
      character*320 input_L2_file
      character*15 starttime
      character*15 endtime
      character*17 nodeCrossingTime
      real*4       nodeLongitude
      integer*4    orbitNumber, cycleNumber, passNumber
      character*512 ancillary_files

      integer i_len             ! typically used to hold the length of a string
      integer lnblnk
      logical*1 i_exist

      ! Variables used to parse command line arguments
      integer optind
      character*320 str_argv
      character*320 string

      ver_n    ='3.0.0'  ! version number

      ! Parse command-line arguments      
      optind = 1
      call getarg(optind,str_argv)
      do while( optind .le. iargc() .and. str_argv(1:1) .eq. '-')
         i_len = lnblnk(str_argv)

         if( str_argv(1:i_len) .eq. '-input_L2_file' )then
            optind = optind + 1
            call getarg(optind, input_L2_file )
        
         elseif( str_argv(1:i_len) .eq. '-read_ancillary_dir' )then
            optind = optind + 1
            call getarg(optind, read_ancillary_dir )
        
         elseif( str_argv(1:i_len) .eq. '-outfile' )then
            optind = optind + 1
            call getarg(optind, outfile )          
   
         endif
         optind = optind + 1
         call getarg(optind,str_argv)
      enddo

      ! Check if input L2 file exists
      INQUIRE(FILE=input_L2_file, EXIST=I_EXIST)
      if (i_exist.eqv..false.) then
         i_len = lnblnk(input_L2_file)
         write(*,*) 'Problem opening ' // input_L2_file(1:i_len)
         call exit(1)
      endif

      write(*,*) 'Reading soil data ...'
      open(1,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/Aquarius_sand_p05d.bin',status='old', &
           convert='little_endian',form='unformatted',access='direct',recl=28800)
      do i=1,3600
         read(1,rec=i) slt1d
         do j=1,7200
            snd2d(j,i) = slt1d(j)
         end do
      end do
      close(1)
      open(1,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/Aquarius_clay_p05d.bin',status='old', &
           convert='little_endian',form='unformatted',access='direct',recl=28800)
      do i=1,3600
         read(1,rec=i) slt1d
         do j=1,7200
            cly2d(j,i) = slt1d(j)
         end do
      end do
      close(1)
      open(1,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/Aquarius_bd_p05d.bin',status='old',  &
           convert='little_endian',form='unformatted',access='direct',recl=28800)
      do i=1,3600
         read(1,rec=i) bkd1d
         do j=1,7200
            bkd2d(j,i) = bkd1d(j)
         end do
      end do
      close(1)
      write(*,*) 'Reading soil data is done.'

      write(*,*) 'Reading land cover ...'
      open(1,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/Aquarius_land_cover_5km.dat', &
           convert='little_endian',status='old',form='unformatted',access='direct',recl=7200)
      do i=1,3600
         read(1,rec=i) lct1d
         do j=1,7200
            lct2d(j,i) = lct1d(j)
         end do
      end do
      close(1)
      write(*,*) 'Reading land cover is done.'

      write(*,*) 'Reading MODIS/NDVI max ...'
      open(4,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/ndvi_max.dat', &
           convert='little_endian',status='old',form='unformatted',access='direct',recl=14400)
      do ii=1,3600
         read(4,rec=ii) ndv1d
         do jj=1,7200
            mdv2d(jj,ii)=ndv1d(jj)
         enddo
      enddo
      close(4)
      write(*,*) 'Reading MODIS/NDVI max is done.'

      write(*,*) 'Reading MODIS/NDVI min ...'
      open(4,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/ndvi_min.dat', &
           convert='little_endian',status='old',form='unformatted',access='direct',recl=14400)
      do ii=1,3600
         read(4,rec=ii) ndv1d
         do jj=1,7200
            minndvi2d(jj,ii)=ndv1d(jj)
         enddo
      enddo
      close(4)
      write(*,*) 'Reading MODIS/NDVI min is done.'

      ! Read Aquarius data from unpacked HDF5 data files listed:
      call H5open_f(hdferr)
      i_len = INDEX(input_L2_file, '/', .true.) 

      ser0='00_0'
      nser0='00_0'

      inff=input_L2_file(i_len+1:i_len+26)
      inyr=inff(2:5)
      indt=inff(6:8)
      write(*,*) inff

      read(inyr,'(i4)') nyr
      read(indt,'(i3)') doy
      jd=jday(nyr,1,doy)
      call gdate(jd-4,yy,mm,dd)
      if (mm.eq.1) mmstr='01'
      if (mm.eq.2) mmstr='02'
      if (mm.eq.3) mmstr='03'
      if (mm.eq.4) mmstr='04'
      if (mm.eq.5) mmstr='05'
      if (mm.eq.6) mmstr='06'
      if (mm.eq.7) mmstr='07'
      if (mm.eq.8) mmstr='08'
      if (mm.eq.9) mmstr='09'
      if (mm.eq.10) mmstr='10'
      if (mm.eq.11) mmstr='11'
      if (mm.eq.12) mmstr='12'

      if (dd.le.10) then
         tend='1'
         ddy=dd-1
      endif
      if (dd.gt.10.and.dd.le.20) then
         tend='2'
         ddy=dd-11
      endif
      if (dd.gt.20) then
         tend='3'
         ddy=dd-21
      endif

      ser=mmstr//'_'//tend
      if (ser0.ne.ser) then
         write(*,*) 'Reading MODIS/NDVI '//ser//' ...'
         open(4,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/ndvi_cmg_'//ser//'.dat', &
              convert='little_endian',status='old',form='unformatted',access='direct',recl=14400)
         do ii=1,3600
            read(4,rec=ii) ndv1d
            do jj=1,7200
               ndv2d(jj,ii)=ndv1d(jj)
            enddo
         enddo
         close(4)
         open(5,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/flag_cmg_'//ser//'.dat', &
              convert='little_endian',status='old',form='unformatted',access='direct',recl=7200)
         do ii=1,3600
            read(5,rec=ii) flg1d
            do jj=1,7200
               flg2d(jj,ii)=flg1d(jj)
            enddo
         enddo
         close(5)
         write(*,*) 'Reading MODIS/NDVI is done.'
      endif
      ser0=ser

      jd=jday(nyr,1,doy)
      call gdate(jd+6,yy,mm,dd)
      if (mm.eq.1) mmstr='01'
      if (mm.eq.2) mmstr='02'
      if (mm.eq.3) mmstr='03'
      if (mm.eq.4) mmstr='04'
      if (mm.eq.5) mmstr='05'
      if (mm.eq.6) mmstr='06'
      if (mm.eq.7) mmstr='07'
      if (mm.eq.8) mmstr='08'
      if (mm.eq.9) mmstr='09'
      if (mm.eq.10) mmstr='10'
      if (mm.eq.11) mmstr='11'
      if (mm.eq.12) mmstr='12'
      if (dd.le.10) tend='1'
      if (dd.gt.10.and.dd.le.20) tend='2'
      if (dd.gt.20) tend='3'
      
      nser=mmstr//'_'//tend
      if (nser0.ne.nser) then
         write(*,*) 'Reading MODIS/NDVI '//nser//' ...'

         open(4,file=read_ancillary_dir(1:lnblnk(read_ancillary_dir))//'/ndvi_cmg_'//nser//'.dat', &
              convert='little_endian',status='old',form='unformatted',access='direct',recl=14400)
         do ii=1,3600
            read(4,rec=ii) ndv1d
            do jj=1,7200
               nndv2d(jj,ii)=ndv1d(jj)
            enddo
         enddo
         close(4)
         write(*,*) 'Reading MODIS/NDVI is done.'
      endif
      nser0=nser

      write(*,*) outfile(1:lnblnk(outfile))
      open(3,file=outfile(1:lnblnk(outfile)),status='unknown', &
           form='formatted',action='write')

      ! Open an existing file.
     call H5Fopen_f(input_L2_file, H5F_ACC_RDWR_F, file_id, hdferr)

     ! Read Start Time, Stop Time attributes
      dattrlen = 13
      attrlen = 13
      attrdims(1) = 1
      call H5Tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
      call H5Tset_size_f(atype_id, dattrlen, hdferr)

      call H5Aopen_name_f(file_id, 'Start Time', attr_id, hdferr)
      call H5Aread_f(attr_id, atype_id, starttime, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 
      call H5Aopen_name_f(file_id, 'End Time', attr_id, hdferr)
      call H5Aread_f(attr_id, atype_id, endtime, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 

      dattrlen = 16
      attrlen = 16
      call H5Tset_size_f(atype_id, dattrlen, hdferr)
      call H5Aopen_name_f(file_id, 'Node Crossing Time', attr_id, hdferr)
      call H5Aread_f(attr_id, atype_id, nodeCrossingTime, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 

      call H5Aopen_name_f(file_id, 'Nominal Navigation', attr_id, hdferr)
      call h5aget_storage_size_f(attr_id, attrlen, hdferr)
      call H5Tset_size_f(atype_id, dattrlen, hdferr)
      call H5Aread_f(attr_id, atype_id, nominalNav, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 

      anomaly_status = 'N/A'
      call h5aexists_by_name_f(file_id, '.', 'Anomaly Status String', attr_exists, hdferr)
      if (attr_exists.eqv..TRUE.) then
         call H5Aopen_name_f(file_id, 'Anomaly Status String', attr_id, hdferr)
         call h5aget_storage_size_f(attr_id, attrlen, hdferr)
         call H5Tset_size_f(atype_id, dattrlen, hdferr)
         call H5Aread_f(attr_id, atype_id, anomaly_status, attrdims, hdferr)
         call H5Aclose_f(attr_id, hdferr) 
      endif

      call H5Aopen_name_f(file_id, 'Orbit Node Longitude', attr_id, hdferr)
      CALL h5aread_f(attr_id, H5T_NATIVE_REAL, nodeLongitude, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 

      call H5Aopen_name_f(file_id, 'Orbit Number', attr_id, hdferr)
      CALL h5aread_f(attr_id, H5T_STD_I32LE, orbitNumber, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 
      call H5Aopen_name_f(file_id, 'Cycle Number', attr_id, hdferr)
      CALL h5aread_f(attr_id, H5T_STD_I32LE, cycleNumber, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 
      call H5Aopen_name_f(file_id, 'Pass Number', attr_id, hdferr)
      CALL h5aread_f(attr_id, H5T_STD_I32LE, passNumber, attrdims, hdferr)
      call H5Aclose_f(attr_id, hdferr) 

      ancillary_files = 'Aquarius_sand_p05d.bin,Aquarius_clay_p05d.bin,Aquarius_bd_p05d.bin,'
      ancillary_files = trim(ancillary_files)//'Aquarius_land_cover_5km.dat,ndvi_max.dat,'
      ancillary_files = trim(ancillary_files)//'ndvi_min.dat,ndvi_cmg_'//ser//'.dat,'
      ancillary_files = trim(ancillary_files)//'flag_cmg_'//ser//'.dat,ndvi_cmg_'//nser//'.dat'

      ! Open an existing dataset.
      call H5Dopen_f(file_id, "/Block Attributes/secGPS", sec_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Navigation/beam_clat", lat_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Navigation/beam_clon", lon_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Navigation/zang", zng_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Navigation/celphi", phi_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Navigation/att_ang", rpy_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
!      call H5Dopen_f(file_id, "/Navigation/acs_mode", acs_id, hdferr)
!      if (hdferr.ne.0) then
!         write(*,*) inff//' is a bad file!'
!         call exit(1)
!      end if
      call H5Dopen_f(file_id, "/Aquarius Data/rad_TbH", tbh_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/rad_TbV", tbv_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/scat_HH_toa", hh_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/scat_HV_toa", hv_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/scat_VH_toa", vh_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/scat_VV_toa", vv_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/anc_sm", vsm_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/anc_surface_temp", lkt_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/anc_subsurf_temp", lst_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/anc_swe", swe_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/rad_land_frac", lfr_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Data/rad_ice_frac", ifr_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dopen_f(file_id, "/Aquarius Flags/radiometer_flags", rflag_id, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      ! The same dimentions of Tb and LatLon data.
      call H5Dget_space_f(lat_id, space, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Sget_simple_extent_dims_f(space, dims, ndims, hdferr)

      rflagdims(1) = 4
      rflagdims(2) = 3
      rflagdims(3) = dims(2)

      secdims(1) = dims(2)
      !       if (hdferr.ne.0) then
      !               write(*,*) inff//' is a bad file!','12'
      !                     !       end if
      ! Read the dataset.
      call H5Dread_f(sec_id, H5T_NATIVE_DOUBLE, sec1d, secdims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(lat_id, H5T_NATIVE_REAL, lat2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(lon_id, H5T_NATIVE_REAL, lon2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(zng_id, H5T_NATIVE_DOUBLE, zang1d, secdims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(phi_id, H5T_NATIVE_REAL, phi2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(rpy_id, H5T_NATIVE_REAL, rpy2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
!      call H5Dread_f(acs_id, H5T_STD_I8LE, acs1d, secdims, hdferr)
!      if (hdferr.ne.0) then
!         write(*,*) inff//' is a bad file!'
!         call exit(1)
!      end if
      call H5Dread_f(tbh_id, H5T_NATIVE_REAL, tbh2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) 'Could not read rad_TbH data in '//inff
      end if
      call H5Dread_f(tbv_id, H5T_NATIVE_REAL, tbv2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(hh_id, H5T_NATIVE_REAL, hh2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(hv_id, H5T_NATIVE_REAL, hv2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(vh_id, H5T_NATIVE_REAL, vh2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
      end if
      call H5Dread_f(vv_id, H5T_NATIVE_REAL, vv2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(vsm_id, H5T_NATIVE_REAL, vsm2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(lkt_id, H5T_NATIVE_REAL, lkt2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(lst_id, H5T_NATIVE_REAL, lst2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(swe_id, H5T_NATIVE_REAL, swe2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(lfr_id, H5T_NATIVE_REAL, lfr2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(ifr_id, H5T_NATIVE_REAL, ifr2d, dims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      call H5Dread_f(rflag_id, H5T_STD_I32LE, rflag3d, rflagdims, hdferr)
      if (hdferr.ne.0) then
         write(*,*) inff//' is a bad file!'
         call exit(1)
      end if
      ! Close the dataset.
      call H5Dclose_f(sec_id, hdferr)
      call H5Dclose_f(lat_id, hdferr)
      call H5Dclose_f(lon_id, hdferr)
      call H5Dclose_f(zng_id, hdferr)
      call H5Dclose_f(phi_id, hdferr)
      call H5Dclose_f(rpy_id, hdferr)
!      call H5Dclose_f(acs_id, hdferr)
      call H5Dclose_f(tbh_id, hdferr)
      call H5Dclose_f(tbv_id, hdferr)
      call H5Dclose_f(hh_id, hdferr)
      call H5Dclose_f(hv_id, hdferr)
      call H5Dclose_f(vh_id, hdferr)
      call H5Dclose_f(vv_id, hdferr)
      call H5Dclose_f(vsm_id, hdferr)
      call H5Dclose_f(lkt_id, hdferr)
      call H5Dclose_f(lst_id, hdferr)
      call H5Dclose_f(swe_id, hdferr)
      call H5Dclose_f(lfr_id, hdferr)
      call H5Dclose_f(ifr_id, hdferr)
      call H5Dclose_f(rflag_id, hdferr)
      ! Close the file.
      CALL H5Fclose_f(file_id, hdferr)

      ! default parameters:
      ssa   = 0.05
      bfac  = 0.12
      hr    = 0.1 
      ntau  = 0.0

      ! Note: The number following 'a' is one less than the character
      ! buffer size.
      write(3,'(a319)') input_L2_file
      write(3,'(a13)') starttime
      write(3,'(a13)') endtime
      write(3,'(a16)') nodeCrossingTime
      write(3,'(f10.5)') nodeLongitude
      write(3,'(i6)') orbitNumber
      write(3,'(i6)') cycleNumber
      write(3,'(i4)') passNumber
      write(3,'(a)') nominalNav
      write(3,'(a)') anomaly_status

      i_len = lnblnk(ancillary_files)
      write(3,'(a)') ancillary_files(1:i_len)

      do k=1, dims(2)
         do kk=1, dims(1)

            do n=1,16
               bflags(n)=0                ! 0 for retrieval and 1 for flagging
            end do
            
            secGPS= sec1d(k)            ! Sec
            lat   = lat2d(kk,k)         ! Latitude
            lon   = lon2d(kk,k)         ! Longitude
            zang  = zang1d(k)           ! Intra-Orbit Angle
            phi   = phi2d(kk,k)         ! Azimuth
            roll  = rpy2d(1,k)          ! spacecraft roll
            pitch  = rpy2d(2,k)         ! spacecraft pitch
            yaw  = rpy2d(3,k)           ! spacecraft yaw
!            acs  = acs1d(k)             ! ACS Mode
            tbh    = tbh2d(kk,k)        ! Aquarius observation (3 beams)
            tbv    = tbv2d(kk,k)        ! Aquarius observation (3 beams)
            hh    = hh2d(kk,k)        ! Aquarius observation (3 beams)
            hv    = hv2d(kk,k)        ! Aquarius observation (3 beams)
            vh    = vh2d(kk,k)        ! Aquarius observation (3 beams)
            vv    = vv2d(kk,k)        ! Aquarius observation (3 beams)
            vsm   = vsm2d(kk,k)        ! NCEP Soil Moisture
            if (vsm.gt.1.0) vsm=1.0
            if (lst2d(kk,k).lt.500) then
               te    = (lkt2d(kk,k)+lst2d(kk,k))/2.0         ! NCEP soil temperature
            else
               te=lkt2d(kk,k)
            endif

            lfr   = lfr2d(kk,k)         ! land fraction
            ifr   = ifr2d(kk,k)         ! ice fraction
            ang   = angs(kk)
            
            swe   = swe2d(kk,k)         ! SWE

            nav   = or(rflag3d(1,kk,k), rflag3d(2,kk,k))
            nav   = or(rflag3d(3,kk,k), nav)
            nav   = or(rflag3d(4,kk,k), nav)

!           Removed in V3.0
!            if (ver_n.eq.'2.0.0') then
!               if (kk.eq.1) tbh=(tbh+2.0714)/1.0416
!               if (kk.eq.2) tbh=(tbh+1.0819)/1.0390
!               if (kk.eq.3) tbh=(tbh+0.8105)/1.0396
!            endif

            if ((lat.le.90).and.(lat.ge.-90).and.(lon.le.180).and.(lon.ge.-180)) then
               sandp = snd2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)*100.0
               clayp = cly2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)*100.0
               bukd  = bkd2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)
               poro  = (1-bukd/2.65)*100.0
               lct   = lct2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)

               ft=0.0
               if (lct.eq.1) ft=15.96
               if (lct.eq.2) ft=19.15
               if (lct.eq.3) ft=7.98
               if (lct.eq.4) ft=12.77
               if (lct.eq.5) ft=12.77
               if (lct.eq.6) ft=3.0
               if (lct.eq.7) ft=1.5
               if (lct.eq.8) ft=4.0
               if (lct.eq.9) ft=3.0
               if (lct.eq.10) ft=1.5
               if (lct.eq.11) ft=4.0
               if (lct.eq.12) ft=3.5
               if (lct.eq.13) ft=6.49
               if (lct.eq.14) ft=3.25
               if (lct.eq.15) ft=0.0
               if (lct.eq.16) ft=0.0
               if ((lct.eq.7).or.(lct.eq.10).or.(lct.eq.16)) hr=0.3
               if ((lct.eq.7).or.(lct.eq.10).or.(lct.eq.16)) bfac=0.13

!               if (ft.lt.2.0) bfac=0.1

               ndv   = ndv2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)*0.0001
               nndv  = nndv2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)*0.0001
               mdv   = mdv2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)*0.0001
               minndvi = minndvi2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)*0.0001
               flg   = flg2d(int((lon+180)/0.05)+1,int((90-lat)/0.05)+1)
               
               if (ndv.lt.minndvi) ndv=minndvi
               if (nndv.lt.minndvi) nndv=minndvi
               
               ndv=ndv+(nndv-ndv)*ddy/10.0

               if (lct.le.5) ft=ft/2.0
               vwc   = 1.9134*ndv*ndv-0.3215*ndv+ft*(mdv-0.1)/0.9
!               if (vwc.gt.5.0) vwc=5.0+(vwc-5.0)/2.0
               if ((lct.eq.10).or.(lct.eq.12).or.(lct.eq.14)) vwc=1.9134*ndv*ndv-0.3215*ndv+ft*(ndv-0.1)/0.9
               if (vwc.lt.0) vwc=0.0
               ntau=vwc*bfac
                        
            else
               bflags(1)=1
               sandp=-999.0
               clayp=-999.0
               bukd=-999.0
               poro=-999.0
               lct=255
               ndv=-999.0
               mdv=-999.0
               flg=255
               vwc=-999.0
            endif

            if (ndv.lt.0) then
               vwc=-999.0
               ntau=-999.0
            endif
                        
            if (roll.gt.1 .or. pitch.gt.1 .or. yaw.gt.5) bflags(3)=1  ! spacecraft roll pitch yaw
!            if (acs.ne.5) bflags(3)=1  ! spacecraft ACS Mode
            if (btest(nav,12) .eqv. .true.) bflags(3)=1

            if (ndv.lt.0) bflags(9)=1                                 ! invalid ndvi
!            if (flg.ne.0) bflags(9)=1
            if (sandp.gt.100 .or. sandp.lt.0) bflags(12)=1            ! invalid soil data
            if (clayp.gt.100 .or. clayp.lt.0) bflags(12)=1            ! invalid soil data  
            if (poro.gt.100.0 .or. poro.lt.0) bflags(12)=1            ! invalid soil data
            if (lct.eq.0) bflags(13)=1                   ! water
            if (lct.eq.13) bflags(11)=1                  ! urban
            if (lct.eq.15) bflags(7)=1                                ! snow and ice from landcover
            if (lfr.lt.0.99) bflags(13)=1                             ! land fraction from Aq
            if (ifr.gt.0.1) bflags(8)=1                               ! ice
            if ((tbh.le.100).or.(tbh.gt.320)) bflags(4)=1                   ! invalid data
            if ((tbh.le.100).or.(tbh.gt.320)) bflags(2)=1                   ! invalid data
            if (te.le.273.15) bflags(6)=1                             ! invalid data

            if (lkt2d(kk,k).le.273.15) bflags(6)=1                    ! invalid data
            if (lst2d(kk,k).le.273.15) bflags(6)=1   
            if (lst2d(kk,k).ge.1000.00) bflags(14)=1                             ! invalid data

!            if (vwc.gt.5.0) bflags(10)=1                             ! Dense Vegetation
            if (te.le.tbh) bflags(5)=1                             ! invalid Tsurf data
            if (swe.gt.10.0) bflags(7)=1                             ! Snow

            sm_flag=0
            do n=1,16
!               sm_flag=sm_flag+bflags(n)*(2**n) ! convert binary bites to decimal
               sm_flag=sm_flag+bflags(n)*(2**(n-1)) ! convert binary bites to decimal (RB 08/26/13)
            enddo

            if (sm_flag.ne.0) then
               mv=-9999.00
            else
               call tb2vsm(ang,tbh,te,sandp,clayp,poro,vwc,bfac,ntau,ssa,hr,mv)
               xporo = poro/100.
               if (mv.gt.xporo) mv = xporo ! validness check
            endif
            
            if (vwc.gt.5.0) then
               bflags(10)=1                       ! Dense Vegetation
!               sm_flag=sm_flag+bflags(10)*(2**10) (RB 08/26/13)
               sm_flag=sm_flag+bflags(10)*(2**(10-1))
            endif

            if (mv.lt.0) bflags(1)=1
            if (mv.lt.0) sm_flag=sm_flag+1

            if (lst2d(kk,k).ge.1000.00) lst2d(kk,k)=-999.00                             ! invalid data
            if (vsm.ge.1000.00) vsm=-999.00                             ! invalid data
            
            !                       if (mv.gt.0.0) then
            write(3,'(i5,f15.3,4f10.3,5f12.4,i8)') kk,secGPS,lat,lon,zang,phi, &
                 tbh,tbv,lkt2d(kk,k),lst2d(kk,k),mv,sm_flag
            write(3,'(5f10.5,6f10.3)') roll,pitch,yaw,lfr,ifr,swe,vsm,hv,vv,vh,hh

            ! kk - Beam Number
            ! lat - Latitude
            ! lon - Longitude
            ! zang - intra-orbit angle
            ! phi - Azimuth angle (used to determine the direction of the orbit (Asc/Dsc) phi > 270 degrees is Asc
            ! tbh - Aquarius h pol radiometer observations (K)
            ! tbv - Aquarius v pol radiometer observations (K)
            ! lkt2d - NCEP surface temperature (K)
            ! lst2d - NCEP subsurface (0-10 cm) (temperature (K)
            ! mv - Volumetric Soil Moisture (m3/m3)
            ! sm_flag - Bit flag for soil moisture retrievals
            ! vsm - NCEP soil moisture (m3/m3)
            ! roll,pitch,yaw - Spacecraft Roll,Pitch,Yaw
            ! lfr - gland fraction
            ! ifr - gice fraction
            ! swe - snow water equ
            ! vsm - ancillary soil moisture
            ! hh,hv,vh,vv - TOA Scatterometer NRCS

         enddo
      enddo

      close(3)

      ! moves on to next swath data
      99    close(1)
      
      write(*,*) '----------'

end program sca_main

