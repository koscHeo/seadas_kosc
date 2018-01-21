module test_nan
  implicit none
  private
  public::is_nan
  interface is_nan
     module procedure is_nan_r4,is_nan_r8
  end interface is_nan
contains
  pure logical function is_nan_r4(a)
    implicit none
    real, intent(in)::a
    if (a.ne.a) then
       is_nan_r4=.true.
    else
       is_nan_r4=.false.
    endif
    return
  end function is_nan_r4
  
  pure logical function is_nan_r8(a)
    implicit none
    real (kind=8), intent(in)::a
    if (a.ne.a) then
       is_nan_r8=.true.
    else
       is_nan_r8=.false.
    endif
    return
  end function is_nan_r8
end module test_nan

module pos_vel_quat
  use params, only: clen
  implicit none
  private


  real(kind=8),allocatable,dimension(:,:),public,protected::pvq_data  
  type header
     character (len=clen)::name, source, subset, outname
     character (len=24)::epoch
     integer::req_pix
     character (len=8)::units
     character (len=5)::central_body,coordinate_system
     real (kind=8)::theta_from_stowed
     character (len=72)::code_name, code_version, code_date, code_author, &
          code_runon, code_idl_family, code_idl_os, code_idl_osname,  &
          code_idl_version, code_start_time, code_end_time,code_runby
     real (kind=8):: exposure_interval
     character (len=4)::iss_orientation
     real (kind=8)::trigger_pulse_width, odrc_gps_seconds
  end type header
  public::header
  public::read_pos_vel_quat


contains
  
  subroutine read_pos_vel_quat(filename,info,lines)
    use params, only:clen
    implicit none
    ! arguments
    character (len=clen),intent(in)::filename
    type (header),intent(out)::info
    integer,intent(out)::lines
!other globals (i.e., shared via modules)

! locals
    integer::ios,icom,ll,i,lun
    character (len=clen)::junk,key,val
    
    ! get lun from my lun finder (tafkaa)
    open(newunit=lun,file=trim(filename),access='sequential',action='read',&
         status='old',iostat=ios)
    
    do i=1,35
       junk=''
       key=''
       val=''
       read(lun,fmt='(a)')junk
       junk=adjustl(junk)
       ll=len_trim(junk)
       if (ll.eq.0) cycle   ! empty
       icom=index(junk,',')
       if (icom.eq.0) cycle ! non valid header
       key(:icom-1)=junk(:icom-1)
       val=trim(adjustl(junk(icom+1:ll)))
       select case(trim(key))
          case ('This file name')
             info%name=val
          case ('Source CSV file')
             info%source=val
          case('Subset CSV file')
             info%subset=val
          case ('Expected Lon/Lat/View angle filename')
             info%outname=val
          case ('Epoch')
             info%epoch=val
          case ('Requested number of pixels')
             read(val,fmt='(i4)')info%req_pix
          case ('Distance Unit')
             info%units=val
          case('Central Body')
             info%central_body=val
          case('CoordinateSystem')
             info%coordinate_system=val
          case('Theta (degrees from stowed position)')
             read(val,*)info%theta_from_stowed
          case('Code name')
             info%code_name=val
          case('Code version')
             info%code_version=val
          case('Code date')
             info%code_date=val
          case('Code author')
             info%code_author=val
          case('Code executed on computer')
             info%code_runon=val
          case('Code executed by username')
             info%code_runby=val
          case ('Code run under IDL osfamily')
             info%code_idl_family=val
          case('Code run under IDL os')
             info%code_idl_os=val
          case('Code run under IDL osname')
             info%code_idl_osname=val
          case('Code run under IDL version')
             info%code_idl_version=val
          case('Code start time')
             info%code_start_time=val
          case('Code end time')
             info%code_end_time=val
          case('Exposure interval (frame time)')
             read(val,*)info%exposure_interval
          case('ISS orientation')
             info%iss_orientation=val
          case('Trigger pulse width (s)')
             read(val,*)info%trigger_pulse_width
          case('ODRC broadcast time - gps time (s)')
             read(val,*)info%odrc_gps_seconds
             if (abs(info%odrc_gps_seconds).ge.2.0) then
                print*,'ODRC Broadcast time - gps time too large?'
                print*,' = ',info%odrc_gps_seconds
                print*,'Exiting...'
                call exit(110)
             endif
       end select
    end do
    read(lun,fmt='(a)')junk ! header
    
    ios=0
    i=0
    junk=''
    do 
       read(lun,fmt='(a)',iostat=ios)junk
       if (ios.ne.0) exit
       i=i+1
    end do
    lines=i
    rewind(lun)
    do i=1,36
       read(lun,fmt='(a)')junk
    end do
    allocate(pvq_data(15,lines))
    do i=1,lines
       read(lun,*)pvq_data(:,i)
    end do
    close(lun)
   

  end subroutine read_pos_vel_quat
end module pos_vel_quat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program proc_hico
  use test_nan
  use params,only:clen
  use pos_vel_quat,only: read_pos_vel_quat,header,pvq_data
  use quat_to_rot,only: q_to_r
  use sol_geo,only:solar_geometry
  use ecef2llh,only:ecef2latlon
  use envi_header_mod,only:write_envi_header
  use auxtools,only:frac_to_deg_min_sec
  use generic_info
  use get_ast_data,only:get_astro_data
  use intercept,only:wgs84_intercept
  use bore_site
  implicit none
  include 'astmath.cmn'
  integer::nargs,status,lh,i,ih,alen,iq,ir
  character (len=clen)::filename,hold
  logical::iex
  type(header)::finfo
  real (kind=8)::xp,yp,LOD,dat,dut1,ut1,tut1,jdut1,utc,tai,ttt,jdtt,tdb,ttdb,&
       jdtdb,ddpsi,ddeps,deltapsi,trueeps,meaneps,omega,frac,tt
  integer::gps_seconds_2009Jan01_00_00_00,ls_2012_06_30
  integer::lun,iolen,imult
  integer::year,mon,day,hour,min,lines
  real*8::sec
  real (kind=8),dimension(:),allocatable::ctmp ! used for temporary values
  real (kind=8),dimension(512)::theta,view_zen,view_az
  real (kind=8),dimension(15)::pvq_line
  real (kind=8),dimension(3)::r_iss,omegaearth,r_ecef
  real (kind=8),dimension(4)::q_iss
  real (kind=8),dimension(3,3)::rot_body,r_hico_to_iss,prec,nut,st,stdot,&
       t_hico_to_ecef,t_eci_to_ecef,pm,r_hico_bs
  real (kind=8), parameter::theta_stowed=-60.d0
  real (kind=8),dimension(3,512)::svec,llh,v_ecef
  character(len=12)::error
  character (len=clen)::hdr_name
  integer::ilen,err
  real (kind=8),dimension(2,512)::sol_zen_az
!  character (len=clen)::ast_data_dir
  character (len=clen)::leapsec_file
  character (len=clen)::earth_orient_file
  real (kind=8),dimension(3)::bsx
! header items
  integer,dimension(3)::icd
  real,dimension(3)::ict,lat,lon,ivz,iva
  character (len=1)::lat_hem,lon_hem
  character (len=5),parameter::code_version='2.0.0'
  character (len=11),parameter::code_date='2012/Nov/06'
  character (len=56),parameter::code_author='Marcos Montes (marcos.montes@nrl.navy.mil, 202-767-7308)'
  character (len=9),parameter::code_name='proc_hico'
  character (len=8)::grun_date
  character (len=10)::grun_time
  character (len=5)::grun_zone
  character (len=45)::username,computername
  character (len=clen),dimension(-1:37)::history
  real (kind=8),dimension(3,1)::ssllh ! subsat lat lon and iss height
  character (len=10)::direction=''
  character (len=clen)::junk
!  external convtime
  interface convtime
     SUBROUTINE CONVTIME    ( Year, Mon, Day, Hr, minute, SEC, &
          TimeZone, TypeUTIn, DUT1, DAT, xp, yp, &
          UT1, TUT1, JDUT1, UTC, TAI, TT, TTT, & 
          JDTT, TDB, TTDB, JDTDB, DDPSi, DDEps, &
          LOD, Error )
       INTEGER Year, Mon, Day, Hr, minute, TimeZone
       REAL*8 DUT1, DAT, xp, yp, UT1, TUT1, JDUT1, UTC, TAI, TT, TTT,&
            Sec, JDTT, TDB, TTDB, JDTDB, DDEPs, DDPsi, LOD, MFME, MJD
       CHARACTER TypeUTIn
       CHARACTER*(*)  Error
     end SUBROUTINE CONVTIME
  end interface convtime
  real (kind=8),parameter::ff=1.d0/298.257223563d0 !WGS-84
  real (kind=8),parameter::re=6378137d0    ! WGS-84
  real (kind=8),parameter::lower=re + 2.5d5, upper=re + 5.d5 ! upper and lower bounds for ISS in m wrt WGS-84
  character (len=clen)::vardir

!!!!!!!!!!!!!!!!!!!!!!! Start !

! First, get the argument name. The first argument should be the "pos_vel_quat" file that is output by 
! interpolate_fields_to_hico_times.

!!!!!!! Begin processing command line arguments !!!!!!!!!!!!!!!
  error=''
  nargs=command_argument_count()
  if (nargs.eq.3.or.nargs.eq.6) then
  else
     print*,'wrong number of args!'
     print*,'either 3 args or 6'
     print*,'Usage: '
     print*,'geogen_hico /path/to/iss*pos_vel_quat.csv leapsec_file earth_orient_file'
     print*,'assumes bore sight offsets are all 0.d0'
     print*,' OR'
     print*,'geogen_hico /path/to/iss*pos_vel_quat.csv leapsec_file earth_orient_file 0.900000 0.1 -0.01'
     print*,'where the last three are XYZ bore sight offsets in degrees (if one is present, all three must be)'
     call exit(101)
  endif

  filename=''
  call get_command_argument(1,filename,status=status,length=alen)
  if (status.ne.0) then 
     print*,'Bad status for obtaining first (filename) argument; stopping'
     print*,'status = ',status
     print*,'filename = ',trim(filename)
     call exit(101)
  end if

!  ast_data_dir=''
  leapsec_file=''
  call get_command_argument(2,leapsec_file,status=status,length=alen)
  leapsec_file=trim(leapsec_file)
  if (status.ne.0) then 
     print*,'Bad status for obtaining 2nd (leapsec_file) argument; stopping'
     print*,'status = ',status
     print*,'leapsec_file = ',trim(leapsec_file)
     call exit(101)
  end if

  earth_orient_file=''
  call get_command_argument(3,earth_orient_file,status=status,length=alen)
  earth_orient_file=trim(earth_orient_file)
  if (status.ne.0) then 
     print*,'Bad status for obtaining 3rd (earth_orient_file) argument; stopping'
     print*,'status = ',status
     print*,'earth_orient_file = ',trim(earth_orient_file)
     call exit(101)
  end if

  inquire(file=trim(filename),exist=iex)
  if (iex.eqv..false.) then
     print*,'file does not exist: ',filename
     print*,'stopping'
     call exit(101)
  endif

  if (nargs.ne.6) then ! PROCESs bore sight parameters
     bsx=0.d0 ! radians or degrees, doesn't matter, they're not present!
  else
     junk=''
     call get_command_argument(4,junk,status=status,length=alen)
     read(junk,*)bsx(1)
     junk=''
     call get_command_argument(5,junk,status=status,length=alen)
     read(junk,*)bsx(2)
     junk=''
     call get_command_argument(6,junk,status=status,length=alen)
     read(junk,*)bsx(3)
     bsx=bsx*deg2rad
  end if
  call bore_sight(bsx(1),bsx(2),bsx(3),r_hico_bs)
  r_hico_bs=transpose(r_hico_bs)

!!!!!!!!! Done processing command line arguments !!!!!!!!!!!!!!!!!!

  call read_pos_vel_quat(filename,finfo,lines)
  allocate(ctmp(lines))

!!!!!! Begin checking magnitudes in pvq file !!!!!!!
! Check validity of pvq interpolated data POSITION magnitude; require 200 km above WGS-84 ellipsoid
! NOTE: this may still produce altitudes too low, but not bad for a first check
!  if (any(ieee_is_nan(pvq_data(2:4,:)))) then 
!     print*,'Retruning; NaN in at least one position in the input file'
!     stop
!  endif

  do iq=1,lines
     do ir=2,4
        if (is_nan(pvq_data(ir,iq))) then
           print*,'Returning; NaN detected in input position'
           print*,'PLEASE CHECK FILE: ',trim(finfo%name)
           print*,'And its source files: ',trim(finfo%source), trim(finfo%subset)
           print*,'STOPPING'
           call exit(110) 
        end if
     end do
     do ir=8,11
        if (is_nan(pvq_data(ir,iq))) then
           print*,'Returning; NaN detected in input quaternion'
           print*,'PLEASE CHECK FILE: ',trim(finfo%name)
           print*,'And its source files: ',trim(finfo%source), trim(finfo%subset)
           print*,'STOPPING'
           call exit(110) 
        end if
     end do
  end do

  ctmp=sqrt(sum((pvq_data(2:4,:)*0.3048)**2,dim=1))
  if (any(ctmp.lt.lower.or.ctmp.gt.upper)) then
     print*,'ERROR: POSITIONS ARE BELOW 200 km above WGS-84 or above 500 km above WGS-84!!!!'
     print*,'PLEASE CHECK FILE: ',trim(finfo%name)
     print*,'And its source files: ',trim(finfo%source), trim(finfo%subset)
     print*,'STOPPING'
     call exit(110) 
  endif

! Check valididty of interpolated velocity magnitude: want at least 6.5 km/s. Pretty much checks against 0s.  
  ctmp=sqrt(sum((pvq_data(5:7,:)*0.3048)**2,dim=1))
  if (any(ctmp.lt.6500.d0.or.ctmp.gt.9000.d0)) then
     print*,'ERROR: ISS velocity out of range!'
     print*,'PLEASE CHECK FILE: ',trim(finfo%name)
     print*,'And its source files: ',trim(finfo%source), trim(finfo%subset)
     print*,'STOPPING'
     call exit(110)
  endif

! Check magnitude of interpolated USGNC quanterions
  ctmp=sqrt(sum((pvq_data(8:11,:))**2,dim=1))
  if (any(abs(ctmp - 1.d0).gt.0.01)) then
     print*,'ERROR: Interpolated ISS USGNC quaternions are far from normalized!'
     print*,'PLEASE CHECK FILE: ',trim(finfo%name)
     print*,'And its source files: ',trim(finfo%source), trim(finfo%subset)
     print*,'STOPPING'
     call exit(110)
  endif
  deallocate(ctmp) ! currently ctmp only used for above error checks
!!!!!!!!! Done checking magnitudes in PVQ file !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now, the times in the HICO files are in seconds since 2009-Jan01-00:00:00 UTC.
! N.B. The original files have GPS times in them, and the file interpolate_fields_to_hico_times.pro
! calculates the GPS second for the above UTC date, and subtracts it form the GPS times in the files. 
! Now, I've verified that the times in the original CSV files are GPS. The YYYYMMDD hhmmss in the files
! are not UTC, but are instead just GPS time converted. So they differ from UTC by the number of leap
! seconds there are. So, while my EPOCH date is UTC, NO LEAP SECONDS are included since the epoch. Which means
! that after 2012/06/30 23:59:59, I have an added second (leap second) that I will need to account for. 
! So, I need to properly count seconds and properly account for leap seconds in order to really get/stay in UTC. 
  if (finfo%epoch.ne.'2009 Jan 01 00:00:00 UTC') then
     print*,'epoch = ',hold(:len_trim(hold))
     print*,'In file : ',trim(filename)
     print*,'Epoch not expected epoch'
     print*,'Stopping'
     call exit(110)
  endif
! GPS Week 0 began at 00:00:00 UTC 1980 Jan 06
! sec in day*( (days in 4yr) * 7 sets +
! yr2008 - first five days of 1980)
!          + leap seconds since 1980 Jan 6
  gps_seconds_2009Jan01_00_00_00 = 86400*((3 * 365 + 366)*7 + &
       366 - 5 ) + 15  
! All seconds in file need to be added to that.
! Seconds since 2009/0/0 until leap second at end of 2012/06/30:
! all of regular years 2009, 2010, 2011
  ls_2012_06_30 = 86400*((365*3)+(31+29+31+30+31+30))+1

! Read nutation data for nutation calculations
  call get_environment_variable("OCVARROOT", vardir)
  call initReduc(trim(vardir)//'/hico/nutation.dat')
  
! xp, yp, LOD are in a file that can be automatically downloaded from the USNO. 
  xp=0.d0  ! Not needed for 100 m pixels; can be read from file
  yp=0.d0  ! Not needed for 100 m pixels; can be read from file 
  LOD=0.d0 ! we don't use this
! dat, dut1 can also be obtained from file
  dat=34.d0
  dut1=-0.4984030d0  

  call increment_seconds_to_date(year,mon,day,hour,min,sec,pvq_data(1,1))

! Get some USNO provided parameters NOTE: USNO FILE NEEDS TO BE UPDATED PERIODICALLY
  call get_astro_data(trim(leapsec_file),trim(earth_orient_file),year,mon,day,xp,yp,LOD,dut1,dat,err)
  if (err.ne.0) then
     print*,'Error: could not read astro data for deltaUT1, etc'
!     print*,'Directory searched was:',trim(ast_data_dir)
     print*,'Stopping'
     call exit(101)
  end if
  
! Define the rotation from HICO frame of reference to ISS frame of reference
! NOTE: This assumes perfect alignment. 
  r_hico_to_iss(:,:)= 0.d0
  r_hico_to_iss(1,1)=-1.d0
  r_hico_to_iss(2,2)=-1.d0
  r_hico_to_iss(3,3)= 1.d0
! include bore sight offsets
  r_hico_to_iss=matmul(r_hico_to_iss,r_hico_bs)

!!!!! BEGIN HICO pointing angles:
! I believe the 1/2 FOV is 3.461515 degrees (from Mike C?)
! A document by Bob Lucke states that the relationship is as follows;
! this was in MJM's HICO email folder in an email from Bill Snyder.
! UNITS ARE DEGREES

  if (trim(finfo%iss_orientation).eq.'+XVV') then
     imult =  1
  else
     imult = -1
  endif
  do i=1,512
     frac= imult*(i-256.5d0)/(255.5d0)
     theta(i) = ( -3.487d0  + 0.035d0 * frac**2)*frac
  end do
  theta= (theta + theta_stowed +finfo%theta_from_stowed)*deg2rad

! Now theta is in radians and accounts for angle from stowed position
! Now I need to convert to pointing angles in HICO's frame of reference.
! Note that HICO is in the YZ plane of HICO. The stowed angle (+ is right
! hand rule about X, measured from +Z) is in the +Y+Z quadrant, even 
! for negative values of  theta...
  svec(1,:) =  0.d0
  svec(2,:) = -sin(theta)
  svec(3,:) =  cos(theta)

!!!! End of HICO pointing angles

! Determine output record length
  inquire(iolength=iolen)llh(1:2,:),sol_zen_az,view_zen,view_az

! Open output file
  open(newunit=lun,file=trim(finfo%outname),action='write',access='direct',&
       form='unformatted',recl=iolen,status='unknown')

! Polar parameters
  call polarm(xp,yp,pm) ! parameters read from file above; not time calculation

  do i=1,lines

     pvq_line=pvq_data(:,i)

     r_iss=pvq_line(2:4)*0.3048d0 ! convert FEET to METERS

     q_iss=[pvq_line(11),pvq_line(8),pvq_line(9),pvq_line(10)]

     call q_to_r(q_iss,rot_body) ! make rotation matrix from USGNC quaternions

     call increment_seconds_to_date(year,mon,day,hour,min,sec,pvq_line(1))

     call convtime(year,mon,day,hour,min,sec,0,'0',dut1,dat,xp,yp,&
          ut1,tut1,jdut1,utc,tai,tt,ttt,jdtt,tdb,ttdb,jdtdb,ddpsi,ddeps,LOD,error)

     call precession(ttdb,prec)

     call nutation(ttdb, deltapsi,trueeps,meaneps,omega,nut)

     call sidereal(jdut1,deltapsi,trueeps,omega,lod,st,stdot,OmegaEarth,2)

     ! from here, call a routine that just makes the transformation matrices
     t_eci_to_ecef=matmul(pm,matmul(st,matmul(nut,prec))) ! for positions
     t_hico_to_ecef=matmul(t_eci_to_ecef,matmul(rot_body,r_hico_to_iss)) ! for attitudes from HICO frame to ECEF

     ! The position does not depend on the attitude transformation
     r_ecef=matmul(t_eci_to_ecef,r_iss)

     ! Here, v_ecef are the 512 pointing vectors. 
     ! v_ecef=matmul(t_hico_to_ecef,v_hico) ! v_hico=[3,NS]
     v_ecef=matmul(t_hico_to_ecef,svec)

     ! use R Reynolds formula for distance along ray
     ! calls to convert ECEF XYZ to ECEF LLH
     ! Note: all math involving Earth surface is done there, including view geometry
     call wgs84_intercept(r_ecef,v_ecef,llh,view_zen,view_az)
     
! and solar zenith and azimuth (note, for solar calc, input long has W>0, E<0 which is why there is neg for longitude arg)
     sol_zen_az=rad2deg*solar_geometry(year,mon,day,hour,min,sec,llh(2,:),-llh(1,:))
! to write BIL I need to transpose the "bip-like" variables here....
     write(lun,rec=i)transpose(llh(1:2,:)),view_zen,modulo(360.d0+view_az,360.d0),transpose(sol_zen_az)

! Below are some items for output header, grabbed near center of line
     if (i.eq.lines/2) then 
        icd=[year,mon,day]
        ict=[real(hour),real(min),real(sec)]
        call frac_to_deg_min_sec(real(abs(llh(2,255))),0.,0.,lat(1),lat(2),lat(3))
        call frac_to_deg_min_sec(real(abs(llh(1,255))),0.,0.,lon(1),lon(2),lon(3))
        if (llh(1,255).ge.0) then 
           lon_hem='E'
        else
           lon_hem='W'
        end if
        if (llh(2,255).ge.0) then 
           lat_hem='N'
        else
           lat_hem='S'
        end if
        ivz=[real((view_zen(255)+view_zen(256)))/2,0.,0.]
        iva=[real(view_az(255)),0.,0.]
        call frac_to_deg_min_sec(real(view_zen(255)),0.,0.,ivz(1),ivz(2),ivz(3))
        call frac_to_deg_min_sec(modulo(360.+real(view_az(255)),360.),0.,0.,iva(1),iva(2),iva(3))
        call ecef2latlon(spread(r_iss,dim=2,ncopies=1),ssllh)
        if (pvq_line(7).ge.0) then ! a bit crude, based on ECI velocity
           direction='ascending'
        else
           direction='descending'
        endif
     endif

  enddo ! loop over lines
  close(lun) ! output file

! write header (history, etc....)
  call date_and_time(grun_date,grun_time,grun_zone)
  call get_username_s(username)
  call get_hostname_s(computername)
  ilen=len_trim(finfo%outname)
  hdr_name=finfo%outname(:ilen-4)//'.hdr'
  history(-1)='history = begins'
  history(0) =''
  history(1) ='geometry_code_version = '//code_version 
  history(2) ='geometry_code_date = '//code_date
  history(3) ='geometry_code_author = '//code_author
  history(4) ='geometry_code_name = '//code_name
  history(5) ='geometry_code_SPosVelQuatCsv = '//trim(finfo%name)
  history(6) ='geometry_source_csv_filename = '//trim(finfo%source)
  history(7) ='geometry_subset_csv_filename = '//trim(finfo%subset)
  history(8)='geometry_run_date = '//grun_date
  history(9)='geometry_run_time = '//grun_time//' '//grun_zone
  history(10)='geometry_run_on = '//trim(computername)
  history(11)='geometry_run_by = '//trim(username)
  history(12) ='pvq_code_name = '//trim(finfo%code_name)
  history(13) = 'pvq_code_version = '//trim(finfo%code_version)
  history(14) = 'pvq_code_date = '//trim(finfo%code_date)
  history(15) = 'pvq_code_author = '//trim(finfo%code_author)
  history(16) = 'pvq_code_computer = '//trim(finfo%code_runon)
  history(17) = 'pvq_code_username = '//trim(finfo%code_runby)
  history(18) = 'pvq_code_IDL_family = '//trim(finfo%code_idl_family)
  history(19) = 'pvq_code_IDL_os = '//trim(finfo%code_idl_os)
  history(20) = 'pvq_code_IDL_osname = '//trim(finfo%code_idl_osname)
  history(21) = 'pvq_code_IDL_version = '//trim(finfo%code_idl_version)
  history(22) = 'pvq_code_IDL_start_time = '//trim(finfo%code_start_time)
  history(23) = 'pvq_code_IDL_end_time = '//trim(finfo%code_end_time)
  write(history(24),fmt='(a,f10.6)')'geometry_hico_angle_in_degrees_from_stowed = ',finfo%theta_from_stowed
  history(25) ='geometry_sensor_orientation = '//trim(finfo%iss_orientation)
  history(26) ='geometry_ISS_Z_direction = '//trim(direction)
  write(history(27),fmt='(a,f15.11)')'hico_exposure_interval_s = ',finfo%exposure_interval
  write(history(28),fmt='(a,f15.11)')'hico_trigger_pulse_width_s = ',finfo%trigger_pulse_width
  write(history(29),fmt='(a,f15.11)')'iss_ODRC_broadcast_time-gps_time_s = ',finfo%odrc_gps_seconds
  write(history(30),fmt='(a,f15.11)')'geometry_deltaUT1_s = ',dut1
  write(history(31),fmt='(a,f10.6)')'geometry_deltaT_s = ',dat
  write(history(32),fmt='(a,f15.11)')'geometry_LOD = ',LOD
  write(history(33),fmt='(a,f15.11)')'geometry_xp_rad = ',xp
  write(history(34),fmt='(a,f15.11)')'geometry_yp_rad = ',yp
  write(history(35),fmt='(a,2(f15.11,","),f15.11,a)')'geometry_bore_sight_params = { ',bsx,'}'
  history(36)=''
  history(37)='history = ends'
  
  open(newunit=lun,file=trim(hdr_name),access='sequential',action='write',status='replace')
  call write_envi_header(lun=lun,nsamples=512,nlines=lines,nbands=6,sort_order=2,&
       data_type=5,desc_string='HICO geometry file, with calculated positions, and solar & view geometry, on the ground',&
       band_names='longitude, latitude, view zenith, view azimuth, sol zenith, '//&
       'sol azimuth',date=icd,time=ict,lat=lat,long=lon,&
       lat_hem=lat_hem,long_hem=lon_hem,zenith=ivz,azimuth=iva,&
       units=[' degrees ',' degrees ',' degrees ',' degrees ',' degrees ',' degrees '],&
       history = history,altitude=real(ssllh(3,1)/1.e3),sensor_type='HICO-ISS')
  close(lun) ! header file
  
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  subroutine wgs84_intercept(rsc,u,out,view_zen,view_az)
!!    use ecef2llh
!!    ! Thanks to beautiful form of equations from Reid Reynolds.
!!    implicit none
!!    real (kind=8),intent(in),dimension(3)::rsc ! spacecraft radius
!!    real (kind=8),intent(in),dimension(:,:)::u ! unit vector pointing toward earth, a lot of them
!!    real (kind=8),intent(out),dimension(size(u,1),size(u,2))::out !lon,lat,height (m)
!!    real (kind=8),intent(out),dimension(size(u,2)),optional::view_zen,view_az
!!    ! Parameters
!!    real (kind=8),parameter::ff=1.d0/298.257223563d0 !WGS-84
!!    real (kind=8),parameter::re=6378137d0    ! WGS-84
!!    ! F is really the diagonal elements of a 3x3 matrix, all rest are 0
!!    ! but we only need the diagonal for the WGS-84 ellipsoid
!!    real (kind=8),parameter,dimension(3)::F=[1.d0,1.d0,1.d0/(1.d0-ff)**2]
!!    !more parameters
!!    real (kind=8), parameter:: Pi      =  4.0D0 * DATAN(1.0D0)
!!    real (kind=8), parameter:: HalfPi  =  0.5D0*pi
!!    real (kind=8), parameter:: TwoPi   =  2.0D0*pi
!!    real (kind=8), parameter:: Rad2Deg = 180.0D0/pi
!!    real (kind=8), parameter:: Deg2Rad = pi/180.0D0
!!    ! Locals
!!    real (kind=8)::c
!!    real (kind=8),dimension(size(u,2))::a,b,s,uDotE,uDotN,det
!!    real (kind=8),dimension(size(u,1),size(u,2))::rout
!!    real (kind=8),dimension(size(u,1),size(u,2))::snorm
!!
!!! If Rg is a vector from the center to the surface of the earth, then
!!! the equation of the ellipsoid is transpose(Rg)F*Rg=Re*Re
!!
!!! The equation of a ray with unit direction u from the spacecraft to the earth is:
!!!  Rg = Rsc + s*u ; s is length of ray; substitue this equation into the one above, 
!!! and use the form of the various vectors and matrices to get
!!! the equations below, which yields a simple quadratic equation for s.  
!!! For us, Rsc is one location, and we have 512 
!!! directions (elements of u). 
!!    
!!! precalculate a few things, using the magic of matmul's rules
!!    c=sum(F*rsc**2)-re**2 ! scalar
!!    b=2.d0*matmul(F*Rsc,u) ! ns elements
!!    a=matmul(F,u**2)  ! ns elements
!!
!!    det=b**2 - 4*a*c
!!    if (any(det.lt.0)) then
!!       print*,'ERROR IN WGS84_intercept: invalid answer. no intercept'
!!       print*,'CHECK INPUT!'
!!       print*,minval(det),minloc(det)
!!       print*,maxval(det),maxloc(det)
!!       stop
!!    endif
!!
!!! Now the distance along the ray, from the spacecraft is:   
!!    s= (-b -sqrt(det))/(2.d0 * a ) ! ns elements
!!
!!! Once you know s, rout is the the location in ECEF of the intercept with the Earth,
!!! traveling a distance s from rsc in the direction(s) u. 
!!    rout= spread(rsc,dim=2,ncopies=512) + spread(s,dim=1,ncopies=3)*u ! 3xns
!!
!!! Now we need an ECEF->LATLONH conversion
!!    call ecef2latlon(rout,out) 
!!
!!! The below essentially use terms from the ECEF -> ENU conversion on the surface of an oblate spheroid,
!!! our WGS-84 ellipsoid
!!! The normal to the ellipsoid has the normal snorm
!!    if (present(view_zen)) then 
!!       snorm(1,:)=cos(llh(1,:))*cos(llh(2,:))
!!       snorm(2,:)=sin(llh(1,:))*cos(llh(2,:))
!!       snorm(3,:)=sin(llh(2,:))
!!! The cos(view zenith) is formed by the pointing vector from the ground to the spacecraft
!!! dotted into the surface normal; this is (for one location) sum(-u(:,i)*snorm(:,i)) 
!!       view_zen=acos(sum(-u*snorm,dim=1))*rad2deg ! deg from zenith
!!    endif
!!
!!    if (present(view_az)) then 
!!       uDotE=sin(llh(1,:))*u(1,:) - cos(llh(1,:))*u(2,:)
!!       uDotN=cos(llh(1,:))*sin(llh(2,:))*u(1,:) + sin(llh(1,:))*sin(llh(2,:))*u(2,:) - &
!!            cos(llh(2,:))*u(3,:)
!!       view_az=atan(udotE,uDotN)*rad2deg ! deg from N, clockwise
!!    end if
!!    out=out*rad2deg
!!
!!  end subroutine wgs84_intercept
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine increment_seconds_to_date(yyyy,mm,dd,hour,min,sec,delta)
! MJM 2012/10/24 Adjusting the time from the epoch (2009-01-01-00:00:00 UTC).
! Input is GPS seconds since the epoch. Output is UTC date/time. 
    implicit none
    integer, intent(out)::yyyy,mm,dd,hour,min
    real (kind=8),intent(out)::sec
    real (kind=8),intent(in)::delta
    integer::delta_d
    integer,parameter,dimension(12)::mlen=[31,28,31,30,31,30,31,31,30,31,30,31]
    integer::idinc,mmlen,nd
    real (kind=8)::dinc

    dinc=delta
    idinc=floor(delta)
    ! Seconds since epoch of the (only) leap second so far
    ls_2012_06_30 = 86400*((365*3)+(31+29+31+30+31+30))+1 

    if (idinc.gt.ls_2012_06_30) then ! date after leap second, so adjust this to new epoch
       yyyy=2012
       mm=07
       dd=01
       idinc=idinc-ls_2012_06_30
       dinc=dinc-ls_2012_06_30
       hour=0
       min=0
       sec=0.d0
    else ! else use epoch definition as start
       yyyy=2009
       mm=01
       dd=01
       hour=0
       min=0
       sec=0.d0
    endif

! Subtract off full months    
    do
       if (mm.eq.2) then
          if (is_leap_year(yyyy).eq.1) then 
             mmlen=mlen(mm)+1
          else
             mmlen=mlen(mm)
          endif
       else
          mmlen=mlen(mm)
       endif
       if (idinc.lt.86400*mmlen) exit ! cannot advance by full months
       idinc=idinc-86400*mmlen
       dinc=dinc - 86400*mmlen
       if (mm.ne.12) then
          mm=mm+1
       else
          mm=1
          yyyy=yyyy+1
       endif
    end do
! at this point, we have less a full month. Let's figure out how many whole days
    nd=(idinc/86400)
    dd=dd+nd
    idinc=idinc - 86400*nd
    dinc=dinc - 86400*nd
! Less than 24 hours left
    hour=(idinc/3600)
    idinc=idinc - hour*3600
    dinc=dinc - hour*3600
! Less than 1 hour left    
    min=(idinc/60)
    idinc=idinc - min*60
    dinc=dinc - min*60
! Less than 1 minute left
    sec=dinc
    
  end subroutine increment_seconds_to_date
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure integer function is_leap_year(yyyy)
    implicit none
    integer,intent(in)::yyyy
    integer::m4,m400,m100,leap
    m4=modulo(yyyy,4)
    if (m4.eq.0) then
       m100=modulo(yyyy,100)
       if (m100.ne.0) then
          leap=1
       else !m100==0
          m400=modulo(yyyy,400)
          if (m400.eq.0) then 
             leap=1
          else
             leap=0
          endif
       endif
    else
       leap=0
    endif
    is_leap_year=leap
  end function is_leap_year
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
end program proc_hico
