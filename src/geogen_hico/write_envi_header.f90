module envi_header_mod
  
  ! 2008-08-27 Added 'image_data_quantity' kwd
  ! 2004-02-05 Removed all references to 'ATREM'
  !            Removed most references to 'ENVI'
  ! 2004-Aug-31 MJM Added name_len to various key structures
  ! 2005-Nov-29 MJM Added use of wl_units, z_plot_labels, and 
  !                 mapping parameters.
  use envi_header
  use check_range_mod
  use default_character_length,only:charlen
  implicit none
  private
  character (len=21),parameter::file_name='write_envi_header.f90'
  character (len=15),parameter::module_name='envi_header_mod'
  integer,parameter::BSQ=0, BIP=1, BIL=2
  public:: write_envi_header, read_envi_header,newread_envi_header
  public:: read_envi_header_nbarrays

  !Notes:
  !    Byte-order: 0=Least Significant Byte First (DEC & MSDOS)
  !                1=Most Significant Byte First (all others)
  !    Data Type:  1=byte, 2=2 byte integer, 3=4byte integer,
  !                4=4byte real, 5=8byte real, 6=8byte complex
  !                Note: **ALL** data types are SIGNED
  type int_key
     character (len=charlen)::name
     integer, pointer::value
     integer::name_len
  end type int_key
  integer, parameter::n_wanted_int_keys=8
  type (int_key), dimension(n_wanted_int_keys)::wanted_int_keys
  
  type rel_key
     character (len=charlen)::name
     real, pointer::value
     integer::name_len
  end type rel_key
  integer, parameter::n_wanted_rel_keys=1
  type (rel_key),dimension(n_wanted_rel_keys)::wanted_rel_keys
  
  type str_key
     character (len=charlen)::name
     character (len=charlen), pointer::value
     integer::vlen !character length of 'value'
     integer::name_len
  end type str_key
  integer, parameter::n_wanted_str_keys=5
  type (str_key), dimension(n_wanted_str_keys)::wanted_str_keys
  
  type fltarr_key
     character (len=charlen)::name
     real, pointer, dimension(:):: value
     integer::name_len
  end type fltarr_key
  integer, parameter::n_wanted_aNB_keys=3
  type (fltarr_key), dimension(n_wanted_anb_keys)::wanted_anb_keys
  integer, parameter::n_wanted_ar3_keys=5
  type (fltarr_key), dimension(n_wanted_ar3_keys)::wanted_ar3_keys
  
  type intarr_key
     character (len=charlen)::name
     integer, pointer, dimension(:)::value
     integer::name_len
  end type intarr_key
  integer, parameter::n_wanted_iar_keys=2
  type (intarr_key), dimension(n_wanted_iar_keys)::wanted_iar_keys
  integer, parameter::n_wanted_inb_keys=1
  type (intarr_key), dimension(n_wanted_iNB_keys)::wanted_inb_keys
  
  type charr_key
     character (len=charlen)::name
     character (len=charlen), pointer, dimension(:)::value
     integer::name_len
  end type charr_key
  integer, parameter::n_wanted_char_keys=2
  type (charr_key), dimension(n_wanted_char_keys)::wanted_char_keys
  integer,parameter::n_wanted_map_keys=4
  type (charr_key),dimension(n_wanted_map_keys)::wanted_map_keys

  ! Mapping parameters
  character (len=charlen),dimension(:),target,allocatable::map_info,pixel_size,&
       geo_points,projection_info
  
  
contains
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine write_envi_header(lun, nsamples, nlines, nbands, sort_order, &
       & data_type, desc_string, time, date, lat, lat_hem, long, long_hem,&
       & zenith, azimuth, altitude,imscale,units,&
       & wavobs, fwhm,bbl, band_names,history,sensor_type,default_stretch,&
       & x_start, y_start,default_bands, wl_units, z_plot_titles,image_data_quantity)
    ! ?? Added FWHM and vavelength headers
    !1999-02-08 Added band names to header. -MJM 
    !1999-09-20 Added various optional headers-MJM
    !1999-09-24 Added optional history section - MJM
    !1999-10-19 Added sensor_type, default_stretch. -MJM
    !2000-02-23 Fixed write for wavobs, fwhm, and bbl when there is only
    !           a single item on the last line. -MJM
    !2000-04-13 Added x_start & y_start. -MJM
    !2002-04-16 Fixed bug in line to write multiple units; only failed on Sun.-MJM
    !2002-08-01 Added support for default bands on output. Default bands is 
    !           required to have three elements. If they are all identical, 
    !           then the one-element version is printed out; otherwise a 
    !           three element version is written.
    !2002-08-09 Fixed bug in write of bbl since I changed bbl to more compact format.
    !           ifmt needs I2, not I1.
    !2003-05-28 Fixed formatting for bbl when < 24 bands
    !2005-11-29 Added wl_units, z_plot_titles (passed in) as well as just passing
    !           through a variety of mapping keyword/values the ENVI may recognize.
    !2008-08-29 Added image_data_quantity

    use file_io_mod,only:byteorder

    implicit none

    integer, intent(in)::lun, nsamples, nlines, nbands,sort_order, data_type
    character (len=*), intent(in)::desc_string
    integer, intent(in),dimension(3),optional::date,default_bands
    real, intent(in), dimension(3), optional::time,lat,long,zenith,azimuth
    character (len=1),intent(in),optional::lat_hem,long_hem
    real,intent(in),optional::altitude
    real,dimension(:), intent(in),optional::imscale
    character (len=*), dimension(:),intent(in),optional::units
    real, dimension(:), intent(in), optional::wavobs,fwhm
    integer, dimension(:), intent(in), optional::bbl
    character (len=*), intent(in), optional::band_names
    character (len=*), dimension(:), intent(in), optional::history
    character (len=*), intent(in), optional::sensor_type
    character (len=*), intent(in), optional::default_stretch
    integer, intent(in), optional::x_start, y_start
    character (len=*), intent(in), optional::wl_units
    character (len=*), intent(in),dimension(2),optional::z_plot_titles
    character (len=*), intent(in),optional:: image_data_quantity

    integer::clen,lholdm,i,nbhold,nbdiff,band_name_size,ngroup,nleft
    integer::bhold,nelem,n_his,l_his,sen_len,defstr_len
    integer::nbbl1,nbbl2,zpl1,zpl2,dim
    integer,allocatable,dimension(:)::mlens
    character (len=15)::holdm,aholdm
    character (len=40)::ifmt
    character (len=1)::extra
    !

    clen=len_trim(desc_string)
    write(lun,'(A)')'ENVI'
    write(lun,'(A)')'description = {'
    write(lun,'(A)')'  '//desc_string(:clen)//' }'
    write(holdm,'(I5)')nsamples
    aholdm=adjustl(holdm)
    lholdm=len_trim(aholdm)
    write(lun,'(A)')'samples = '//aholdm(1:lholdm)
    holdm='               '
    write(holdm,'(I5)')nlines
    aholdm=adjustl(holdm)
    lholdm=len_trim(aholdm)
    write(lun,'(A)')'lines   = '//aholdm(1:lholdm)
    holdm='               '
    write(holdm,'(I5)')nbands
    aholdm=adjustl(holdm)
    lholdm=len_trim(aholdm)
    write(lun,'(A)')'bands   = '//aholdm(1:lholdm)
    write(lun,'(A)')'header offset = 0'  ! No header in output data files
    write(lun,'(A)')'file type = ENVI Standard'
    write(lun,'(A,I1)')'data type = ',data_type  ! 1 for land mask, 2 otherwise
    if (sort_order.eq.bsq) then
       write(lun,'(A)')'interleave = bsq'
    else if (sort_order.eq.bip) then 
       write(lun,'(A)')'interleave = bip'
    else if (sort_order.eq.bil) then
       write(lun,'(A)')'interleave = bil'
    endif
    !MS F90 seems to barf unless I protect the variable by assigning 
    !byteorder() to a variable:
    bhold=byteorder()
    write(lun,'(A,I1)')'byte order = ',bhold ! see byteorder()

    !Various optional, but extremely useful quantities.
    if (present(sensor_type)) then
       sen_len=len_trim(sensor_type)
       write(lun,'(A)')'sensor type = '//sensor_type(:sen_len)
    else
       write(lun,'(A)')'sensor type = Unknown'
    end if
    if (present(x_start)) then
       holdm=''
       write(holdm,'(I6)')x_start
       aholdm=adjustl(holdm)
       lholdm=len_trim(aholdm)
       write(lun,'(A)')'x start = '//aholdm(1:lholdm)
    end if
    if (present(y_start)) then
       holdm=''
       write(holdm,'(I6)')y_start
       aholdm=adjustl(holdm)
       lholdm=len_trim(aholdm)
       write(lun,'(A)')'y start = '//aholdm(1:lholdm)
    end if
    if (present(default_stretch)) then
       defstr_len=len_trim(default_stretch)
       write(lun,'(A)')'default stretch = '//default_stretch(:defstr_len)
    end if

    if (present(default_bands)) then
       if (all(default_bands.eq.default_bands(1))) then ! only one unique value
          write(lun,'(A,I4.1,A)') &
               &'default bands = {',default_bands(1),' }'
       else ! three values
          write(lun,'(A,I4.1,",",I4.1,",",I4.1,A)') &
               &'default bands = {',default_bands,'}'
       endif
    endif
    if (present(date)) write(lun,'(A,I5.4,",",I3.2,",",I3.2,A)') &
         &'image_center_date = {',date,'}'
    if (present(time)) write(lun,'(A,F4.0,",",F4.0,",",f7.3,A)')&
         &'image_center_time = {',time,'}'
    if (present(lat)) write(lun,'(A,F4.0,",",F4.0,",",f7.3,A)')&
         &'image_center_lat = {',lat,'}'
    if (present(lat_hem)) write(lun,'(A)')'image_center_lat_hem = '//lat_hem
    if (present(long)) write(lun,'(A,F5.0,",",F4.0,",",f7.3,A)')&
         &'image_center_long = {',long,'}'
    if (present(long_hem)) write(lun,'(A)')'image_center_long_hem = '//long_hem
    if (present(zenith)) write(lun,'(A,F4.0,",",F4.0,",",f7.3,A)')&
         &'image_center_zenith_ang = {',zenith,'}'
    if (present(azimuth)) write(lun,'(A,F5.0,",",F4.0,",",f7.3,A)')&
         &'image_center_azimuth_ang = {',azimuth,'}'
    if (present(altitude)) write(lun,'(A,F9.3)')'sensor_altitude = ',altitude
    if (present(imscale)) then
       nelem=size(imscale) ! Must be 1 or nbands
       if (nelem==1) then
          write(lun,'(A,F12.4,A)')'image_scale_factor = { ',imscale(1),'}'
       else !there must be nbands of them present
          nbhold=nbands/8  ! how many groups of 8
          if (mod(nbands,8).eq.0) nbhold=nbhold-1 !always a last line, so subtract 1
          nbdiff=nbands-8*nbhold ! There will be nbdiff elements in the last line
          ifmt=''
          write(ifmt,'(A,I1,A,A)')'(',nbdiff-1,'(f12.4,",")',',f12.4,"}")'
          write(lun,'(A)')'image_scale_factor = {'
          if (nbhold.gt.0) write(lun,'(8(f12.4,","))')(imscale(i),i=1,8*nbhold)
          write(lun,ifmt)(imscale(i),i=8*nbhold+1,nbands)
       end if
    end if
    if (present(units)) then
       nelem=size(units) ! Must be 1 or nbands
       nbhold=len_trim(units(1))
       if (nelem==1) then
          write(lun,'(A)')'image_unscaled_units = { '//units(1)(:nbhold)//'}'
       else !must be nbands elements
          !is there anyway to predict the sizes? Probably not.
          !Also, probably if here, there are not too many bands, so we'll
          !try this for now:
          ifmt=''
          write(ifmt,'(A,I4,A)')'(A,',nelem-1,'(A,","),A,"}")'
          write(lun,ifmt)'image_unscaled_units = { ',(units(i)(:len_trim(units(i))),i=1,nelem)
       end if
    endif

    if (present(image_data_quantity)) then
       write(lun,'(2a)')'image_data_quantity = ',image_data_quantity
    end if

    if (present(z_plot_titles)) then
       if (z_plot_titles(1).ne.'') then
          zpl1=len_trim(z_plot_titles(1))
          zpl2=len_trim(z_plot_titles(2))
          write(lun,'(5a)')'z plot titles = { ',z_plot_titles(1)(:zpl1),', ',&
               z_plot_titles(2)(:zpl2),'}'
       endif
    endif

    if (allocated(map_info)) then ! "passed through" if in input header file
       dim=size(map_info)
       if (allocated(mlens)) deallocate(mlens)
       allocate(mlens(dim))
       mlens=len_trim(map_info)
       ifmt=''
       write(ifmt,'(A,I2,A)')'(A,',dim-1,'(A,", "),A,"}")'
       write(lun,ifmt)'map info = { ',(map_info(i)(:mlens(i)),i=1,dim)
    endif

    if (allocated(projection_info)) then  ! "passed through" if in input header file
       dim=size(projection_info)
       if (allocated(mlens)) deallocate(mlens)
       allocate(mlens(dim))
       mlens=len_trim(projection_info)
       ifmt=''
       write(ifmt,'(A,I2,A)')'(A,',dim-1,'(A,", "),A,"}")'
       write(lun,ifmt)'projection info = { ',(projection_info(i)(:mlens(i)),i=1,dim)
    endif

    if (allocated(pixel_size)) then   ! "passed through" if in input header file
       dim=size(pixel_size)
       if (allocated(mlens)) deallocate(mlens)
       allocate(mlens(dim))
       mlens=len_trim(pixel_size)
       write(ifmt,'(A,I2,A)')'(A,',dim-1,'(A,", "),A,"}")'
       write(lun,ifmt)'pixel size = { ',(pixel_size(i)(:mlens(i)),i=1,dim)
    endif

    if (allocated(geo_points)) then
       dim=size(geo_points) 
       if (allocated(mlens)) deallocate(mlens)
       allocate(mlens(dim))
       mlens=len_trim(geo_points) 
       ! This is a bit different. There should be 4N elements, N an integer. 
       ! Easier to read if put 4 to a line
       write(lun,'(A)')'geo points = { '
       nbhold=(dim/4)-1 ! how many groups of 4, minus 1 (last line)
       ! All but last line:
       write(lun,'(4(" ",A,","))')(geo_points(i)(:mlens(i)),i=1,4*nbhold)
       ! Last line:
       write(lun,'(3(A,", "),A,"}")')(geo_points(i)(:mlens(i)),i=4*nbhold+1,dim)
    endif

    if (present(wl_units)) then
       if (wl_units.ne.'') then
          write(lun,'(2a)')'wavelength units = ',wl_units
       endif
    endif

    !     Now let's try to write the wavelengths and FWHM if they're 
    !     supplied. They're optional so we can write the H_20 and/or
    !     aerosol ouput file headers, too. 

    !     Observed wavelengths

    if (present(wavobs).or.present(fwhm)) then
       !     a specially constructed format statement for the last line
       nbhold=nbands/8  ! how many groups of 8
       if (mod(nbands,8).eq.0) nbhold=nbhold-1 !always a last line, so subtract 1
       nbdiff=nbands-8*nbhold ! There will be nbdiff elements in the last line
       ifmt=''
       if (nbdiff.gt.1) then
          write(ifmt,'(A,I1,A,A)')'(',nbdiff-1,'(f8.5,",")',',f8.5,"}")'
       else ! bad to have 0 elements when nbdiff-1=0, so special case
          ifmt='(f8.5,"}")'
       end if
    endif
    if (present(wavobs)) then
       write(lun,'(A)')'wavelength = {'
       write(lun,'(8(f8.5,","))')(wavobs(i),i=1,8*nbhold)
       write(lun,ifmt)(wavobs(i),i=8*nbhold+1,nbands)
    endif
    !     Now FWHM
    if (present(fwhm)) then
       write(lun,'(A)')'fwhm = {'
       write(lun,'(8(f8.5,","))')(fwhm(i),i=1,8*nbhold)
       write(lun,ifmt)(fwhm(i),i=8*nbhold+1,nbands)
    endif
    ! Now bad bands list
    if (present(bbl)) then 
       write(lun,'(A)')'bbl = {'
       nbbl1=nbands/24 ! number of lines with 24 elements
       if (mod(nbands,24).eq.0) nbbl1=nbbl1-1 ! always a last line
       nbbl2=nbands-24*nbbl1 ! number of lines with < 24 elements
       if (nbbl1 > 0) then !2003 May 28: only do this if >24 bands
          write(lun,'(24(I2,","))')(bbl(i),i=1,24*nbbl1)
       endif
       if (nbbl2.gt.1) then
          write(ifmt,'(A,I2,A,A)')'(',nbbl2-1,'(I2,",")',',I2,"}")'
       else
          ifmt='(I2,"}")'
       end if
       write(lun,ifmt)(bbl(i),i=24*nbbl1+1,nbands)
    end if
    !Band names : band_names needs to be formatted by calling program
    if (present(band_names)) then
       write(lun,'(A)')'band names = { '
       !comment: the band_names needed for the reflectance file seem to be too long
       !Break at the end of word before so line <=78 characters long
       band_name_size=len_trim(band_names)
       if (band_name_size.lt.300) then !short enough
          write(lun,'(A)')band_names//'}'
       else
          nleft=1
          !A scheme to break long strings at the end of words with max. line length 
          ! of 78 characters- don't put newlines in band_names
          nleft=1
          extra=''
          do ; if (nleft.ge.band_name_size) exit
             ngroup=index(band_names(nleft:nleft+min(77,band_name_size-nleft)),&
                  &', ',back=.true.)
             if (ngroup.eq.0.or.band_name_size-nleft.le.77) then
                ngroup=band_name_size-nleft+1
                extra='}'
             endif
             write(lun,'(a)')band_names(nleft:nleft+ngroup-1)//extra
             nleft=nleft+ngroup
          enddo
       endif
    end if
    ! Finally, the history,if present, gets written last
    if (present(history)) then
       n_his=size(history)
       write(lun,'(A)')''
       do i=1,n_his
          l_his=len_trim(history(i))
          write(lun,'(A)')history(i)(:l_his)
       end do
    endif
    !Done writing, so end
  end subroutine write_envi_header
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine newread_envi_header(found_keys, comb_file,hdrec, nsamps, nlines, &
       & nbands, storder, iswap, data_type, date, time, lat,lat_hem,&
       & long,long_hem,zenith,azimuth,altitude,sensor_type,x_start,y_start,&
       & default_bands, image_data_quantity)
    !Subroutine written 1999-09-28 by MJM
    !Modified 1999-09-28 by MJM.
    !1999-10-20 MJM: added optional argument sensor_type
    !1999-10-21 MJM: discovered I had not put stop if byte order ne 0 or 1,
    !                so I put a stop there. 
    !1999-10-25 MJM: Broke up sensor type if check into two. -MJM
    !2000-04-13 MJM: Added x_start and y_start optional keywords
    !2002-08-01 MJM: Added default bands optional keyword - can read one or
    !                three elements, but always returns three.
    !2005-11-29 MJM: Added ability to read ENVI mapping parameters. Those 
    !                are not passed out, but are stored for use only in this module
    !                in the output header files. 
    !2008-02-26 MJM: Used chartools eqic
    use errors_mod,only:fatal_error
    use file_io_mod,only:byteorder
    use chartools
    type(key), intent(in), dimension(:)::found_keys !found keys from a run of search_for_keys
    character (len=*), dimension(:), intent(in)::comb_file !The file, as an array
    integer, intent(out),target::hdrec !header size
    integer, intent(out),target::nsamps !Samples
    integer, intent(out),target::nlines   !Lines
    integer, intent(out),target::nbands !bands
    integer, intent(out)::storder   !Storage order (bip,bil,bsq)
    integer, intent(out),target::data_type !data type
    integer, intent(out)::iswap     !1 if need to swap, 0 if don't
    integer, intent(out),dimension(3),optional,target::date,default_bands
    real, intent(out), dimension(3), optional,target::time,lat,long
    real, intent(out), dimension(3), optional,target::zenith,azimuth
    character (len=1),intent(out),optional::lat_hem,long_hem
    real,intent(out),optional,target::altitude
    character (len=*), optional, target, intent(out)::sensor_type
    integer, intent(out), optional, target::x_start, y_start
    character (len=*), optional, target, intent(out)::image_data_quantity

    integer,target::byte_order
    character (len=charlen),target::interleave,lnghem,lathem
    integer::num_found_keys,i,j,tmp_int,bhold
    character (len=19),parameter::subroutine_name='newread_envi_header'

    integer,dimension(1)::tmp
    integer::tmp_i
    
    bhold=byteorder()
    num_found_keys=size(found_keys(:)%int_value)
    wanted_int_keys(:)%name='' !paranoia
    wanted_rel_keys(:)%name='' !paranoia
    wanted_str_keys(:)%name=''  
    wanted_ar3_keys(:)%name=''
    wanted_iar_keys(:)%name=''
    
    wanted_int_keys(1)%name='samples' ; wanted_int_keys(1)%value=>nsamps; nsamps=-1
    wanted_int_keys(2)%name='lines' ; wanted_int_keys(2)%value=>nlines; nlines=-1
    wanted_int_keys(3)%name='bands' ; wanted_int_keys(3)%value=>nbands; nbands=-1
    wanted_int_keys(4)%name='header offset' ; wanted_int_keys(4)%value=>hdrec; hdrec=-1
    wanted_int_keys(5)%name='data type' ; wanted_int_keys(5)%value=>data_type; data_type=-1
    wanted_int_keys(6)%name='byte order' ; wanted_int_keys(6)%value=>byte_order; byte_order=-1
    if (present(x_start)) then
       wanted_int_keys(7)%name='x start' ; wanted_int_keys(7)%value=>x_start; x_start=1
    end if
    if (present(y_start)) then
       wanted_int_keys(8)%name='y start' ; wanted_int_keys(8)%value=>y_start; y_start=1
    end if
    wanted_int_keys(:)%name_len=len_trim(wanted_int_keys(:)%name)

    
    if (present(altitude)) then
       wanted_rel_keys(1)%name='sensor_altitude'; wanted_rel_keys(1)%value=>altitude; altitude=-1.
    end if
    wanted_rel_keys(:)%name_len=len_trim(wanted_rel_keys(:)%name)

    wanted_str_keys(1)%name='interleave' ;wanted_str_keys(1)%value=>interleave ;interleave=''
    if (present(long_hem)) then
       wanted_str_keys(2)%name='image_center_long_hem' ;wanted_str_keys(2)%value=>lnghem; lnghem='' 
    end if
    if (present(lat_hem)) then
       wanted_str_keys(3)%name='image_center_lat_hem' ;wanted_str_keys(3)%value=>lathem; lathem=''
    end if
    if (present(sensor_type)) then
       wanted_str_keys(4)%name='sensor type'; wanted_str_keys(4)%value=>sensor_type; sensor_type=''
    end if
    if (present(image_data_quantity)) then
       wanted_str_keys(5)%name='image_data_quantity'
       wanted_str_keys(5)%value=>image_data_quantity
       image_data_quantity=''
    endif
    wanted_str_keys(:)%name_len=len_trim(wanted_str_keys(:)%name)

    if (present(time)) then
       wanted_ar3_keys(1)%name='image_center_time' ; wanted_ar3_keys(1)%value=>time; time=-1.
    end if
    if (present(long)) then
       wanted_ar3_keys(2)%name='image_center_long' ;wanted_ar3_keys(2)%value=>long; long=-1.
    end if
    if (present(lat)) then
       wanted_ar3_keys(3)%name='image_center_lat'  ;wanted_ar3_keys(3)%value=>lat; lat=-1.
    end if
    if (present(zenith)) then
       wanted_ar3_keys(4)%name='image_center_zenith_ang' ;wanted_ar3_keys(4)%value=>zenith; zenith=-1.
    end if
    if (present(azimuth)) then
       wanted_ar3_keys(5)%name='image_center_azimuth_ang';wanted_ar3_keys(5)%value=>azimuth; azimuth=-1.
    end if
    wanted_ar3_keys(:)%name_len=len_trim(wanted_ar3_keys(:)%name)

    if (present(date)) then
       wanted_iar_keys(1)%name='image_center_date' ; wanted_iar_keys(1)%value=>date; date=-1
    end if

    if (present(default_bands)) then
       wanted_iar_keys(2)%name='default bands' ; wanted_iar_keys(2)%value=>default_bands; default_bands=-1
    end if
    wanted_iar_keys(:)%name_len=len_trim(wanted_iar_keys(:)%name)
    
    do i=1,n_wanted_int_keys
       if (wanted_int_keys(i)%name.eq.'') cycle
       jloop:  do j=1,num_found_keys
          if (wanted_int_keys(i)%name(:wanted_int_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             wanted_int_keys(i)%value=found_keys(j)%int_value
             exit jloop
          end if
       end do jloop
    end do

    if (byte_order.lt.0.or.byte_order.gt.1) then
       print*,'byte order is not 0 or 1, so STOP'
       print*,'byte order = ',byte_order
       if (byte_order.eq.-1) print*,'byte order not found'
       call fatal_error(file_name,module_name,subroutine_name,&
            &'stopping: no valid order found')
    else
       if (byte_order.eq.bhold) then
          iswap=0
       else
          iswap=1
       end if
    end if

    do i=1,n_wanted_rel_keys
       if (wanted_rel_keys(i)%name.eq.'') cycle
       j2loop:  do j=1,num_found_keys
          if (wanted_rel_keys(i)%name(:wanted_rel_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             wanted_rel_keys(i)%value=found_keys(j)%rel_value
             exit j2loop
          end if
       end do j2loop
    end do

    do i=1,n_wanted_str_keys
       if (wanted_str_keys(i)%name.eq.'') cycle
       j3loop: do j=1,num_found_keys
          if (wanted_str_keys(i)%name(:wanted_str_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             wanted_str_keys(i)%value=''
             wanted_str_keys(i)%vlen=found_keys(j)%str_len
             wanted_str_keys(i)%value=found_keys(j)%str_value(:found_keys(j)%str_len)
             exit j3loop
          end if
       end do j3loop
    end do

    !The following two ifs need to split as they are since fortran does not 
    !guarantee the evaluation order.
    if (present(sensor_type)) then
       if (sensor_type.eq.'') sensor_type='unknown'
    end if

    call check_range('interleave',interleave(:3),(/'bip','bil','bsq'/))
    if (interleave(:3).eq.'bsq') then
       storder=bsq
    else if (interleave(:3).eq.'bip') then
       storder=bip
    else !bil since check_range above implies only choices are bip,bil,bsq
       storder=bil
    end if

    if (present(image_data_quantity)) then ! default data is radiance
       if (image_data_quantity.eq.'') image_data_quantity='at-sensor radiance'
    endif

    if (present(lat_hem)) lat_hem=lathem(:1)
    if (present(long_hem)) long_hem=lnghem(:1)

    do i=1,n_wanted_ar3_keys
       if (wanted_ar3_keys(i)%name.eq.'') cycle
       j5loop: do j=1,num_found_keys
          if (wanted_ar3_keys(i)%name(:wanted_ar3_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             call read_envi_value(found_keys(j)%first_line, &
                  &comb_file(found_keys(j)%first_line),&
                  & wanted_ar3_keys(i)%value,char_array=comb_file,lastrec=tmp_int)            
             exit j5loop
          end if
       end do j5loop
    end do

    ! Need more complicated structure since default bands can have one or three 
    ! elements.
    do i=1,n_wanted_iar_keys
       if (wanted_iar_keys(i)%name.eq.'') cycle
       j6loop: do j=1,num_found_keys
          if (wanted_iar_keys(i)%name(:wanted_iar_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             if (found_keys(j)%arr_size==1) then
                call read_envi_value(found_keys(j)%first_line, &
                     & comb_file(found_keys(j)%first_line),&
                     & tmp,char_array=comb_file,lastrec=tmp_int) !read 1 element array
                tmp_i=tmp(1) ! turn into scalar
                wanted_iar_keys(i)%value=tmp ! 'cuz a scalar can be assigned to array
             else if (found_keys(j)%arr_size==3) then
                call read_envi_value(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     & wanted_iar_keys(i)%value,char_array=comb_file,lastrec=tmp_int)
             else
                print*,'In the input or image header file the key ',found_keys(j)%name
                print*,'has ',found_keys(j)%arr_size,' elements.'
                print*,'The required number of elements is either 1 or 3 .'
                call fatal_error(file_name,module_name,subroutine_name,&
                     &'wrong number of elements')
             endif
             exit j6loop
          end if
       end do j6loop
    end do

    ! Now read geometric parameters. 
    ! pixel size = { } 2 or 3 elements, read as character array
    ! geo points = { } 4*n elements, x_i, y_i, lat_i, lon_i, real array
    ! map info = {} several elements, can vary, read as character array
    ! projection info = {} several elements, can vary, read as character array 
    wanted_map_keys(:)%name=''
    wanted_map_keys(1)%name='pixel size'
    wanted_map_keys(2)%name='map info'
    wanted_map_keys(3)%name='geo points'
    wanted_map_keys(4)%name='projection info'

    wanted_map_keys(:)%name_len=len_trim(wanted_map_keys(:)%name)
    
    do i=1,n_wanted_map_keys
       if (wanted_map_keys(i)%name.eq.'') cycle
       j7loop: do j=1,num_found_keys
          if (wanted_map_keys(i)%name(:wanted_map_keys(i)%name_len).eq.&
               found_keys(j)%name(:found_keys(j)%name_len)) then
             
             if (found_keys(j)%arr_size.le.1) cycle

             select case ((wanted_map_keys(i)%name(:wanted_map_keys(i)%name_len)))
                case('pixel size')
                   if (allocated(pixel_size)) deallocate(pixel_size)
                   allocate(pixel_size(found_keys(j)%arr_size))
                   wanted_map_keys(i)%value=>pixel_size
                case('map info')
                   if (allocated(map_info)) deallocate(map_info)
                   allocate(map_info(found_keys(j)%arr_size))
                   wanted_map_keys(i)%value=>map_info
                case('geo points')
                   if (allocated(geo_points)) deallocate(geo_points)
                   allocate(geo_points(found_keys(j)%arr_size))
                   wanted_map_keys(i)%value=>geo_points
                case('projection info')
                   if (allocated(projection_info)) deallocate(projection_info)
                   allocate(projection_info(found_keys(j)%arr_size))
                   wanted_map_keys(i)%value=>projection_info
             end select
             !print*,found_keys(j)
             call read_envi_string_array(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     & wanted_map_keys(i)%value,char_array=comb_file,lastrec=tmp_int)
             !print*,wanted_map_keys(i)%value
             exit j7loop
          endif

       end do j7loop
    end do

  end subroutine newread_envi_header
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine read_envi_header_nbarrays(found_keys, comb_file,nbands,&
       &wavobs,imscale,fwhm,units,band_names,bbl)
    !Subroutine written 1999-09-28 by MJM
    !Last modified 1999-09-28 by MJM.
    !2002-08-07 Added support for bbl. MJM
    !2008-02-26 MJM Use chartools eqic
    use errors_mod,only:fatal_error
    use chartools
    type(key), intent(in), dimension(:)::found_keys !found keys from a run of search_for_keys
    character (len=*), dimension(:), intent(in)::comb_file !The file, as an array
    integer, intent(in)::nbands !needed for some tests
    real,dimension(:), intent(out),optional,target::imscale,wavobs,fwhm
    character (len=*), dimension(:),intent(out),optional,target::units,band_names
    integer,dimension(:),intent(out),optional,target::bbl

    integer::i,j,tmp_int,num_found_keys
    real::tmp_real
    real, dimension(1)::tmp
    character (len=charlen),dimension(1)::tmp_ch_ar=''
    character (len=charlen)::tmp_ch=''
    character (len=25),parameter::subroutine_name='read_envi_header_nbarrays'
    integer, dimension(1)::tmpi
    integer::tmpint
    

    num_found_keys=size(found_keys(:)%int_value)

    wanted_anb_keys(:)%name=''
    wanted_char_keys(:)%name=''
    wanted_inb_keys(:)%name=''

    if (present(wavobs)) then
       wanted_anb_keys(1)%name='wavelength'; wanted_anb_keys(1)%value=>wavobs
    end if
    if (present(fwhm)) then
       wanted_anb_keys(2)%name='fwhm';       wanted_anb_keys(2)%value=>fwhm
    end if
    if (present(imscale)) then
       wanted_anb_keys(3)%name='image_scale_factor'; wanted_anb_keys(3)%value=>imscale
    end if
    wanted_anb_keys(:)%name_len=len_trim(wanted_anb_keys(:)%name)
    
    if (present(units)) then
       wanted_char_keys(1)%name='image_unscaled_units'; wanted_char_keys(1)%value=>units
    end if
    if (present(band_names)) then
       wanted_char_keys(2)%name='band names'; wanted_char_keys(2)%value=>band_names
    end if
    wanted_char_keys(:)%name_len=len_trim(wanted_char_keys(:)%name)

    if (present(bbl)) then
       wanted_inb_keys(1)%name='bbl'; wanted_inb_keys(1)%value=>bbl
    endif
    wanted_inb_keys(:)%name_len=len_trim(wanted_inb_keys(:)%name)
    
    !also do units and bands names here
    do i=1,n_wanted_anb_keys
       if (wanted_anb_keys(i)%name.eq.'') cycle 
       j4loop: do j=1,num_found_keys
          if (wanted_anb_keys(i)%name(:wanted_anb_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             if (found_keys(j)%arr_size==1) then !array of one elemnt found
                call read_envi_value(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     &tmp,char_array=comb_file,lastrec=tmp_int)
                tmp_real=tmp(1) !make single element array a scalar
                wanted_anb_keys(i)%value=tmp_real  !fill array with scalar
             else if (found_keys(j)%arr_size==nbands) then
                call read_envi_value(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     & wanted_anb_keys(i)%value,char_array=comb_file,lastrec=tmp_int)            
             else
                print*,'In the input or image header file the key ',found_keys(j)%name
                print*,'has ',found_keys(j)%arr_size,' elements.'
                print*,'The required number of elements is either 1 or ',nbands, '.'
                call fatal_error(file_name,module_name,subroutine_name,&
                     &'wrong number of elements')
             end if
             exit j4loop
          end if
       end do j4loop
    end do

    do i=1,n_wanted_char_keys
       if (wanted_char_keys(i)%name.eq.'') cycle
       j7loop: do j=1,num_found_keys
          if (wanted_char_keys(i)%name(:wanted_char_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             if (found_keys(j)%arr_size.eq.1) then
                call read_envi_value(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     & tmp_ch_ar,char_array=comb_file,lastrec=tmp_int)
                tmp_ch=tmp_ch_ar(1)
                wanted_char_keys(i)%value=tmp_ch
             else
                call read_envi_value(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     & wanted_char_keys(i)%value,char_array=comb_file,lastrec=tmp_int)
             end if
             exit j7loop
          end if
       end do j7loop
    end do
  
    do i=1,n_wanted_inb_keys
       if (wanted_inb_keys(i)%name.eq.'') cycle 
       j8loop: do j=1,num_found_keys
          if (wanted_inb_keys(i)%name(:wanted_inb_keys(i)%name_len).eqic.&
               &found_keys(j)%name(:found_keys(j)%name_len)) then
             if (found_keys(j)%arr_size==1) then !array of one elemnt found
                call read_envi_value(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     &tmpi,char_array=comb_file,lastrec=tmp_int)
                tmpint=tmpi(1) !make single element array a scalar
                wanted_inb_keys(i)%value=tmpint  !fill array with scalar
             else if (found_keys(j)%arr_size==nbands) then
                call read_envi_value(found_keys(j)%first_line, &
                     &comb_file(found_keys(j)%first_line),&
                     & wanted_inb_keys(i)%value,char_array=comb_file,lastrec=tmp_int)            
             else
                print*,'In the input or image header file the key ',found_keys(j)%name
                print*,'has ',found_keys(j)%arr_size,' elements.'
                print*,'The required number of elements is either 1 or ',nbands, '.'
                call fatal_error(file_name,module_name,subroutine_name,&
                     &'wrong number of elements')
             end if
             exit j8loop
          end if
       end do j8loop
    end do

  end subroutine read_envi_header_nbarrays
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine read_envi_header(iout, hdr_size, num_samples, num_lines, &
       & num_bands, storder, iswap, iret)
    use file_io_mod,only:byteorder
    integer, intent(in)::iout ! logical unit number
    !below, in order: return code, header size, number of samples, number of lines, 
    !number of bands, storage order
    integer, intent(out)::iret, hdr_size, num_samples, num_lines
    integer,intent(out)::num_bands, storder, iswap
    character (len=charlen)::junk,junk2 
    integer::ilen,border,ios,ikey,ieq
    integer::data_type, interleave_type, byte_order,bhold
    !Once again, MS Fortran barfs when it shouldn't. In particular, for the
    !below internal reads, if integer value is less than 10, and there's a 
    !preceding space (after the "=" and before the integer; i.e. "lines = 1"),
    !then for '*' formatting, MS fortran _always_ assigns the integer value
    !a value of 0. If the integer is >=10, everything works perfectly. 
    !Because of this, I needed to put in a few hacks in order to ensuure that
    !I am only looking at a number, and not at spaces, when reading from an 
    !internal file. 
    !1999-04-05 MJM
    !MS F90 seems to barf unless I protect the variable by assigning 
    !byteorder() to a variable:
    bhold=byteorder()

    !
    iret=0
    data_type=-999
    interleave_type=-999
    byte_order=-999
    num_bands=-999
    storder=-999
    hdr_size=-999
    num_samples=-999
    num_lines=-999
    iswap=0 !don't swap data when reading input file ; =1 means do swap
    read(iout,'(A4)')junk(:4)
    if (junk(:4).ne.'ENVI') then
       write(*,*)'The .hdr file is not ENVI compatible.'
       write(*,*)''
       iret=-1
       return
    endif
    do  ! loop until end-of-file condition leads to exit
       read(iout,'(A)',iostat=ios)junk
!       print*,junk
       if (ios.lt.0) exit ! end of file, no more processing possible
       if (verify(junk,"samples=0123456789 ") == 0) then
          ikey=index(junk,'samples')
          ieq=index(junk,'=')
          if (ikey*ieq.gt.0.and.ieq.gt.ikey) then
             ilen=len_trim(junk)
             junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen))
             ilen=len_trim(junk2)
             read(junk2(1:ilen),*)num_samples
             if (num_samples.le.0) then
                write(*,*)'My reading of the header file yields'
                write(*,*)'samples = ',num_samples,' .'
                write(*,*)'Either enter header info in input file,'
                write(*,*)'fix header file, or edit header reading module.'
                write(*,*)''
                iret=-1
             endif
          endif
       else if (verify(junk,"lines=0123456789 ") == 0) then
          ikey=index(junk,'lines')
          ieq=index(junk,'=')
          if (ikey*ieq.gt.0.and.ieq.gt.ikey) then
             ilen=len_trim(junk)
             junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen))
             ilen=len_trim(junk2)
             read(junk2(1:ilen),*)num_lines
             if (num_lines.le.0) then
                write(*,*)'My reading of the header file yields'
                write(*,*)'lines = ',num_lines,' .'
                write(*,*)'Either enter header info in input file,'
                write(*,*)'fix header file, or edit header reading module.'
                write(*,*)''
                iret=-1
             endif
          endif
       else if (verify(junk,"bands=0123456789 ") == 0) then
          ikey=index(junk,'bands')
          ieq=index(junk,'=')
          if (ikey*ieq.gt.0.and.ieq.gt.ikey) then
             ilen=len_trim(junk)
             junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen))
             ilen=len_trim(junk2)
             read(junk2(1:ilen),*)num_bands
             if (num_bands.le.0) then
                write(*,*)'My reading of the header file yields'
                write(*,*)'bands = ',num_bands,' .'
                write(*,*)'Either enter header info in input file,'
                write(*,*)'fix header file, or edit header reading module.'
                write(*,*)''
                iret=-1
             endif
          endif
       else if (verify(junk,"headerofst=0123456789 ") == 0) then
          ikey=index(junk,'header offset')
          ieq=index(junk,'=')
          if (ikey*ieq.gt.0.and.ieq.gt.ikey) then
             ilen=len_trim(junk)
             junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen))
             ilen=len_trim(junk2)
             read(junk2(1:ilen),*)hdr_size
             if (hdr_size.ne.0) then
                write(*,*)'The .hdr file reports that the data ', &
                     & 'cube has a non-zero header length, thus this',&
                     & 'file cannot be processed by this version of ',&
                     & 'Tafkaa.'
                write(*,*)''
                iret=-1
                return
             endif
          endif
       else if (verify(junk,"dataype=0123456789 ") == 0) then
          ikey=index(junk,'data type')
          ieq=index(junk,'=')
          if (ikey*ieq.gt.0.and.ieq.gt.ikey) then
             ilen=len_trim(junk)
             junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen))
             ilen=len_trim(junk2)
             read(junk2(1:ilen),*)data_type
             if (data_type.ne.2) then ! not returned, just tested
                write(*,*)'The header reports that the data',&
                     &' type is NOT I*2, so Tafkaa will not properly',&
                     &' read the input file.'
                write(*,*)'Reported data type=',data_type
                write(*,*)''
                iret=-1
                return
             endif
          endif
       else if (verify(junk,"interlavbsqp= ") == 0) then
          ikey=index(junk,'interleave')
          ieq=index(junk,'=')
          if (ikey*ieq.gt.0.and.ieq.gt.ikey) then
             ilen=len_trim(junk)
             if (index(junk(ieq+1:ilen),'bsq').gt.0) then
                storder=BSQ
             else if (index(junk(ieq+1:ilen),'bip').gt.0) then
                storder=BIP
             else if (index(junk(ieq+1:ilen),'bil').gt.0) then 
                storder=BIL
             else
                write(*,*)'The header header reports an unsupported',&
                     & ' value for the interleave parameter.'
                write(*,*)''
                iret=-1
                return
             endif
          endif
       else if (verify(junk,"byteord=01 ") == 0) then
          ikey=index(junk,'byte order')
          ieq=index(junk,'=')
          if (ikey*ieq.gt.0.and.ieq.gt.ikey) then
             ilen=len_trim(junk)
             junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen))
             ilen=len_trim(junk2)
             read(junk2(1:ilen),*)byte_order
             if (byte_order.eq.1) then
                border=1
             else if (byte_order.eq.0) then
                border=0
             else
                border=-1
             endif
             if (border.eq.-1) then
                write(*,*)'I cannot properly read the byte order from ',&
                     &'the header file. Please check the header file',&
                     &' and verify the byte order is either 0 or 1.'
                write(*,*)''
                iret=-1
             else if (border.ne.bhold) then
                write(*,*)'The header reports the byteorder is',&
                     &' not the native byte order for this machine. '
                write(*,*)'Tafkaa will attempt to swap bytes on the fly. '
                write(*,*)'There will be a small timing penalty because of this,'
                write(*,*)'but you won''t need two copies of the data lying around.'
                write(*,*)''
                iswap=1
                return
             endif
          endif
       endif
    enddo  ! whew. Long do loop; will exit only at end of .hdr file

    if (num_bands.le.0) then
       write(*,*)'num_bands = ',num_bands,' is out of range.'
       write(*,*)''
       iret=-1
    endif
    
    if (num_samples.le.0) then
       write(*,*)'num_samples = ',num_samples,' is out of range.'
       write(*,*)''
       iret=-1
    endif
    
    if (num_lines.le.0) then
       write(*,*)'num_lines = ',num_lines,' is out of range.'
       write(*,*)''
       iret=-1
    endif
    
    if (storder.eq.-999) then
       write(*,*)'I could not determine the storage order from the header file.'
       write(*,*)''
       iret=-1
    endif

    if (hdr_size.eq.-999) then
       write(*,*)'I could not determine the header size from the header file.'
       write(*,*)'Therefore: use the results with caution; verify the header file'
       write(*,*)'is correct; verify the header size is, indeed, 0.'
       write(*,*)''
    endif

    if (data_type.eq.-999) then
       write(*,*)''
       write(*,*)'I could not verify the data type to be I*2 from the header file.'
       write(*,*)'Therefore: use the results with caution; verify the header file'
       write(*,*)'is correct; verify the data is, indeed, I*2.'
    endif

    if (byte_order.eq.-999) then
       write(*,*)''
       write(*,*)'I could not verify the byte order to be ',bhold,' from the header file.'
       write(*,*)'Therefore: use the results with caution; verify the header file.'
       write(*,*)'is correct; verify the byte order is indeed ',bhold,'.'
    endif


    
  end subroutine read_envi_header
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module envi_header_mod
