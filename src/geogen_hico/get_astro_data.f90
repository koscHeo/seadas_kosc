module get_ast_data
  implicit none
  real (kind=8), parameter:: Pi      =  4.0D0 * DATAN(1.0D0)
  real (kind=8), parameter:: HalfPi  =  0.5D0*pi
  real (kind=8), parameter:: TwoPi   =  2.0D0*pi
  real (kind=8), parameter:: Rad2Deg = 180.0D0/pi
  real (kind=8), parameter:: Deg2Rad = pi/180.0D0
  real (kind=8), parameter:: as2deg = 1.d0/3600.d0
  real (kind=8), parameter:: as2rad = as2deg*deg2rad

  character (len=3),dimension(12),parameter::MMM=['JAN','FEB','MAR',&
       'APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

  private
  public::get_astro_data

contains
  
  subroutine get_astro_data(leapsec_file,earth_orient_file,year,mon,day,xp,yp,lod,dut1,dat,err)
    use params,only:clen
    implicit none
!    character (len=*),intent(in)::ast_data_dir
    character (len=*),intent(in)::leapsec_file
    character (len=*),intent(in)::earth_orient_file
    integer, intent(in)::year,mon,day
    real (kind=8),intent(out)::xp,yp,lod,dut1,dat
    integer, intent(out)::err
    !
    character (len=clen)::filename
    integer::lun,yf,mf,df,cc,yy,i,yfo,dfo,mfo,idate,mdate,done,ios
    real (kind=8)::mjd,pmx,pmx_e,pmy,pmy_e,ut1mutc,ut1mutc_e,loda,loda_e,taimutc,taimutco
    character (len=1)::pmflag,utflag
    character (len=3)::mc,junk,mco
    character (len=8)::junk2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! only partially read lines, and stop when find needed date
    cc=year/100
    yy=modulo(year,100)

    open(newunit=lun,file=trim(earth_orient_file),access='sequential',action='read',status='old',form='formatted')
    err=1 ! assume not found
    ios=0
    do 
       read(lun,fmt='(3i2,1x,f8.2,1x,a1,1x,2f9.6,1x,2f9.6,2x,a1,2f10.7,1x,2f7.4)', &
            iostat=ios)& 
            yf,mf,df,mjd,pmflag,pmx,pmx_e,pmy,pmy_e,utflag,ut1mutc,ut1mutc_e,loda,loda_e
    
       if (ios.ne.0) then
          print*,'Returning: Error in read in get_astro_data'
          print*,'Last good line: ',yf,mf,df,mjd
          err=ios
          return
       endif

       if ((yf.eq.yy).and.(mf.eq.mon).and.(df.eq.day)) then 
          xp=pmx*as2rad ! convert arcsec to rad
          yp=pmy*as2rad ! ditto
          lod=loda/1000.d0 ! convert ms to s
          dut1=ut1mutc
          err=0
          exit
       endif
    enddo
    close(lun)
    if ((yf.ne.yy).or.(mf.ne.mon).or.(df.ne.day)) then
       err=10
       print*,'Date not found in file '
       print*,'File searched: ',trim(filename)
       print*,'year, mon, day: ',year, mon, day
       return
    endif

    mdate=year*10000 + mon*100 + day
    open(newunit=lun,file=trim(leapsec_file),access='sequential',action='read',status='old',form='formatted')
! skip the first line as it is a header
    read(lun,*)
    read(lun,fmt='(1x,i4,1x,a3,1x,i2,1x,a3,1x,f9.1,2x,a8,2x,f10.7)')yfo,mco,dfo,junk,mjd,junk2,taimutco
    do i=1,12
       if (mco.eq.mmm(i)) mfo=i
    end do
    done=0
    do
       read(lun,fmt='(1x,i4,1x,a3,1x,i2,1x,a3,1x,f9.1,2x,a8,2x,f10.7)',end=20)yf,mc,df,junk,mjd,junk2,taimutc
       do i=1,12
          if (mc.eq.mmm(i)) mf=i
       end do
       idate=10000*yf+100*mf+df
       if (mdate.gt.idate) then
          taimutco=taimutc
          done=0
          cycle
       else
          dat=taimutco
          done=1
          exit
       endif
    end do
20  continue
    if (done.ne.1) dat=taimutco

    close(lun)
    
  end subroutine get_astro_data


end module get_ast_data

