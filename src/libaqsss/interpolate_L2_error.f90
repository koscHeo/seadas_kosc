!   simplified routine for calculating randon and systematic uncertainties
!   based on climatology
!   T. Meissner, RSS, April 20, 2015

 
    subroutine interpolate_L2_error (err_map_file, idayjl, xlat, xlon, iasc,  err_ran, err_sys)
      USE HDF5                  ! This module contains all necessary modules 
      implicit none

      include "aquarius.fin"

      character(len=150), intent(in) :: err_map_file

      integer(4), intent(in)      ::  idayjl  ! julian day of Lw observation 1,...,366
      real(4), intent(in)         ::  xlat    ! between -90 and +90
      real(4), intent(in)         ::  xlon    ! between -180 and +180
      integer(4), intent(in)      ::  iasc    ! =1 ascending  (0  <=zang <=180)
                                                ! =2 descending (180< zang < 360)
        
      real(4), intent(out)        :: err_ran  ! random     component of L2 salinity error [psu]
      real(4), intent(out)        :: err_sys  ! systematic component of L2 salinity error [psu]
        
      real(4)                                         ::  ylon, ylat, ymon
      integer(4)                                      ::  ilon, ilat, imon        
      integer(4)                                      ::  jlon, jlat, jmon
      real(4)                                         ::  ylon0,ylat0, ymon0      

      real(4)                                         ::  b1, b2, c1, c2, wt(4), bmon
      integer(4)                                      ::  n
      real(8), dimension(0:2)                         ::  tsum
      integer(4)                                      ::  j, k

      real(4)                                         ::  err_ran_a, err_sys_a, err_ran_b, err_sys_b
        

      INTEGER(HID_T) :: file_id       ! File identifier 
      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: dataspace     ! Dataspace identifier 
      INTEGER(HID_T) :: memspace      ! memspace identifier 

      INTEGER(HSIZE_T), DIMENSION(1:4) :: memdims = (/14,362,180,2/)

      INTEGER(HSIZE_T), DIMENSION(1:4) :: offset_in_mem = (/1,1,0,0/)
      INTEGER(HSIZE_T), DIMENSION(1:4) :: count_in_mem = (/12,360,180,2/)

      INTEGER     ::   error ! Error flag


      integer(4), dimension(0:13,0:361,180,2), save   ::  num_map !, num_map2
      real(4),    dimension(0:13,0:361,180,2), save   ::  err_map_sys, err_map_ran
      integer(4), save                                ::  istart=1
      integer(4), parameter                           ::  iunit=3
        
      real(4), parameter                              ::  err_missing = -9999.0 !1.0E30  ! very large value if no entry
        
      ! default = missing
      err_ran = err_missing
      err_sys = err_missing
        
        
      if (istart==1) then
            
         istart=0
            
         ! Initialize FORTRAN interface.
         CALL h5open_f(error) 

         ! Open an existing file.
         CALL h5fopen_f (err_map_file, H5F_ACC_RDONLY_F, file_id, error)
         if (error.ne.0) then
            write(*,*) adjustr(err_map_file)//' not found.'
            call exit(1)
         endif

         ! Define hyperslab in the memory.
         CALL h5screate_simple_f(4, memdims, memspace, error)
         CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
              offset_in_mem, count_in_mem, error)

         CALL h5dopen_f(file_id, 'l2_num_map', dset_id, error)
         CALL h5dget_space_f(dset_id, dataspace, error)
         CALL h5sselect_all_f(dataspace, error)
         CALL H5dread_f(dset_id, H5T_NATIVE_INTEGER, num_map, memdims, error, &
              memspace, dataspace)
         CALL h5sclose_f(dataspace, error)
         CALL h5dclose_f(dset_id, error)

         CALL h5dopen_f(file_id, 'l2_err_map_ran', dset_id, error)
         CALL h5dget_space_f(dset_id, dataspace, error)
         CALL h5sselect_all_f(dataspace, error)
         CALL H5dread_f(dset_id, H5T_NATIVE_REAL, err_map_ran, memdims, error, &
              memspace, dataspace)
         CALL h5sclose_f(dataspace, error)
         CALL h5dclose_f(dset_id, error)

         CALL h5dopen_f(file_id, 'l2_err_map_sys', dset_id, error)
         CALL h5dget_space_f(dset_id, dataspace, error)
         CALL h5sselect_all_f(dataspace, error)
         CALL H5dread_f(dset_id, H5T_NATIVE_REAL, err_map_sys, memdims, error, &
              memspace, dataspace)
         CALL h5sclose_f(dataspace, error)
         CALL h5dclose_f(dset_id, error)

         CALL h5sclose_f(memspace, error)
         CALL h5fclose_f(file_id, error)

!         open(unit=iunit,file='/disk02/joel/ocssw/wentz/V4.0/L2_uncertainty_maps_AD_V4.0.dat', &
!              access='stream',form='unformatted',status='old',action='read')
!         read(iunit) num_map2     (1:12,1:360,:,:)
!         read(iunit) err_map_ran(1:12,1:360,:,:)
!         read(iunit) err_map_sys(1:12,1:360,:,:)
!         close(iunit)    
            
         ! time wrap
         num_map    ( 0,1:360,:,:)=num_map    (12,1:360,:,:)
         num_map    (13,1:360,:,:)=num_map    ( 1,1:360,:,:)
         err_map_ran( 0,1:360,:,:)=err_map_ran(12,1:360,:,:)
         err_map_ran(13,1:360,:,:)=err_map_ran( 1,1:360,:,:)
         err_map_sys( 0,1:360,:,:)=err_map_sys(12,1:360,:,:)
         err_map_sys(13,1:360,:,:)=err_map_sys( 1,1:360,:,:)
            
!         num_map2    ( 0,1:360,:,:)=num_map2    (12,1:360,:,:)
!         num_map2    (13,1:360,:,:)=num_map2    ( 1,1:360,:,:)

         ! longitude wrap
         num_map    (:,  0,:,:)=num_map    (:,360,:,:)
         num_map    (:,361,:,:)=num_map    (:,  1,:,:)
         err_map_ran(:,  0,:,:)=err_map_ran(:,360,:,:)
         err_map_ran(:,361,:,:)=err_map_ran(:,  1,:,:)
         err_map_sys(:,  0,:,:)=err_map_sys(:,360,:,:)
         err_map_sys(:,361,:,:)=err_map_sys(:,  1,:,:)

!         num_map2    (:,  0,:,:)=num_map    (:,360,:,:)
!         num_map2   (:,361,:,:)=num_map    (:,  1,:,:)               
      endif
        
      ylat = xlat
      ylon = xlon
      if (ylon > 180.)     ylon =ylon-360.
      if (ylon <-179.9999) ylon =-179.9999
      if (ylon > 179.9999) ylon = 179.9999
      if (ylat < -89.4999) ylat = -89.4999
      if (ylat >  89.4999) ylat =  89.4999
                
        
      ilat = floor(ylat - (-89.5)) + 1
      if (ilat <  1)  ilat = 1
      if (ilat >180)  ilat = 180
        
      ilon = floor(ylon - (-179.5)) + 1
      if (ilon < 0 .or. ilon > 360)  then
         write(*,*) xlon,ylon,ilon,' lon oob in interpolate_L2_error'
         call exit(1)
      endif
        
      if(idayjl < 1 .or. idayjl > 366)  then
         write(*,*) idayjl,' idayl  oob in interpolate_L2_error'
         call exit(1)
      endif
 
      ymon=12.*(idayjl - 0.5)/365.25
      if(ymon.gt.11.9999) ymon=11.9999
      imon = nint(ymon) 
        
      if (iasc<1 .or. iasc>2) then
         write(*,*) iasc, ' iasc oob in interpolate_L2_error'
         call exit(1)
      endif

      !if (num_map(imon,ilon,ilat,iasc) < 1) return ! no data     

!      if (num_map(imon,ilon,ilat,iasc) .ne. num_map2(imon,ilon,ilat,iasc)) then
!         stop
!      endif

      jlon=ilon+1
      jlat=ilat+1
      jmon=imon+1
        
      ylon0 = (ilon-1)*1.0 + (-179.5)
      b2  = (ylon - ylon0)/1.0
      if (b2<0.0) b2=0.0
      if (b2>1.0) b2=1.0
      b1  = 1.0-b2
      if (b1<0.0) b1=0.0
      if (b1>1.0) b1=1.0      
        
        
      ylat0 = (ilat-1)*1.0 + (-89.5)
      c2    = (ylat - ylat0)/1.0
      if (c2<0.0) c2=0.0
      if (c2>1.0) c2=1.0      
      c1    = 1.0-c2
      if (c1<0.0) c1=0.0
      if (c1>1.0) c1=1.0
        
      ymon0 = (imon-1)*1.0 + 0.5
      bmon  = (ymon - ymon0)/1.0      
      if (bmon<0) bmon=0.0
      if (bmon>1) bmon=1.0
        
      wt(1)=b1*c1
      wt(2)=b1*c2
      wt(3)=b2*c1
      wt(4)=c2*b2
    
      err_ran_a = err_missing
      err_sys_a = err_missing
      tsum=0.d0
      n=0
      do j = ilon,jlon
         do k = ilat,jlat
            n=n+1
            if (num_map(imon,j,k,iasc) < 1) cycle       
            tsum(0) = tsum(0) + wt(n)
            tsum(1) = tsum(1) + wt(n)*err_map_ran(imon,j,k,iasc)**2
            tsum(2) = tsum(2) + wt(n)*err_map_sys(imon,j,k,iasc)
         enddo
      enddo
      if (tsum(0)>0.0001) then
         err_ran_a = tsum(1)/tsum(0)
         err_sys_a = tsum(2)/tsum(0)
      endif
        

      err_ran_b = err_missing
      err_sys_b = err_missing
      tsum=0.d0
      n=0
      do j = ilon,jlon
         do k = ilat,jlat
            n=n+1
            if (num_map(jmon,j,k,iasc) < 1) cycle       
            tsum(0) = tsum(0) + wt(n)
            tsum(1) = tsum(1) + wt(n)*err_map_ran(jmon,j,k,iasc)**2
            tsum(2) = tsum(2) + wt(n)*err_map_sys(jmon,j,k,iasc)
         enddo
      enddo
      if (tsum(0)>0.0001) then
         err_ran_b = tsum(1)/tsum(0)
         err_sys_b = tsum(2)/tsum(0)
      endif
        

      if (abs(err_ran_a-err_missing)<0.1 .and. abs(err_ran_b-err_missing)<0.1 ) then
         err_ran = err_missing
      endif
    
      if (abs(err_ran_a-err_missing)<0.1 .and. abs(err_ran_b-err_missing)>=0.1 ) then
         err_ran = sqrt(err_ran_b)
      endif
        
      if (abs(err_ran_a-err_missing)>=0.1 .and. abs(err_ran_b-err_missing)<0.1 ) then
         err_ran = sqrt(err_ran_a)
      endif
        
      if (abs(err_ran_a-err_missing)>=0.1 .and. abs(err_ran_b-err_missing)>=0.1 ) then
         err_ran = err_ran_a*(1.0-bmon) + err_ran_b*bmon
         err_ran = sqrt(err_ran)
      endif


      if (abs(err_sys_a-err_missing)<0.1 .and. abs(err_sys_b-err_missing)<0.1 ) then
         err_sys = err_missing
      endif
    
      if (abs(err_sys_a-err_missing)<0.1 .and. abs(err_sys_b-err_missing)>=0.1 ) then
         err_sys = err_sys_b
      endif
        
      if (abs(err_sys_a-err_missing)>=0.1 .and. abs(err_sys_b-err_missing)<0.1 ) then
         err_sys = err_sys_a
      endif
        
      if (abs(err_sys_a-err_missing)>=0.1 .and. abs(err_sys_b-err_missing)>=0.1 ) then
         err_sys = err_sys_a*(1.0-bmon) + err_sys_b*bmon
      endif


      return 
    end subroutine interpolate_L2_error

