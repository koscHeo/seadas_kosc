!
!     !!! Note: To use include "aquarius.fin" in calling routine !!!
!
      subroutine open_hdf5(filename, dsetname, file_id, dset_id)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      character(len=*), intent (in) :: filename
      character(len=*), intent (in) :: dsetname

      INTEGER(HID_T) :: file_id       ! File identifier 
      INTEGER(HID_T) :: dset_id       ! Dataset identifier 

      INTEGER     ::   error ! Error flag

      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error) 

      !
      ! Open an existing file.
      !
      CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(filename)//' not found.'
         stop
      endif

      !
      ! Open an existing dataset. 
      !
      CALL h5dopen_f(file_id, dsetname, dset_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(dsetname)//' not found.'
         stop
      endif

      END


      subroutine create_hdf5(filename, file_id)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      character(len=*), intent (in) :: filename

      INTEGER(HID_T) :: file_id       ! File identifier 

      INTEGER     ::   error ! Error flag

      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error) 

      !
      ! Create a new file using default properties.
      ! 
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      END



      subroutine create_hdf5_real(file_id, dsetname, rank, dims, dset_id)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T), intent (in) :: file_id       ! File identifier 
      character(len=*), intent (in) :: dsetname

      INTEGER :: rank
      INTEGER(HSIZE_T), DIMENSION(:) :: dims

      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

      INTEGER     ::   error ! Error flag

      !
      ! Initialize FORTRAN predefined datatypes.
      !
      CALL h5open_f(error)

      ! 
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)

      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, dset_id, error)

      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(dspace_id, error)

      END



      subroutine create_hdf5_i32le(file_id, dsetname, rank, dims, dset_id)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T), intent (in) :: file_id       ! File identifier 
      character(len=*), intent (in) :: dsetname

      INTEGER :: rank
      INTEGER(HSIZE_T), DIMENSION(:) :: dims

      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier

      INTEGER     ::   error ! Error flag

      !
      ! Initialize FORTRAN predefined datatypes.
      !
      CALL h5open_f(error)

      ! 
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)

      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_STD_I32LE, dspace_id, dset_id, error)

      !
      ! Terminate access to the data space.
      !
      CALL h5sclose_f(dspace_id, error)

      END


      subroutine close_hdf5(file_id, dset_id, dataspace, filespace)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T) :: file_id       ! File identifier 
      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace

      INTEGER     ::   error ! Error flag

      CALL h5sclose_f(dataspace, error) 
      CALL h5sclose_f(filespace, error) 

      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)

      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
     
      !
      ! Close FORTRAN interface.
      !
!      CALL h5close_f(error)

      END


      subroutine close_hdf5_df(file_id)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T) :: file_id       ! File identifier 

      INTEGER     ::   error ! Error flag

      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
     
      !
      ! Close FORTRAN interface.
      !
!      CALL h5close_f(error)

      END


      subroutine close_hdf5_ds(dset_id, dataspace, filespace)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace

      INTEGER     ::   error ! Error flag

      CALL h5sclose_f(dataspace, error) 
      CALL h5sclose_f(filespace, error) 

      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)

      END


      subroutine set_space_hdf5(dset_id, rank, count, filespace, dataspace)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER :: rank
      INTEGER(HSIZE_T), DIMENSION(:) :: count
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace

      INTEGER     ::   error ! Error flag

!      write (*,*) count

      ! Select hyperslab
      call h5dget_space_f(dset_id, filespace, error) 

      ! Get dataset's dataspace identifier.
      CALL h5dget_space_f(dset_id, dataspace, error)

      ! Create memory dataspace.
      CALL h5screate_simple_f(rank, count, dataspace, error)

      END


      subroutine read_hdf5_real(dset_id, filespace, dataspace, offset, count, buffer)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace
      INTEGER(HSIZE_T), DIMENSION(:) :: offset
      INTEGER(HSIZE_T), DIMENSION(:) :: count
      REAL*4, DIMENSION(:) :: buffer

      INTEGER     ::   error ! Error flag

!      write (*,*) offset
!      write (*,*) count

      ! Select hyperslab in the dataset.
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

      ! Read the dataset.
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, buffer, count, error, dataspace, filespace, H5P_DEFAULT_F)

      END



      subroutine write_hdf5_real(dset_id, filespace, dataspace, offset, count, buffer)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace

      INTEGER(HSIZE_T), DIMENSION(:) :: offset
      INTEGER(HSIZE_T), DIMENSION(:) :: count
      REAL*4, POINTER :: buffer

      INTEGER     ::   error ! Error flag

!      write (*,*) offset
!      write (*,*) count

      ! Select hyperslab in the dataset.
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

      ! Write the dataset.
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, buffer, count, error, dataspace, filespace, H5P_DEFAULT_F)

      END



      subroutine write_hdf5_i32le(dset_id, filespace, dataspace, offset, count, buffer)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HID_T) :: filespace
      INTEGER(HID_T) :: dataspace

      INTEGER(HSIZE_T), DIMENSION(:) :: offset
      INTEGER(HSIZE_T), DIMENSION(:) :: count
      INTEGER*4, POINTER :: buffer

      INTEGER     ::   error ! Error flag

!      write (*,*) offset
!      write (*,*) count

      ! Select hyperslab in the dataset.
      CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

      ! Write the dataset.
      CALL h5dwrite_f(dset_id, H5T_STD_I32LE, buffer, count, error, dataspace, filespace, H5P_DEFAULT_F)

      END


      subroutine dump_hdf5(filename, dsetname, data_dims, buffer)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      character(len=*), intent (in) :: filename
      character(len=*), intent (in) :: dsetname

      INTEGER(HID_T) :: file_id       ! File identifier 
      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HSIZE_T), DIMENSION(:) :: data_dims

      INTEGER     ::   error ! Error flag

      REAL*4, POINTER :: buffer

!      write (*,*) filename
!      write (*,*) dsetname

      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error) 

      !
      ! Open an existing file.
      !
      CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(filename)//' not found.'
         call exit(1)
      endif

      !
      ! Open an existing dataset. 
      !
      CALL h5dopen_f(file_id, dsetname, dset_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(dsetname)//' not found.'
         call exit(1)
      endif

!      write(*,*) data_dims(1)
      !
      ! Read the dataset.
      !
      CALL h5dread_f(dset_id, H5T_NATIVE_REAL, buffer, data_dims, error)

      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)

      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
     
      !
      ! Close FORTRAN interface.
      !
!      CALL h5close_f(error)

      END


      subroutine dump_hdf5_int(filename, dsetname, data_dims, buffer)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      character(len=*), intent (in) :: filename
      character(len=*), intent (in) :: dsetname

      INTEGER(HID_T) :: file_id       ! File identifier 
      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HSIZE_T), DIMENSION(:) :: data_dims

      INTEGER     ::   error ! Error flag

      INTEGER*4, POINTER :: buffer

!      write (*,*) filename
!      write (*,*) dsetname

      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error) 

      !
      ! Open an existing file.
      !
      CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(filename)//' not found.'
         call exit(1)
      endif

      !
      ! Open an existing dataset. 
      !
      CALL h5dopen_f(file_id, dsetname, dset_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(dsetname)//' not found.'
         call exit(1)
      endif

      !
      ! Read the dataset.
      !
      CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, buffer, data_dims, error)

      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)

      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
     
      !
      ! Close FORTRAN interface.
      !
!      CALL h5close_f(error)

      END



      subroutine dump_hdf5_dbl(filename, dsetname, data_dims, buffer)

      USE HDF5 ! This module contains all necessary modules 
        
      IMPLICIT NONE

      character(len=*), intent (in) :: filename
      character(len=*), intent (in) :: dsetname

      INTEGER(HID_T) :: file_id       ! File identifier 
      INTEGER(HID_T) :: dset_id       ! Dataset identifier 
      INTEGER(HSIZE_T), DIMENSION(:) :: data_dims

      INTEGER     ::   error ! Error flag

      REAL*8, POINTER :: buffer

!      write (*,*) filename
!      write (*,*) dsetname

      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error) 

      !
      ! Open an existing file.
      !
      CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(filename)//' not found.'
         call exit(1)
      endif

      !
      ! Open an existing dataset. 
      !
      CALL h5dopen_f(file_id, dsetname, dset_id, error)
      if (error.ne.0) then
         write(*,*) adjustr(dsetname)//' not found.'
         call exit(1)
      endif


      !
      ! Read the dataset.
      !
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, data_dims, error)

      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)

      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
     
      !
      ! Close FORTRAN interface.
      !
!      CALL h5close_f(error)

      END

