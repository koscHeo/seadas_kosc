module file_io_mod

  ! 2001-04-30 MJM Added functions for swapping many different ranks,
  !                added generic features, added more basic guts_of_swap
  !                now used by most swapping subroutines here.
  ! 2000-12-07 MJM Added ability to swap 4 byte integers and reals.
  !                Modified swap algorithm slightly to use vector
  !                subscripting and get rid of do variables.
  ! Notes: The functions mbyteorder, swap_scalar_i2, 
  ! and swap_vector_i2 are not strict (that is ANSI standard) F90 because 
  ! in three of the functions below I equivalence integers of 
  ! different kind. 
  ! 
  ! I think that the only way to deal with this in ANSI F90 is
  ! to read the values byte-by-byte from an external file, 
  ! write them in two byte swapped order to a different external file, then 
  ! re-read them as 2 byte integers from the second external file. 
  ! Now, that would probably work quite well for the byteorder() function, 
  ! but may be a bit unwieldy for the swap_i2 parts of the module.
  ! Actually, it works, and without too much of a timing penalty. 
  ! It has been implemented in f_byteorder, f_swap_scalar_i2, and 
  ! f_swap_vector_i2. 
  ! -Marcos Montes 1998 November 4
  implicit none
  private
  public:: byteorder, get_file_unit_number, swap_i2, swap_i4,swap_r4, swap_endian
  
!!$interface swap_i2
!!$   module procedure swap_scalar_i2, swap_vector_i2
!!$end interface
  
  interface byteorder
     module procedure f_byteorder
  end interface
  
  interface swap_endian
     module procedure f_swap_scalar_i2, f_swap_vector_i2, f_swap_rank2_i2, &
          & f_swap_rank3_i2,f_swap_rank4_i2,f_swap_rank5_i2,f_swap_rank6_i2, &
          & f_swap_rank7_i2, f_swap_scalar_i4, f_swap_vector_i4, f_swap_rank2_i4, &
          & f_swap_rank3_i4,f_swap_rank4_i4,f_swap_rank5_i4,f_swap_rank6_i4, &
          & f_swap_rank7_i4, f_swap_scalar_r4, f_swap_vector_r4, f_swap_rank2_r4, &
          & f_swap_rank3_r4,f_swap_rank4_r4,f_swap_rank5_r4,f_swap_rank6_r4, &
          & f_swap_rank7_r4
  end interface

  interface swap_i2
     module procedure f_swap_scalar_i2, f_swap_vector_i2, f_swap_rank2_i2, &
          & f_swap_rank3_i2,f_swap_rank4_i2,f_swap_rank5_i2,f_swap_rank6_i2, &
          & f_swap_rank7_i2
  end interface

  interface swap_i4
     module procedure f_swap_scalar_i4, f_swap_vector_i4, f_swap_rank2_i4, &
          & f_swap_rank3_i4,f_swap_rank4_i4,f_swap_rank5_i4,f_swap_rank6_i4, &
          & f_swap_rank7_i4
  end interface

  interface swap_r4
     module procedure f_swap_scalar_r4, f_swap_vector_r4, f_swap_rank2_r4, &
          & f_swap_rank3_r4,f_swap_rank4_r4,f_swap_rank5_r4,f_swap_rank6_r4, &
          & f_swap_rank7_r4
  end interface
  
  integer, parameter::kind4=selected_int_kind(9) !999,999,999 < (2^(4*8))/2
  integer, PARAMETER::kind2=selected_int_kind(4) !9,999 < (2^(2*8))/2
  integer, parameter::kind1=selected_int_kind(2) !99 <  (2^(1*8))/2
  
  ! Note this won't work with the "wrong" kind of numerical models.
  ! I am not sure how to enforce Real*4 requirements with selected 
  ! real kind. Use default real, assume 4 byte?????
  
  character (len=11),parameter::file_name='file_io.f90'
  character (len=11),parameter::module_name='file_io_mod'
  
contains
  
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine get_file_unit_number(outnum)
    !     returns a valid, unused logical unit number in the range of 10-100
    !     that should be valid on many systems. 
    ! Marcos Montes, 1998 November 4
    use errors_mod,only:fatal_error
    implicit none
    integer,intent(out)::outnum
    integer::i
    logical::iex,iop
    character (len=20),parameter::subroutine_name='get_file_unit_number'

    do i=10,100 !we'll only use units 10 to 100, inclusive
       inquire(unit=i,exist=iex,opened=iop)
       if (iex.and..not.iop) then ! unit exists (ie is legal?) and not open
          outnum=i
          exit ! leave loop
       endif
       if (i.eq.100) then  !if i=100 and you're here, must have failed
          outnum=-1 ! error code
          call fatal_error(file_name,module_name,subroutine_name,&
               &'No valid logical unit numbers available')
          ! stop here since I don't always check for non-negativeness.
       endif
    enddo
  end subroutine get_file_unit_number
!!$  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  !     function to determine native byteorder, and return what ENVI wants. Works
!!$  !     on SGI and PC. Requires no input values, and returns an integer,
!!$  !     either a 0 or a 1. -Marcos Montes 1998 Aug 14, NRL Code 7212
!!$  ! According to the ENVI manual: 0= Least Significant Byte first (DEC and MS-DOS)
!!$  !                               1= Most Significant Byte first (all other systems).
!!$  integer function mbyteorder()  !=1 for SGI, =0 for INTEL; same as ENVI
!!$    implicit none
!!$    integer (kind=kind4)::num=1
!!$    integer (kind=kind2),dimension(2)::a
!!$    equivalence (num,a) ! Valid and useful use of equivalence!
!!$    !NOT STRICT F90 BECAUSE WE'RE EQUIVALENCING VARIABLES OF DIFFERENT KIND TOGETHER.
!!$
!!$    mbyteorder=a(2)
!!$    if (kind2.ge.kind4) then
!!$       write(*,*)'**** VALUES REPORTED IN BYTEORDER WILL BE WRONG ****'
!!$    endif
!!$  end function mbyteorder
!!$  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  ! Marcos Montes, 1998 November 4
!!$  integer (kind=kind2) function swap_scalar_I2(i2_in) result(swapped_i2)
!!$    implicit none
!!$    integer (kind=kind2), intent(in)::i2_in
!!$    integer (kind=kind2),save::dummy_in, i2_swapped
!!$    integer (kind=kind1), dimension(2),save::dummy_out, swapped_out
!!$    equivalence (dummy_in, dummy_out)  !before swap
!!$    equivalence (i2_swapped, swapped_out) !after swap
!!$    ! byte swap an I2; does what dd (on Unix system) does, but slower.
!!$    ! Could (optionally) be called on various systems if disk space is a 
!!$    ! problem, I suppose. Currently not used in any of the processing. 
!!$    ! THIS IS NOT STRICT F90 BECAUSE IN F90 ONE IS NOT SUPPOSED TO 
!!$    ! EQUIVALENCE DIFFERENT KINDs TOGETHER.
!!$
!!$    dummy_in=i2_in ! put the input value into a dummy value equivalenced to an I1(:2)
!!$
!!$    swapped_out(1)=dummy_out(2) !swap first part
!!$    swapped_out(2)=dummy_out(1) !swap second part
!!$
!!$    swapped_i2=i2_swapped  !pass to the result
!!$
!!$  end function swap_scalar_i2
!!$  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  !Efforts to make a dynamically allocated array valued swap_i2:
!!$  !apparently one is not allowed to equivalence arrays without fixed size. Sigh. 
!!$  !SO one needs arrays with fixed sizes. Let's see what we can fix up. 
!!$  ! Marcos Montes 1998 November 4 
!!$  function swap_vector_i2(vector_in) result(vector_out)
!!$    implicit none
!!$    integer (kind=kind2), intent(in), dimension(:):: vector_in
!!$    !make n_biggest bigger than an AVIRIS line
!!$    integer, parameter:: n_biggest=650*250 
!!$    integer (kind=kind2), dimension(size(vector_in))::vector_out
!!$    integer (kind=kind2), dimension(n_biggest),save::before_swap, after_swap
!!$    integer (kind=kind1), dimension(2*n_biggest),save::i1_before_swap,i1_after_swap
!!$    integer::insize, i, insize_2
!!$    !Note: THIS IS NOT STRICT F90 because in F90 you are not supposed to EQUIVALENCE 
!!$    !variables of different KINDs.
!!$    equivalence (before_swap,i1_before_swap)
!!$    equivalence (after_swap,i1_after_swap)
!!$
!!$    insize=size(vector_in)
!!$    insize_2=2*insize
!!$    if (insize.le.n_biggest) then !smaller than n_biggest
!!$       before_swap(:insize)=vector_in(:insize)
!!$       
!!$       i1_after_swap(1:insize_2:2)=i1_before_swap(2:insize_2:2)
!!$       i1_after_swap(2:insize_2:2)=i1_before_swap(1:insize_2:2)
!!$
!!$       vector_out(:insize)=after_swap(:insize)
!!$    else ! fall back to scalar swap, explicitly called, if size>n_biggest
!!$       do i=1,insize
!!$          vector_out(i)=swap_scalar_i2(vector_in(i))
!!$       end do
!!$    endif
!!$
!!$  end function swap_vector_i2
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer function f_byteorder() 
    !     function to determine native byteorder, and return what ENVI wants. Works
    !     on SGI and PC. Requires no input values, and returns an integer,
    !     either a 0 or a 1. -Marcos Montes 1998 Aug 14, NRL Code 7212
    ! According to the ENVI manual: 
    !         0= Least Significant Byte first (DEC and MS-DOS)
    !         1= Most Significant Byte first (all other systems).
    ! ANSI standard method of determining byteorder. 
    ! Not as elegant as methods involving equivalences, but it works.
    ! -Marcos Montes 1998 November 4
    implicit none
    integer (kind=kind4)::num=1
    integer (kind=kind2),dimension(2)::a
    integer :: scratchfile, i4len

    call get_file_unit_number(scratchfile)
    inquire(iolength=i4len)num
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=i4len)

    write(unit=scratchfile,rec=1)num
    read(unit=scratchfile,rec=1)a

    close(unit=scratchfile)

    f_byteorder=a(2)
    
    if (kind2.ge.kind4) then
       write(*,*)'**** VALUES REPORTED IN f_BYTEORDER WILL BE WRONG ****'
    endif
  end function f_byteorder
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer (kind=kind2) function f_swap_scalar_I2(i2_in) result(swapped_i2)
    !2000-12-06 MJM Use vector triplet in one line for swap.
    implicit none
    integer (kind=kind2), intent(in)::i2_in
    integer (kind=kind1), dimension(2)::dummy_out,swapped_out
    integer ::scratchfile,i2len
!   
    ! byte swap an I2; does what dd (on Unix system) does, but slower.
    ! Could (optionally) be called on various systems if disk space is a 
    ! problem, I suppose. 
    ! ANSI standard.
    ! Really slow on arrays; use f_swap_vector_i2 for those
    ! Marcos Montes 1998 November 4
    call get_file_unit_number(scratchfile)
    inquire(iolength=i2len)i2_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=i2len)
    write(unit=scratchfile,rec=1)i2_in
    read(unit=scratchfile,rec=1)dummy_out

    swapped_out(1:2:1)=dummy_out(2:1:-1)

    write(unit=scratchfile,rec=1)swapped_out
    read(unit=scratchfile,rec=1)swapped_i2
    close(unit=scratchfile)


  end function f_swap_scalar_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_vector_i2(vector_in) result(vector_out)
    ! 2000-04-30 MJM: use guts_of_swap
    ! 2000-12-07 MJM : Use vector triplet and 2d byte array now
    ! This should be portable. It works on the systems that I use. 
    ! It is ANSI standard. It is not nearly as elegant as the solution using 
    ! equivalences since it involves IO, but it actually executes pretty quickly. 
    ! 600 calls to this function only costs an extra 10-12 seconds on my SGI Octane.
    ! Marcos Montes 1998 November 4
    implicit none
    integer (kind=kind2), intent(in), dimension(:):: vector_in
    integer (kind=kind2), dimension(size(vector_in))::vector_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)vector_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)vector_in !write the I2
    call guts_of_swap(2,size(vector_in),scratchfile)
    read(unit=scratchfile,rec=1)vector_out     !read the swapped I2

    close(unit=scratchfile)

  end function f_swap_vector_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank2_i2(array_in) result(array_out)
    ! 2001-04-30 Written. MJM
    implicit none
    integer (kind=kind2), intent(in), dimension(:,:):: array_in
    integer (kind=kind2), dimension(size(array_in,dim=1), &
         & size(array_in,dim=2))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I2
    call guts_of_swap(2,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I2
    close(unit=scratchfile)

  end function f_swap_rank2_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank3_i2(array_in) result(array_out)
    ! 2001-04-30 Written. MJM
    implicit none
    integer (kind=kind2), intent(in), dimension(:,:,:):: array_in
    integer (kind=kind2), dimension(size(array_in,dim=1), &
         & size(array_in,dim=2),size(array_in,dim=3))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I2
    call guts_of_swap(2,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I2
    close(unit=scratchfile)

  end function f_swap_rank3_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank4_i2(array_in) result(array_out)
    ! 2001-04-30 Written. MJM
    implicit none
    integer (kind=kind2), intent(in), dimension(:,:,:,:):: array_in
    integer (kind=kind2), dimension(size(array_in,dim=1), &
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I2
    call guts_of_swap(2,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I2
    close(unit=scratchfile)

  end function f_swap_rank4_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank5_i2(array_in) result(array_out)
    ! 2001-04-30 Written. MJM
    implicit none
    integer (kind=kind2), intent(in), dimension(:,:,:,:,:):: array_in
    integer (kind=kind2), dimension(size(array_in,dim=1), &
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4),size(array_in,dim=5))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I2
    call guts_of_swap(2,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I2
    close(unit=scratchfile)

  end function f_swap_rank5_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank6_i2(array_in) result(array_out)
    ! 2001-04-30 Written. MJM
    implicit none
    integer (kind=kind2), intent(in), dimension(:,:,:,:,:,:):: array_in
    integer (kind=kind2), dimension(size(array_in,dim=1), &
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4),size(array_in,dim=5), &
         & size(array_in,dim=6))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I2
    call guts_of_swap(2,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I2
    close(unit=scratchfile)

  end function f_swap_rank6_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank7_i2(array_in) result(array_out)
    ! 2001-04-30 Written. MJM
    implicit none
    integer (kind=kind2), intent(in), dimension(:,:,:,:,:,:,:):: array_in
    integer (kind=kind2), dimension(size(array_in,dim=1), &
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4),size(array_in,dim=5), &
         & size(array_in,dim=6),size(array_in,dim=7))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I2
    call guts_of_swap(2,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I2
    close(unit=scratchfile)

  end function f_swap_rank7_i2
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  integer (kind=kind4) function f_swap_scalar_I4(i4_in) result(swapped_i4)
    !2000-12-07 MJM Use vector triplet in one line for swap.
    implicit none
    integer (kind=kind4), intent(in)::i4_in
    integer (kind=kind1), dimension(4)::dummy_out,swapped_out
    integer ::scratchfile,i4len
!   
    ! byte swap an I4; does what dd (on Unix system) does, but slower.
    ! Could (optionally) be called on various systems if disk space is a 
    ! problem, I suppose. 
    ! ANSI standard.
    ! Really slow on arrays; use f_swap_vector_i4 for those
    call get_file_unit_number(scratchfile)
    inquire(iolength=i4len)i4_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=i4len)
    write(unit=scratchfile,rec=1)i4_in
    read(unit=scratchfile,rec=1)dummy_out

    swapped_out(1:4:1)=dummy_out(4:1:-1)

    write(unit=scratchfile,rec=1)swapped_out
    read(unit=scratchfile,rec=1)swapped_i4
    close(unit=scratchfile)

  end function f_swap_scalar_I4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_vector_i4(vector_in) result(vector_out)
    ! 2000-12-07 MJM : Use vector triplet and 2d byte array now
    ! This should be portable. It works on the systems that I use. 
    ! It is ANSI standard. It is not nearly as elegant as the solution using 
    ! equivalences since it involves IO, but it actually executes pretty quickly. 
    ! 600 calls to this function only costs an extra 10-12 seconds on my SGI Octane.
    ! Marcos Montes 1998 November 4
    implicit none
    integer (kind=kind4), intent(in), dimension(:):: vector_in
    integer (kind=kind4), dimension(size(vector_in))::vector_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)vector_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)vector_in !write the I4
    call guts_of_swap(4,size(vector_in),scratchfile)
    read(unit=scratchfile,rec=1)vector_out     !read the swapped I4
    close(unit=scratchfile)

  end function f_swap_vector_i4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank2_i4(array_in) result(array_out)
    ! 2001-04-30 MJM Written.
    implicit none
    integer (kind=kind4), intent(in), dimension(:,:):: array_in
    integer (kind=kind4), dimension(size(array_in,dim=1),&
         & size(array_in,dim=2))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I4
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I4
    close(unit=scratchfile)

  end function f_swap_rank2_i4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank3_i4(array_in) result(array_out)
    ! 2001-04-30 MJM Written.
    implicit none
    integer (kind=kind4), intent(in), dimension(:,:,:):: array_in
    integer (kind=kind4), dimension(size(array_in,dim=1),&
         & size(array_in,dim=2),size(array_in,dim=3))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I4
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I4
    close(unit=scratchfile)

  end function f_swap_rank3_i4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank4_i4(array_in) result(array_out)
    ! 2001-04-30 MJM Written.
    implicit none
    integer (kind=kind4), intent(in), dimension(:,:,:,:):: array_in
    integer (kind=kind4), dimension(size(array_in,dim=1),&
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I4
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I4
    close(unit=scratchfile)

  end function f_swap_rank4_i4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank5_i4(array_in) result(array_out)
    ! 2001-04-30 MJM Written.
    implicit none
    integer (kind=kind4), intent(in), dimension(:,:,:,:,:):: array_in
    integer (kind=kind4), dimension(size(array_in,dim=1),&
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4),size(array_in,dim=5))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I4
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I4
    close(unit=scratchfile)

  end function f_swap_rank5_i4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank6_i4(array_in) result(array_out)
    ! 2001-04-30 MJM Written.
    implicit none
    integer (kind=kind4), intent(in), dimension(:,:,:,:,:,:):: array_in
    integer (kind=kind4), dimension(size(array_in,dim=1),&
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4),size(array_in,dim=5), &
         & size(array_in,dim=6))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I4
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I4
    close(unit=scratchfile)

  end function f_swap_rank6_i4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank7_i4(array_in) result(array_out)
    ! 2001-04-30 MJM Written.
    implicit none
    integer (kind=kind4), intent(in), dimension(:,:,:,:,:,:,:):: array_in
    integer (kind=kind4), dimension(size(array_in,dim=1),&
         & size(array_in,dim=2),size(array_in,dim=3), &
         & size(array_in,dim=4),size(array_in,dim=5), &
         & size(array_in,dim=6),size(array_in,dim=7))::array_out
    integer ::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile, status='scratch', access='direct', action='readwrite',&
         & form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in !write the I4
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out     !read the swapped I4
    close(unit=scratchfile)

  end function f_swap_rank7_i4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_scalar_r4(value_in) result(value_out)
    ! 2000-12-07 MJM Written and added to this module
    implicit none
    real (kind=4), intent(in)::value_in
    real (kind=4)::value_out
    integer (kind=kind1), dimension(4)::i1_before_swap,i1_after_swap
    integer::inlen,scratchfile

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)value_in
    open(unit=scratchfile,status='scratch',access='direct',recl=inlen, &
         action='readwrite',form='unformatted')
    write(unit=scratchfile,rec=1)value_in
    read(unit=scratchfile,rec=1)i1_before_swap

    i1_after_swap(1:4:1)=i1_before_swap(4:1:-1)
    write(unit=scratchfile,rec=1)i1_after_swap
    read(unit=scratchfile,rec=1)value_out
    close(unit=scratchfile)
    
  end function f_swap_scalar_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_vector_r4(vector_in) result(vector_out)
    ! 2001-04-30 MJM extended to use guts_of_swap subroutine
    ! 2000-12-07 MJM Written and added to this module
    implicit none
    real(kind=4), intent(in),dimension(:)::vector_in
    real(kind=4),dimension(size(vector_in))::vector_out
    integer::inlen,scratchfile

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)vector_in
    open(unit=scratchfile,status='scratch',action='readwrite',&
         &access='direct',form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)vector_in
    call guts_of_swap(4,size(vector_in),scratchfile)
    read(unit=scratchfile,rec=1)vector_out
    close(unit=scratchfile)

  end function f_swap_vector_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank2_r4(array_in) result (array_out)
    implicit none
    !2001-04-30 MJM written and added to module
    real(kind=4), intent(in),dimension(:,:)::array_in
    real(kind=4), dimension(size(array_in,dim=1),size(array_in,dim=2))::array_out
    integer::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile,status='scratch',action='readwrite',&
         &access='direct',form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out
    close(unit=scratchfile)

  end function f_swap_rank2_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank3_r4(array_in) result (array_out)
    implicit none
    !2001-04-30 MJM written and added to module
    real(kind=4), intent(in),dimension(:,:,:)::array_in
    real(kind=4), dimension(size(array_in,dim=1),size(array_in,dim=2),&
         & size(array_in,dim=3))::array_out
    integer::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile,status='scratch',action='readwrite',&
         &access='direct',form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out
    close(unit=scratchfile)

  end function f_swap_rank3_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank4_r4(array_in) result (array_out)
    implicit none
    !2001-04-30 MJM written and added to module
    real(kind=4), intent(in),dimension(:,:,:,:)::array_in
    real(kind=4), dimension(size(array_in,dim=1),size(array_in,dim=2),&
         & size(array_in,dim=3),size(array_in,dim=4))::array_out
    integer::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile,status='scratch',action='readwrite',&
         &access='direct',form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out
    close(unit=scratchfile)

  end function f_swap_rank4_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank5_r4(array_in) result (array_out)
    implicit none
    !2001-04-30 MJM written and added to module
    real(kind=4), intent(in),dimension(:,:,:,:,:)::array_in
    real(kind=4), dimension(size(array_in,dim=1),size(array_in,dim=2),&
         & size(array_in,dim=3),size(array_in,dim=4),size(array_in,dim=5))::array_out
    integer::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile,status='scratch',action='readwrite',&
         &access='direct',form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out
    close(unit=scratchfile)

  end function f_swap_rank5_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank6_r4(array_in) result (array_out)
    implicit none
    !2001-04-30 MJM written and added to module
    real(kind=4), intent(in),dimension(:,:,:,:,:,:)::array_in
    real(kind=4), dimension(size(array_in,dim=1),size(array_in,dim=2),&
         & size(array_in,dim=3),size(array_in,dim=4),size(array_in,dim=5),&
         & size(array_in,dim=6))::array_out
    integer::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile,status='scratch',action='readwrite',&
         &access='direct',form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out
    close(unit=scratchfile)

  end function f_swap_rank6_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function f_swap_rank7_r4(array_in) result (array_out)
    implicit none
    !2001-04-30 MJM written and added to module
    real(kind=4), intent(in),dimension(:,:,:,:,:,:,:)::array_in
    real(kind=4), dimension(size(array_in,dim=1),size(array_in,dim=2),&
         & size(array_in,dim=3),size(array_in,dim=4),size(array_in,dim=5),&
         & size(array_in,dim=6),size(array_in,dim=7))::array_out
    integer::scratchfile,inlen

    call get_file_unit_number(scratchfile)
    inquire(iolength=inlen)array_in
    open(unit=scratchfile,status='scratch',action='readwrite',&
         &access='direct',form='unformatted',recl=inlen)
    write(unit=scratchfile,rec=1)array_in
    call guts_of_swap(4,size(array_in),scratchfile)
    read(unit=scratchfile,rec=1)array_out
    close(unit=scratchfile)

  end function f_swap_rank7_r4
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine guts_of_swap(bytes_per_kind,n_elements,scratchfile)
    !2001-04-30 MJM This is the guts of swapping in fortran 90, as far as
    !               I can tell.
    implicit none
    integer, intent(in)::bytes_per_kind,n_elements,scratchfile
    integer (kind=kind1),dimension(bytes_per_kind,n_elements)::i1_before_swap,&
         & i1_after_swap

    read(unit=scratchfile,rec=1)i1_before_swap
    i1_after_swap(1:bytes_per_kind:1,:)=i1_before_swap(bytes_per_kind:1:-1,:)
    write(unit=scratchfile,rec=1)i1_after_swap
        
  end subroutine guts_of_swap
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module file_io_mod

