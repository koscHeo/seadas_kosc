module check_range_mod

  implicit none
  private

  public::check_range
  interface check_range
     module procedure check_range_real, check_range_integer, &
          & check_range_string, check_range_real_array, &
          & check_range_integer_array, check_range_string_array
  end interface

contains
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine check_range_real(var_name, var, min, max)
    !Modified 2008-08-25 MJM: replaced s(:len_trim(s)) with trim(s)
    !Modified 2000-06-29 MJM: Added call to err_stop
    !Modified 2000-06-29 MJM: Added "FATAL ERROR" part of message
    !Modified 1999-09-23 by MJM
    implicit none
    character (len=*), intent(in)::var_name
    real, intent(in)::var, min, max
    !Range should be min<=var<max

    if (var.lt.min) then
       print*,'FATAL ERROR:'
       print*,trim(var_name),' < ',min
       print*,'and is out of range. Please correct input and continue.'
       print*,'STOP'
       call err_stop
    else if (var.ge.max) then
       print*,'FATAL ERROR:'
       print*,trim(var_name),' >= ',max
       print*,'and is out of range. Please correct input and continue.'
       print*,'STOP'
       call err_stop
   end if
  end subroutine check_range_real
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine check_range_real_array(var_name, var,mini,maxi)
    !Modified 2000-06-29 MJM: Added call to err_stop
    !Modified 2000-06-29 MJM: Added "FATAL ERROR" part of message
    !Last modified 1999-09-23 by MJM
    implicit none
    character (len=*), intent(in), dimension(:)::var_name
    real, intent(in), dimension(:)::var,mini,maxi
    integer::n,i

    if (size(mini).ne.size(var)) then
       print*,'FATAL ERROR:'
       print*,'check_range_real_array: dimensions of min and var disagree'
       print*,'                        for variable',var_name
       call err_stop
    else if (size(maxi).ne.size(var)) then
       print*,'FATAL ERROR:'
       print*,'check_range_real_array: dimensions of max and var disagree'
       print*,'                        for variable',var_name
       call err_stop
    end if

    n=size(var)
    do i=1,n
       call check_range_real(var_name(i),var(i),mini(i),maxi(i))
    end do

  end subroutine check_range_real_array
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine check_range_integer(var_name, var, min, max)
    !Modified 2008-08-25 MJM: replaced s(:len_trim(s)) with trim(s)
    !Modified 2000-06-29 MJM: Added call to err_stop
    !Modified 2000-06-29 MJM: Added "FATAL ERROR" part of message
    !Last modified 1999-09-23 by MJM
    implicit none
    character (len=*), intent(in)::var_name
    integer, intent(in)::var, min, max
    !Integers are easier. The allowed range is: min<=var<=max

    if (var.lt.min) then
       print*,'FATAL ERROR:'
       print*,trim(var_name),' < ',min
       print*,'and is out of range. Please correct input and continue.'
       print*,'STOP'
       call err_stop
    else if (var.gt.max) then
       print*,'FATAL ERROR:'
       print*,trim(var_name),' > ',max
       print*,'and is out of range. Please correct input and continue.'
       print*,'STOP'
       call err_stop
    end if
  end subroutine check_range_integer
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine check_range_integer_array(var_name, var,mini,maxi)
    !Modified 2000-06-29 MJM: Added call to err_stop
    !Modified 2000-06-29 MJM: Added "FATAL ERROR" part of message
    !Last modified 1999-09-23 by MJM
    implicit none
    character (len=*), intent(in), dimension(:)::var_name
    integer, intent(in), dimension(:)::var,mini,maxi
    integer::n,i

    if (size(mini).ne.size(var)) then
       print*,'FATAL ERROR:'
       print*,'check_range_integer_array: dimensions of min and var disagree'
       print*,'                        for variable',var_name
       call err_stop
    else if (size(maxi).ne.size(var)) then
       print*,'FATAL ERROR:'
       print*,'check_range_integer_array: dimensions of max and var disagree'
       print*,'                        for variable',var_name
       call err_stop
    end if

    n=size(var)
    do i=1,n
       call check_range_integer(var_name(i),var(i),mini(i),maxi(i))
    end do

  end subroutine check_range_integer_array
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine check_range_string(var_name, var, set,trim_front,match_index)
    !2008-Aug-27 MJM: optionally return index of value that matches
    !                 Added option to trim leading balnks before comparison 
    !                 default is trim_front=.true.
    !Modified 2000-06-29 MJM: Added call to err_stop
    !Modified 2000-06-29 MJM: Added "FATAL ERROR" part of message
    !1999-10-26 MJM: checked for vlen=0 in case we're sent a null value.
    !1999-10-26 MJM: fixed bug: vhold was declared len=len(var_name),
    !                which is wrong
    ! Note that default (as originally written) is to trim leading and trailing spaces
    ! (blanks) for both 'set' and 'var'
    ! Note that the default Fortran comparison ignores trailing blanks, i.e., the 
    ! if a and b are character variables of different lengths, then a==b is true 
    ! as long as a and b are identical except for trailing blanks (the shorter one
    ! is padded with blanks and then the two strings (now of identical length) are compared.
    !Created 1999-09-23 by MJM
   implicit none
    character (len=*), intent(in)::var_name
    character (len=*), intent(in)::var
    character (len=*), dimension(:), intent(in)::set
    logical, intent(in),optional::trim_front
    !logical, intent(in),optional::trim_trail
    integer, intent(out),optional:: match_index
    integer::N,vlen,i
    integer, dimension(size(set))::ilen
    character (len=len(var))::vhold
    character (len=len(set)), dimension(size(set))::shold
    character (len=len(set))::junk
    logical::trim_f,trim_t
    ! For the string case, we want to verify that var is a member of set. 
    ! We'll trim leading and trailing spaces on both var and set 
    ! (actually, on their local copies: neither var nor set is changed). 

    if (present(trim_front)) then
       trim_f=trim_front
    else
       trim_f=.true. ! default is to ignore leading blanks
    endif
    !if (present(trim_trail)) then
    !   trim_t=trim_trail
    !else
    !   trim_t=.true. ! default is to ignore trailing blanks
    !endif
    

    vlen=0
    N=size(set) ! need to get the correct number of elements
    vhold=''
    vlen=len_trim(var)
    !if (vlen.eq.0) then
    !   print*,'FATAL ERROR:'
    !   print*,trim(var_name),' = ',''
    !   print*,' and was apparently not set.'
    !   print*,'stopping'
    !   call err_stop
    !end if
    if (trim_f) then 
       vhold(:vlen)=adjustl(var(:vlen)) ! trim leading blanks
    else
       vhold(:vlen)=var(:vlen)
    endif
    vlen=len_trim(vhold) !now vhold(:vlen) has no leading or trailing spaces
    
    if (present(match_index)) match_index=0
    
    do i=1,N
       ilen(i)=len_trim(set(i))
       shold(i)=''
       
       if (trim_f) then 
          shold(i)(:ilen(i))=adjustl(set(i)(:ilen(i))) ! trim leading blanks
       else
          shold(i)(:ilen(i))=set(i)(:ilen(i))
       endif
       !The Microsoft Power Fortran 4.0 compiler does not like the next line-
       !at least with the compiler flags I use- so it will be replaced.
       ! ilen=len_trim(shold(i)) !Now shold(i)(:ilen) has no leading or trailing sp
       junk=''; junk(1:ilen(i))=shold(i)(1:ilen(i)); ilen(i)=len_trim(junk)

       if (vhold(:vlen).eq.shold(i)(:ilen(i))) then
          if (present(match_index)) match_index=i
          return !successful match
       end if
    end do
    !If you get here, you didn't match anything in set.That means you should 
    !print an error message and stop
    print*,'FATAL ERROR:'
    print*,trim(var_name),' = ',vhold(:vlen)
    print*,'was not an allowed member of the set:'
    write(*,'(A)')(shold(i)(1:ilen(i)),i=1,N)
    print*,'Please correct input, and try again.'
    print*,'STOPPING'
    call err_stop
  end subroutine check_range_string
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine check_range_string_array(var_name, var, set,trim_front,match_index)
    !2008 Aug 27 MJM: optionally return array (size(var)) of match indices
    !                 Added control for trimming front (default is true) which ignores
    !                 leading blanks
    !Last modified 1999-09-23 by MJM
    implicit none
    character (len=*), intent(in)::var_name
    character (len=*), dimension(:), intent(in)::var
    character (len=*), dimension(:), intent(in)::set
    logical,optional,intent(in)::trim_front
    integer, intent(out),dimension(:), optional::match_index

    integer:: N, i,vlen
    integer,dimension(size(var))::match_hold
    logical::trim_f
    
    if (present(trim_front)) then
       trim_f=trim_front
    else
       trim_f=.true.
    endif
 
 
    !In this version, var_name is a scalar
    n=size(var)
    vlen=len_trim(var_name)
    
    if (present(match_index)) then
       if (size(match_index).ne.n) then ! Incorrectly written caller
          print*,'In subroutine check_range_string_array'
          print*,'var_name = ',trim(var_name)
          print*,'dimensions of var = ',n
          print*,'dimensions of match_index = ',size(match_index)
          print*,'dimensions do not match - a programming issue'
          print*,'Please report error'
          print*,'STOPPING'
          call err_stop
       endif 
    endif
    match_hold=0 ! always has the write dimensions, always initialized
    
    ! Always call check_range_string here with match hold. 
    do i=1,n
       call check_range_string(trim(var_name), var(i), set, &
              trim_front=trim_f,match_index=match_hold(i))
    end do
    ! Only set match_index (and pass values) if present
    if (present(match_index)) match_index=match_hold

  end subroutine check_range_string_array
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine err_stop
    ! Created 2000-06-29 by MJM
    implicit none
    character (len=1)::ans

    write(*,*)''
    write(*,*)'PLEASE COPY THE ABOVE MESSAGES BEFORE TYPING ANY KEYS AT ALL:'
    write(*,*)'Depending on the implementation and how this program was called,'
    write(*,*)'the window may disappear following your response to the question'
    write(*,*)'below, so this _may_ be the last chance to copy the messages.'
    write(*,*)''
    write(*,*)'To stop, ENTER any character'
    read(*,*)ans
    stop
  end subroutine err_stop
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
end module check_range_mod
