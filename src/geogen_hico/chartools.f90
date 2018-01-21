module chartools

implicit none

! Hmm, can't get universal private and specific publics to work and still
! get the .eqic. function to be recognized. So everything is public. 
!private
!public:: to_upper, to_lower, eqic, .eqic.

interface operator (.eqic.)
   module procedure eqic
end interface

contains

  subroutine crush(string)
    !2012/11/14-16 MJM
    ! Crushes two spaces to one, removes trailing spaces.
    ! output string should have no more than one space separating "words"
    ! based on routine I wrote a long time ago.
    ! Truly a brute force solution
    implicit none
    character (len=*),intent(inout):: string
    ! locals
    character(len=512):: copy
    integer::i,ilength,ilen
    
    ilen=len_trim(string) ! no spaces at end
    copy=''
    copy(1:ilen)=string(1:ilen)
    ilength=ilen
    do i=1,ilen-1
       if (i.eq.ilength) exit
       if (string(i:i).eq.achar(9)) string(i:i)=char(32) ! convert tab to space
       do while (verify(string(i:i+1),achar(9)//achar(32)).eq.0)
          string(i:ilen-1)=copy(i+1:ilength)
          ilength=ilength-1
          copy='' ! paranoia?
          copy(1:ilength)=string(1:ilength)
       end do
    enddo
  end subroutine crush

  pure function to_upper(in) result(out)
    ! 2004Aug31 MJM Assumes standard ascii collating sequence. 
    ! In the ASCII sequence, [A:Z] = 65:90, [a:z]=[97:122]
    ! Does not change other characters
    implicit none
    character (len=*),intent(in):: in
    character (len=len(in)):: out
    character (len=1)::c
    integer::i,n
    
    N=len(in)
    if (N.eq.0) then !strinctly paranoia
       out=in
       return
    end if
    
    do i=1,N
       c=in(i:i)
       select case(c)
       case ('a':'z') ! lowercase letters
          out(i:i)=char(iachar(c)-32)
       case default   ! everything else
          out(i:i)=c
       end select
    end do
    
  end function to_upper
  
  pure function to_lower(in) result(out)
    ! 2004Aug31 MJM Assumes standard ascii collating sequence. 
    ! In the ASCII sequence, [A:Z] = 65:90, [a:z]=[97:122]
    ! Does not change other characters.
    implicit none
    character (len=*),intent(in):: in
    character (len=len(in)):: out
    character (len=1)::c
    integer::i,n
    
    N=len(in)
    if (N.eq.0) then ! strictly paranoia
       out=in
       return
    end if
    
    do i=1,N
       c=in(i:i)
       select case(c)
       case ('A':'Z') ! uppercase letters
          out(i:i)=char(iachar(c)+32)
       case default   ! everything else
          out(i:i)=c
       end select
    end do
    
  end function to_lower

  pure function eqic(a,b) result(output)
    ! 2004Aug31: character equivalence ignoring case
    ! output = .true. if match when case is ignored. Otherwise
    ! output will be false.
    implicit none
    character (len=*),intent(in)::a,b
    logical::output
    integer::na,nb
    
    na=len(a)
    nb=len(b)
    select case (nb-na)
       case (0)     ! same length, so try
         ! Might as well switch them both to lower case, and then test.
          output= to_lower(a) == to_lower(b)
       case default ! not the same length, so won't match
          output=.false.
    end select
  end function eqic
  
!  function has_letters(a) result(output)
!     implicit none
!     character(len=*),intent(in)::a
! 
!     
!  end function has_letters
  
end module chartools

!program testchar
!  use chartools
!  character (len=80)::junk=''
!!  
!  junk='ABCDEFGHIJKLMNOPQRSTUVWXYZ_1234,./-abcdefghijklmnopqrstuvwxyz'
!  
!  write(*,*)to_lower(junk)
!  write(*,*)to_upper(junk)
!
!  write(*,*)to_lower(junk) .eqic. to_upper(junk)
!
!  
!end program testchar
