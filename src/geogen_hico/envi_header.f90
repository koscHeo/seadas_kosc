!!!!program envi_test
!!!!
!!!!use envi_header
!!!!
!!!!implicit none
!!!!
!!!!integer::samples,lines,bands,header,data_type,byte_order,ios
!!!!character (len=132)::junk
!!!!character (len=3)::interleave
!!!!character (len=13)::file_type
!!!!real, dimension(224)::wavelengths
!!!!character (len=80), dimension(224)::band_names
!!!!integer::i,num_text_lines,i_last
!!!!integer, dimension(224)::bbl
!!!!character (len=132), allocatable, dimension(:)::text_file
!!!!character (len=80)::out_key
!!!!integer::out_type,arr_size,first_rec
!!!!logical::out_scalar
!!!!character (len=3), dimension(3), parameter::type=(/'int','rel','str'/)

!scratch file
!!!!open(10,status='scratch',access='sequential',form='formatted',&
!!!!     &position='asis',action='readwrite')

! First count lines
!!!!i=0
!!!!do 
!!!!   junk=''
!!!!   read(5,'(A)',iostat=ios)junk
!!!!   if (ios.lt.0) exit
!!!!   write(10,'(A)')junk
!!!!   i=i+1
!!!!enddo
!!!!num_text_lines=i

!Allocate space
!!!!allocate(text_file(num_text_lines))
!!!!text_file=''

!read file into memory
!!!!rewind(10)
!!!!read(10,'(A)')(text_file(i),i=1,num_text_lines)
!!!!close(10)

! Test subroutine determine_key
!!!!i=0
!!!!do
!!!!   i=i+1
!!!!   out_key=''
!!!!   first_rec=i
!!!!   call determine_key(text_file,i,out_key,out_type,out_scalar,arr_size)
!!!!   if (out_key.ne.'') then
!!!!      print*,'key=','^'//out_key(:len_trim(out_key))//'^'
!!!!      print*,'first_rec=',first_rec,'   last_rec=',i
!!!!      print*,'scalar=',out_scalar
!!!!      print*,'arr_size=',arr_size
!!!!      print*,'out_type=',type(out_type)
!!!!      print*;
!!!!   end if
!!!!   if (i.eq.num_text_lines) exit
!!!!end do

!Test everything used by read_envi_*
!!!!i=0
!!!!do
!   read(5,'(A)',iostat=ios)junk
!!!!   i=i+1
!!!!   junk=''
!!!!   junk=text_file(i)
!   if (ios.lt.0) exit
!!!!   if (match_key(junk,'samples')) then
!!!!      call read_envi_value(junk,samples)
!!!!   else if (match_key(junk,'lines')) then
!!!!      call read_envi_value(junk,lines)
!!!!   else if (match_key(junk,'bands')) then
!!!!      call read_envi_value(junk,bands)
!!!!   else if (match_key(junk,'header offset')) then
!!!!      call read_envi_value(junk,header)
!!!!   else if (match_key(junk,'data type')) then
!!!!      call read_envi_value(junk,data_type)
!!!!   else if (match_key(junk,'byte order')) then
!!!!      call read_envi_value(junk,byte_order)
!!!!   else if (match_key(junk,'interleave')) then
!!!!      call read_envi_value(junk,interleave)
!!!!   else if (match_key(junk,'file type')) then
!!!!      call read_envi_value(junk,file_type)
!!!!   else if (match_array_key(junk,'wavelength')) then
!!!!      call read_envi_value(i,junk,wavelengths,text_file,i_last)
!!!!      i=i_last
!!!!   else if (match_array_key(junk,'band names')) then
!!!!      call read_envi_value(i,junk,band_names,text_file,i_last)
!!!!      i=i_last
!!!!   else if (match_array_key(junk,'bbl')) then
!!!!      call read_envi_value(i,junk,bbl,text_file,i_last)
!!!!      i=i_last
!!!!   end if
!!!!   if (i == num_text_lines) exit
!!!!end do

!!!!print*,'samples=',samples
!!!!print*,'bands=',bands
!!!!print*,'lines=',lines
!!!!print*,'header=',header
!!!!print*,'data type=',data_type
!!!!print*,'byte order=',byte_order
!!!!print*,'interleave=',interleave(:len_trim(interleave))
!!!!print*,'file type=',file_type(:len_trim(file_type))
!!!!write(*,'(A,7(F9.6,","))')'wavelengths=',wavelengths
!!!!do i=1,224
!!!!   print*,i,band_names(i)(:len_trim(band_names(i)))
!!!!end do
!!!!write(*,'(A,7(I3,","))')'bbl=',bbl


!!!!end program envi_test
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
module envi_header

  !This module contains subroutines that read values from ENVI like headers.
  !When called, it means that 'key' has already been found and we also know
  !whether we are looking for an integer, real, string, integer array, real 
  !array, or string array. 
  !
  ! Arrays have the structure:
  ! key = { value(1), value(2), ..., value(n)}
  ! where 1) elements of value ARE separated by commas; 
  !       2) single elements CAN NOT be split across multiple lines
  !       3) '}' can appear on a line by itelf [the line following value(n)],
  !          or at the end of the line with value(n).
  ! 2004-08-31 MJM Made determine_key and match_array_key private
  ! 2004-09-28 MJM Fixed get_data_type in order to be more robust.
  ! 2005-11-29 MJM Made read_envi_string_array public - used in newread_envi_header
  ! 2006-06-07 MJM Now uses default_character_length

  use default_character_length,only:charlen
  implicit none
  private

  type key
     integer:: first_line
     character (len=charlen)::name
     integer::type
     logical::scalar
     integer::arr_size
     integer::int_value
     real::rel_value
     character (len=charlen)::str_value
     integer::name_len, str_len
  end type key
  
  public::match_key !, determine_key, match_array_key !Key matching
  public::search_for_keys, key
  
  public::read_envi_value,read_envi_string_array
  interface read_envi_value ! Use one name to call the following
     module procedure read_envi_integer, read_envi_real,&
          & read_envi_string, read_envi_integer_array,  &
          & read_envi_real_array, read_envi_string_array
  end interface
  
  character (len=15),parameter::file_name='envi_header.f90'
  character (len=11),parameter::module_name='envi_header'

  ! If this is the first character on a line, then the line is ignored...
  ! Safest to comment all lines in an array
  character (len=6),parameter::comment='!$#%*;'
  
  character (len=7),dimension(3),parameter:: type_c=['integer','float  ','string ']
  
contains
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine search_for_keys(found_keys,comb_file,num_found_keys)
    ! Modified 2003-11-05 MJM Don't attempt to determine value if RHS is empty
    !                         do not treat as a found key at all!
    ! Modified 1999-09-23 by MJM
    type(key),dimension(:),intent(inout)::found_keys
    character (len=*), dimension(:),intent(in)::comb_file
    integer,intent(out)::num_found_keys
    
    integer::i,first_rec,out_type,arr_size,j,num_comb_file_lines
    logical::out_scalar
    character (len=charlen)::out_key

    num_comb_file_lines=size(comb_file)
    i=0; j=0
    found_keys(:)%name=''
    search: do 
       i=i+1
       out_key=''
       first_rec=i
       ! if and only if the first character on a line with  a keyword
       ! is one of the comment characters, then that keyword is 
       ! skipped. Because of how arrays are processed, 
       ! the comment character is IGNORED inside arrays. 
       if (scan(comb_file(i),comment).eq.1) then
          if (i.eq.num_comb_file_lines) exit search
          cycle search
       endif
       call determine_key(comb_file,i,out_key,out_type,out_scalar, arr_size)
       if (out_key.ne.''.and.out_type.gt.0) then
          j=j+1
          found_keys(j)%first_line=first_rec
          found_keys(j)%name_len=len_trim(out_key)
          found_keys(j)%name(:found_keys(j)%name_len)=out_key(:found_keys(j)%name_len)
          found_keys(j)%type=out_type
          found_keys(j)%scalar=out_scalar
          found_keys(j)%arr_size=arr_size
          if (out_scalar) then
             select case(out_type)
             case(1) !If integer, could actually be real without decimal; so read twice
                call read_envi_value(comb_file(i),found_keys(j)%int_value) !read as int
                call read_envi_value(comb_file(i),found_keys(j)%rel_value) !read as real
             case(2); call read_envi_value(comb_file(i),found_keys(j)%rel_value)
             case(3)
                call read_envi_value(comb_file(i),found_keys(j)%str_value,found_keys(j)%str_len)
             end select
          end if
       end if
       if (i.eq.num_comb_file_lines) exit search
    end do search
    num_found_keys=j
  end subroutine search_for_keys
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  logical function match_key(junk,key)
    !writen 1999-Sep by MJM
    character (len=*), intent(in)::junk
    character (len=*), intent(in)::key
    
    character (len=len(junk))::hold
    integer::ieq,ilen,ikey
    
    
    match_key=.false.
    ieq=index(junk,'=')
    ikey=len_trim(key)
    if (ieq.gt.0) then
       hold(:ieq-1)=adjustl(junk(:ieq-1))
       ilen=len_trim(hold(:ieq-1))
       if (hold(:ilen).eq.key(:ikey)) match_key=.true.
    end if
    
  end function match_key
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  logical function match_array_key(junk,key)
   !written 1999-Sep by MJM
    character (len=*), intent(in)::junk
    character (len=*), intent(in)::key
    
    character (len=len(junk))::hold
    integer::ieq,ilen,ikey,ibr
    
    
    match_array_key=.false.
    ieq=index(junk,'=')
    ibr=index(junk,'{')
    ikey=len_trim(key)
    if (ieq.gt.0.and.ibr.gt.ieq) then
       hold(:ieq-1)=adjustl(junk(:ieq-1))
       ilen=len_trim(hold(:ieq-1))
       if (hold(:ilen).eq.key(:ikey)) match_array_key=.true.
    end if
    
  end function match_array_key
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine read_envi_integer(junk,value)
    !Modified to be nice to MS Power Fortran 4.o compiler by MJM 1999-10-01
    character (len=*), intent(in)::junk
    integer, intent(out)::value
    
    integer::ieq,ilen
    character (len=len(junk))::junk2
    character (len=5)::fmt
    
    junk2=''; fmt=''
    ieq=index(junk,'=')
    ilen=len_trim(junk)
    junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen)) !mv leading spaces to rt
    ilen=len_trim(junk2)                      !trim trailing spaces
    !MS seems to barf again even though no leading or trailing spaces,
    !so need the following line:
    !Actually, this is overkill: I should use the old form everywhere
    !except when ilen=1, since that is the only time MS has a problem.
    !old form is: read(junk2(1:ilen),*)value
    !when ilen=1: read(junk2(1:ilen),'(I1)')value
    !which works for the reals case, below.
    write(fmt,'("(I",I1,")")')ilen !should handle everything less I10
    read(junk2(1:ilen),fmt)value
    
  end subroutine read_envi_integer
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine read_envi_real(junk,value)
    !written 1999-Sep
    !modified to be nice to MS compilers 1999-10-01
    !Since we sometimes use this to read an int, need
    !fix for MS compiler. -1999-10-21 MJM
    character (len=*), intent(in)::junk
    real, intent(out)::value
    
    integer::ieq,ilen
    character (len=len(junk))::junk2
    
    junk2=''
    ieq=index(junk,'=')
    ilen=len_trim(junk)
    junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen)) !move leading spaces to rt
    ilen=len_trim(junk2)                       !trim trailing spaces
    !Looks like I need to be more careful here. For those times when
    !I want to read an integer and interpret it as a rel, the MS list-directed
    !read on integers less that 10 bug jumps up and bites you in the face. 
    !I can, at this point, construct a correct format statement. 
    !Make ilen=w, then find out where the decimal is from the right, and if 
    !nonzero, make it Fw.d-1 format. I shouldn't need to this, but oh well. 
    !Also, I'll need to test for the presence of an "e" or "E" and the exponent. 
    !Even better: this only bothers us when the value is <10, which means ilen=1;
    !furthermore, it only bothers us it's an integer! So we're set:
    if (ilen.ne.1) then
       read(junk2(1:ilen),*)value
    else
       read(junk2(1:ilen),'(f1.0)')value
    end if

  end subroutine read_envi_real
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine read_envi_string(junk,value,vlen)
    !written 1999-Sep
    !modified to be nice to MS compilers 1999-10-01
    character (len=*), intent(in)::junk
    character (len=*), intent(out)::value
    integer, intent(out)::vlen
    
    integer::ieq,ilen
    character (len=len(junk))::junk2
    
    value=''
    ieq=index(junk,'=')
    ilen=len_trim(junk)
    junk2=''
    junk2(:ilen-ieq)=adjustl(junk(ieq+1:ilen)) !move leading spaces to rt
    ilen=len_trim(junk2)                       !trim trailing spaces
    vlen=ilen
    value(1:ilen)=junk2(1:ilen)

  end subroutine read_envi_string
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine read_envi_integer_array(lun,junk,value,char_array,lastrec)
    !written by MJM 1999-Sep
    !modified to be nice to MS compilers 1999-10-01
    !1999-10-05 More modifications to be even nicer to MS compilers. -MJM
    use errors_mod,only:fatal_error
    integer, intent(in)::lun
    character (len=*), intent(in)::junk
    integer, intent(out), dimension(:)::value
    character (len=*), dimension(:), intent(in), optional::char_array
    integer, intent(out), optional::lastrec
    
    integer::ieq,ilen,ilb,irb,nvals
    integer::fi, li,icom,cfi,cli
    character (len=len(junk))::junk2,junk3
    character (len=23),parameter::subroutine_name='read_envi_integer_array'
    
    ! If char_array is present, it means lun is not lun, but the current
    ! index, i.e. line number of the text file being processed. 
    ! In that case, we don't do a read, we get the next element of the array.
    ! On output, the lastrec is the last index that was used.

    !This is a bit more complicated. We'll get the first line and
    !we'll know that we're expecting an integer array. We will 
    !read until we find a line with '}'. 
    ! How many values are there on a line? If there is not '}', then
    ! the number of integers = number of ','; if there is '}',
    ! then number of integers=number of ','+1; also need to check to 
    ! make sure '{' and '}' might be on line without numbers
    ieq=index(junk,'=')
    ilb=index(junk,'{')
    fi=1
    if (present(lastrec)) lastrec=lun

    junk2=''
    junk3=''
    ilen=len_trim(junk)
    if (ilen.eq.ilb) then !line ends in '{', so read next line
       if (present(char_array)) then
          lastrec=lastrec+1
          junk3=char_array(lastrec)
       else
          read(lun,'(A)')junk3
       end if
       ilen=len_trim(junk3)
       junk2=adjustl(junk3(:ilen))
    else
       junk2(:ilen-ilb)=adjustl(junk(ilb+1:ilen)) !move leading spaces to rt
    end if
    !OK, now we're past '{' go into a loop and see what we find
    
    loop: do
       ilen=len_trim(junk2)                       !trim trailing spaces
       if (ilen.gt.0) then !ilen=0 => line with only '{'
          irb=index(junk2(:ilen),'}')
          if (irb.eq.1.and.ilen.eq.1) exit loop !=> line with only '}'
          nvals=count_occurence(junk2(:ilen),',')
          ! the following is used if last doesn't end in ',' or '}', 
          ! presumably because next line is '}'
          if (verify(junk2(ilen:ilen),'.0123456789').eq.0) nvals=nvals+1
          if (irb.gt.0) then
             nvals=nvals+1 !comma replaced with '}'
             ilen=ilen-1   !'}' is last because of adjustl and len_trim
          end if
          if (junk2(ilen:ilen).eq.',') then
             ilen=ilen-1    !some fortrans don't want dangling comma
          endif
          li=fi-1+nvals
          if (li.gt.size(value)) then
             print*,'ERROR: read_envi_integer_array'
             print*,' Check envi headers: parser finds more values'
             print*,' than size of array pased.'
             print*,' Last line read:',junk2(:ilen)
             call fatal_error(file_name,module_name,subroutine_name,'')
          end if
          !At this point, MS Fortran and SGI Fortran are having 
          !severe disagreements. MS seems to miss the last value for 
          !a list directed read of the character variable, so it needs a format
          !specified. MS seems to work okay with I5 (long enough for all integers
          !we'll use) but SGI F90 doesn't like it. 
          !So, I'll see if I can call read_envi_integer from here
          !for each piece of the string.
!          fmt=''
!          write(fmt,'("(",I2,"I5)",)')nvals
!          read(junk2(:ilen),fmt)value(fi:li)
          cli=-1
          do 
             cfi=cli+2
             icom=index(junk2(cfi:ilen),',')
             if (icom.ne.0) then
                cli=(cfi-1)+icom-1 !cfi-1 is zero for icom search; icom-1 is char before comma
             else
                cli=ilen
             end if
             call read_envi_integer(junk2(cfi:cli),value(fi))
             li=fi !not sure we need this anymore
             if (icom.eq.0) exit !reached end of line
             fi=fi+1 !next index to read
          end do
          if (irb.gt.0) exit loop
          fi=li+1
       endif
       junk3=''
       if (present(char_array)) then
          lastrec=lastrec+1
          junk3=char_array(lastrec)
       else
          read(lun,'(A)')junk3
       end if
       ilen=len_trim(junk3)
       junk2=''
       junk2=adjustl(junk3(:ilen))
    enddo loop

  end subroutine read_envi_integer_array
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine read_envi_real_array(lun,junk,value,char_array,lastrec)
    !written 1999-Sep by MJM
    !modified to be nice to MS compiler 1999-10-01 by MJM
    !Since we sometimes use this to read an int, need
    !fix for MS compiler. -1999-10-21 MJM
    use errors_mod,only:fatal_error
    integer, intent(in)::lun
    character (len=*), intent(in)::junk
    real, intent(out), dimension(:)::value
    character (len=*), dimension(:), intent(in), optional::char_array
    integer, intent(out), optional::lastrec
    
    integer::ieq,ilen,ilb,irb,nvals
    integer::fi, li,icom,cfi,cli
    character (len=len(junk))::junk2,junk3
    character (len=20),parameter::subroutine_name='read_envi_real_array'
    
    !This is a bit more complicated. We'll get the first line and
    !we'll know that we're expecting an real array. We will 
    !read until we find a line with '}'. 
    ! How many values are there on a line? If there is not '}', then
    ! there the number of real = number of ','; if there is '}',
    ! then number of real=number of ','+1; also need to check to 
    ! make sure '{' and '}' might be on line without numbers
    ieq=index(junk,'=')
    ilb=index(junk,'{')
    fi=1
    junk2=''
    junk3=''
    if (present(lastrec)) lastrec=lun

    ilen=len_trim(junk)

    if (ilen.eq.ilb) then !line ends in '{', so read next line
       if (present(char_array)) then
          lastrec=lastrec+1
          junk3=char_array(lastrec)
       else
          read(lun,'(A)')junk3
       end if
       ilen=len_trim(junk3)
       junk2=adjustl(junk3(:ilen))
    else
       junk2(:ilen-ilb)=adjustl(junk(ilb+1:ilen)) !move leading spaces to rt
    end if
    !OK, now we're past '{' go into a loop and see what we find
    loop: do
       ilen=len_trim(junk2)                       !trim trailing spaces
       if (ilen.gt.0) then !ilen=0 => line with only '{'
          irb=index(junk2(:ilen),'}')
          if (irb.eq.1.and.ilen.eq.1) exit loop !=> line with only '}'
          nvals=count_occurence(junk2(:ilen),',')
          ! the following is used if last doesn't end in ',' or '}', 
          ! presumably because next line is '}'
          if (verify(junk2(ilen:ilen),'.0123456789').eq.0) nvals=nvals+1
          if (irb.gt.0) then
             nvals=nvals+1 !comma replaced with '}'
             ilen=ilen-1   !'}' is last because of adjustl and len_trim
          end if
          if (junk2(ilen:ilen).eq.',') then
             ilen=ilen-1    !some fortrans don't want dangling comma
          endif
          li=fi-1+nvals
          if (li.gt.size(value)) then
             print*,'ERROR: read_envi_real_array'
             print*,' Check envi headers: parser finds more values'
             print*,' than size of array pased.'
             print*,' Last line read:',junk2(:ilen)
             call fatal_error(file_name,module_name,subroutine_name,'')
          end if
!!!!          read(junk2(:ilen),*)value(fi:li)
          ! The darn MS list-directed-read-of-integer<=9 bug!
          ! Because of it, I need a kluge. The same one from the 
          ! read_integer_array section above. 
          ! start of fix
          cli=-1
          do 
             cfi=cli+2
             icom=index(junk2(cfi:ilen),',')
             if (icom.ne.0) then
                cli=(cfi-1)+icom-1 !cfi-1 is zero for icom search; icom-1 is char before comma
             else
                cli=ilen
             end if
             call read_envi_real(junk2(cfi:cli),value(fi))
             li=fi !not sure we need this anymore
             if (icom.eq.0) exit !reached end of line
             fi=fi+1 !next index to read
          end do 
          ! end of fix.
          if (irb.gt.0) exit loop
          fi=li+1
       endif
       junk3=''
       if (present(char_array)) then
          lastrec=lastrec+1
          junk3=char_array(lastrec)
       else
          read(lun,'(A)')junk3
       end if
       ilen=len_trim(junk3)
       junk2=''
       junk2=adjustl(junk3(:ilen))
    enddo loop
    
  end subroutine read_envi_real_array
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine read_envi_string_array(lun,junk,value,char_array,lastrec)
    !written 1999-Sep by MJM
    !modified to be nice to MS compilers 1999-10-01 by MJM
    use errors_mod,only:fatal_error
    integer, intent(in)::lun
    character (len=*), intent(in)::junk
    character (len=*), intent(out), dimension(:)::value
    character (len=*), dimension(:), intent(in), optional::char_array
    integer, intent(out), optional::lastrec
    
    integer::ieq,ilen,ilb,irb
    integer::fi,it,ic
    character (len=len(junk))::junk2,junk3
    character (len=21),parameter::subroutine_name='read_envi_string_array'
    
    !This is a bit more complicated. We'll get the first line and
    !we'll know that we're expecting an string array. We will 
    !read until we find a line with '}'. 
    ! How many values are there on a line? If there is not '}', then
    ! there the number of strings = number of ','; if there is '}',
    ! then number of real=number of ','+1; also need to check to 
    ! make sure '{' and '}' might be on line without numbers
    ieq=index(junk,'=')
    ilb=index(junk,'{')
    fi=1
    if (present(lastrec)) lastrec=lun

    junk2=''; junk3=''
    ilen=len_trim(junk)
    if (ilen.eq.ilb) then !line ends in '{'
       if (present(char_array)) then
          lastrec=lastrec+1
          junk3=char_array(lastrec)
       else
          read(lun,'(A)')junk3
       end if
       ilen=len_trim(junk3)
       junk2=adjustl(junk3(:ilen))
    else
       junk2(:ilen-ilb)=adjustl(junk(ilb+1:ilen)) !move leading spaces to rt
    endif
    !OK, now we're past '{' go into a loop and see what we find
    loop: do
       ilen=len_trim(junk2)   !trim trailing spaces
       if (ilen.gt.0) then !ilen=0 => line with only '{'
          irb=index(junk2(:ilen),'}')
          if (irb.eq.1.and.ilen.eq.1) exit loop !=> line with only '}'
          ic=0; it=0
          inner: do 
             it=scan(junk2(ic+1:ilen),',}')
             irb=scan(junk2(ic+1:ilen),'}')
             if (it.eq.0.and.irb.eq.0) it=ilen-ic+1 !no comma or '}'
             if (it.eq.0.and.irb.gt.0) it=irb
             if (fi.gt.size(value)) then
                print*,'ERROR: read_envi_string_array'
                print*,' Check envi headers: parser finds more values'
                print*,' than size of array pased.'
                print*,' Last line read:',junk2(:ilen)
                call fatal_error(file_name,module_name,subroutine_name,'')
             end if
             value(fi)=adjustl(junk2(ic+1:ic+it-1))
             ic=ic+it   !move up start look
             fi=fi+1    !increment index
             if (it.eq.0.and.irb.gt.0) exit inner
             if (ic.ge.ilen) exit inner 
          enddo inner
          if (irb.gt.0) exit loop
       endif
       junk2=''; junk3=''
       if (present(char_array)) then
          lastrec=lastrec+1
          junk3=char_array(lastrec)
       else
          read(lun,'(A)')junk3
       end if
       ilen=len_trim(junk3)
       junk2=adjustl(junk3(:ilen))
    enddo loop
    
  end subroutine read_envi_string_array
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  pure integer function count_occurence(string,set)
    !Function to count number of occurences of _set_ in _string_
    character (len=*),intent(in)::string
    character (len=*),intent(in)::set

    integer::scan_out,lhs
    character (len=len(string))::hold

    
    count_occurence=0
    lhs=1
    hold=string
    do
       scan_out=scan(hold(lhs:),set)     !find position of desired character
       if (scan_out.eq.0) exit           !no more desired chars, exit
       count_occurence=count_occurence+1 !found one, increment count
       lhs=lhs+scan_out                  !start next search just past last found
    end do

  end function count_occurence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function list_all_occurences(string,set) result(locations)
  ! MJM 2010 May 18, for better processing of header files
    character (len=*),intent(in)::string
    character (len=*),intent(in)::set
    integer,dimension(:),allocatable::locations
    
    integer::scan_out,lhs
    character (len=len(string))::hold
    integer::n,i
    
    n=count_occurence(string,set)
!    print*,'n = ',n
    if (n.eq.0) then
      if (allocated(locations)) deallocate(locations)
      allocate(locations(1))
      locations(1)=-1 ! signals error
      return 
    endif
    
    if (allocated(locations)) deallocate(locations)
    allocate(locations(n))
    
    lhs=1
    hold=trim(adjustl(string))
    do i=1,n
       scan_out=scan(hold(lhs:),trim(set))!find position of desired character
       locations(i)=lhs+scan_out-1  !If 1 is returned, then first character in set. 
       lhs=lhs+scan_out
    end do    
   
  end function list_all_occurences
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
!  function split_string(string,set)
!  
!  
!  end function split_string
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  subroutine determine_key(in_char,rec,key,type,scalar,arr_size)
    !Written 1999 Sep by MJM
    !Modified to be nice to MS compiler 1999-10-01 by MJM
    ! 2010 May 18 - check type of each element of an array
    !   if all(types eq 1) then integer
    !   if any(types eq 
    character (len=*), intent(in), dimension(:)::in_char
    integer, intent(inout)::rec
    character (len=*), intent(out)::key
    integer, intent(out)::type
    logical, intent(out)::scalar
    integer, intent(out)::arr_size

    character (len=len(in_char(1)))::junk,rhs,junk2,junk3
    integer:: ieq,key_len, ilbr,ilen,irbr,irhs,icom,inlg,ilen2,i,start_loc,stop_loc,lj2
    integer, dimension(:),allocatable::types ! used to check the type of each element of an array
    integer, dimension(:),allocatable::locations
    
    ! key = value
    ! key = {value, ..., value}
    ! Read in_char and determine if there is a key
    ! if so, then determine scalar or array AND int|real|string.
    ! If array, determine how many elements there are.
    ! On entry, rec is the initial index of the character to check;
    ! on exit, rec is the last index of the character array that has been checked. 
    ! type=1: integer; type=2: real; type=3:string
    ! scalar=true, scalar; scalar=false, array
    ! arr_size=-1 if scalar; the array size otherwise. 

    arr_size=-1
    scalar=.true.
    junk=''
    key=''
    junk=in_char(rec)
    ilen=len_trim(junk)

    ! First, find key name, which is stuff before "="
    ieq=index(junk,'=')
    if (ieq==1) then
       type=-1 ! No key before "="!
       return
    else if (ieq ==0) then
       type=0  ! No "="!
       return
    endif
    key(:ieq-1)=adjustl(junk(:ieq-1))
    key_len=len_trim(key)
    rhs=''; rhs(:ilen-ieq)=adjustl(junk(ieq+1:ilen))

    ! Now, determine key type (scalar or array)
    ilbr=index(junk(ieq+1:ilen),'{')
    irbr=index(junk(ieq+1:ilen),'}')
    if (ilbr == 0 ) then !scalar
       scalar=.true.
    else
       scalar=.false.
    end if
    
    ! Now data type: integer, real, string
    if (scalar) then
       irhs=len_trim(rhs(:ilen-ieq))
       type=get_data_type(rhs(:irhs))
!       write(*,'(5a)')trim(type_c(type)),' : ',trim(key),' = ',trim(adjustl(rhs(:ilen-ieq)))
       return
    else if (.not.scalar.and.irbr>0) then !one line array, so it is short
       rhs=''
       rhs(:irbr-ilbr-1)=adjustl(junk(ieq+ilbr+1:ieq+irbr-1)) !mv sp to rt
       irhs=len_trim(rhs(:irbr-ilbr-1))                       !trunc rt sp
       arr_size=1+count_occurence(rhs(:irhs),',')
       if (arr_size.eq.1) then
          type=get_data_type(rhs(:irhs))
       else
!       icom=index(rhs(:irhs),',') ! get first element
!       if (icom.ne.0) irhs=len_trim(rhs(:icom-1))
          if (allocated(locations)) deallocate(locations)
          allocate(locations(arr_size-1))
          locations=list_all_occurences(rhs(:irhs),',')
          if (allocated(types)) deallocate(types)
          allocate(types(arr_size))
          start_loc=1
          do i=1,arr_size
             if (i.ne.arr_size) then 
                stop_loc=locations(i)-1
             else
                stop_loc=len_trim(rhs)
             endif
             junk2=''
             junk2=trim(adjustl(rhs(start_loc:stop_loc)))
             lj2=len_trim(junk2)
             types(i)=get_data_type(junk2(:lj2))
             if (i.lt.arr_size) start_loc=stop_loc+2 ! get prepared for next: space to after comma
          enddo
!       type=get_data_type(rhs(:irhs))
          if (any(types.eq.3)) then ! if there's one character, they're all characters
             type=3
          elseif (all(types.eq.1)) then ! at this point we're either all ints, ints & reals, or reals
             type=1
          else ! one or more reals (and no chars) means real
             type=2
          endif
          deallocate(types,locations)
       endif
       return
    else !multi-line array
       arr_size=0 !test others, this last case is harder
       !since we're multi-line, are we guaranteed to find something on the 
       !next line? Not necessarily! Next line could be '}'. 
       !So, we'll see what we can learn from this line:
       !If we find a comma, then there's an element on this line:
       icom=index(junk(ieq+ilbr:ilen),',') !see if comma this line
       junk2=''
       !In MS Fortran PwerStation 4.0, with the compiler flags I use, 
       !the following statement produces problems, probably with optimization.
       !The fix follows immediately
       !junk2=adjustl(in_char(rec+1)) !mv ld sp to rt
       junk3=''; junk3=in_char(rec+1); junk2=adjustl(junk3)
       ilen2=len_trim(junk2)         !count to rtmost nonspace
       inlg=verify(junk2(:ilen2),'}') !>0 if more than '}' on line
       !Next line in case array is one elem this line, '}' next line
       if (inlg.eq.0.and.icom.eq.0) icom=ilen-(ieq+ilbr)+1+1 
       if (icom.gt.0) then
          !look at quantity after '{' and before ','
          rhs(:icom-2)=adjustl(junk(ieq+ilbr+1:ieq+ilbr+icom-1-1))
          irhs=len_trim(rhs(:icom-2)) !just first element-no leading or trailing sp
          type=get_data_type(rhs(:irhs))
       else !first line was '{', so use junk2 to determine type
            !No leading or trailing spaces, so:
          irbr=index(junk2(:ilen2),'}')
          if (irbr.gt.0) ilen2=irbr-1 !go to character before '}'
          icom=index(junk2(:ilen2),',') 
          if (icom.gt.0) ilen2=icom-1 !go to character before ','
          irhs=len_trim(junk2(:ilen2)) !trim any uncovered trailing spaces 
          !Now we should have one element w/o leading/trailing spaces
          rhs(:irhs)=junk2(:irhs)
          type=get_data_type(rhs(:irhs))
       end if
       !Now move forward to end of array
       ! and count arr_size: arr_size=number of ',' +1
       ilen=len_trim(junk(ilbr+1:))
       arr_size=count_occurence(junk(ilbr+1:ilbr+ilen),',')
       do_loop:do
          junk=''
          rec=rec+1
          junk=in_char(rec)
          ilen=len_trim(junk)
          arr_size=arr_size+count_occurence(junk(:ilen),',')
          if (index(junk(:ilen),'}').gt.0) exit do_loop
       end do do_loop
       arr_size=arr_size+1 ! '}' not ',' after last element
    end if

  end subroutine determine_key
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
  integer function get_data_type(in_text)
  ! Output value:
  ! type = 0: no value (i.e. empty string was passed
  ! type = 1: integer (may or may not be representable)
  ! type = 2: real    (may or may not be representable)
  ! type = 3: character (may or may not be representable!)
    !written 1999 Sep by MJM
    !1999-Oct-21: Neglected "+-" in integer class definition, so fixed.
    !2003-Nov-05 MJM: Return type 0 when there is no value. This is 
    !                 used to keep from being a key/value pair
    !2004-Sep-28 MJM: Modified a lot in order to use the ability of the 
    !                 program to determine whether or not it can read the
    !                 input as the specified type.
    !2005-Feb-08 MJM: Sigh. Of course, who says compiler writers get it correct? 
    !                 Another modification. Check if strictly character: i.e., 
    !                 no numbers. If so, then automatically classified as character. 
    !                 Some compilers use "E" or "e" as valid reals. I require a number. 
    !                 Sigh.
    !2010-May-18 MJM  Tried to make more foolproof against "long" integers or floats.
    !                 If it cannot be represented as two I*4 and it is int, call it char
    use chartools
    character (len=*), intent(in)::in_text
!o    character (len=12), parameter::integer_class='+-0123456789'
!o    character (len=15), parameter::real_class='0123456789.+-Ee'
!o    integer::ic,rc
    character (len=10),parameter::number_class='0123456789'
    character (len=8)::ifmt, ffmt
    integer::il,ijunk,iios,fios,ichar,idot,ie
    real::fjunk

    get_data_type=0
    il=len(in_text)
    
!o    if (len(in_text).ne.0) then
    ifmt = ''; ffmt = ''
    if (il.ne.0) then

!o       ic=verify(string=in_text,set=integer_class)
!o       if (ic==0) then !integer
!o          get_data_type=1
!o       else 
!o          rc=verify(string=in_text,set=real_class)
!o          if (rc == 0) then !real
!o             get_data_type=2
!o          else  !must be string 
!o             get_data_type=3
!o          end if
!o       end if

       ichar=scan(in_text,number_class)
       if (ichar.eq.0) then ! must be character since there are no numbers
          get_data_type = 3
       else 
          idot=scan(in_text,'.ed') ! Is there a decimal point or d or e? => could be float
          if (il.le.9) then
             write(ifmt,'(A,I1,A)')'(I',il,')'
             write(ffmt,'(A,I1,A)')'(F',il,'.0)'
          else if (il.le.99) then
             write(ifmt,'(A,I2,A)')'(I',il,')'
             write(ffmt,'(A,I2,A)')'(F',il,'.0)'
          else !il ge 100
             write(ifmt,'(A,I3,A)')'(I',il,')'
             write(ffmt,'(A,I3,A)')'(F',il,'.0)'
          end if
          iios=-1 ; fios=-1
          read(in_text,fmt=ifmt,iostat=iios)ijunk
! Here we need discriminate from "this is an integer" vs.
! "this is an integer that can be represented with Fortran 2 byte signed integers"            
          if (iios.eq.0) then ! read correctly as integer, so integer
             get_data_type = 1
             ! we need to focus on representable integers... since I don't have I8 yet, this means
             ! that if il >=10, need to treat as char
             if (il.ge.10) get_data_type=3 ! may never get here since tests 4 byte ints.
          else ! Not an integer
             read(in_text,fmt=ffmt,iostat=fios)fjunk
             ! Need to be careful here; do we want more than 10 digits that looks like integer
             ! to be an integer (fails iios test above) or a float (successfully read as float here).
             ! I don't think we need to do math with such values (yet) but we also don't want to convert
             ! them.....
             if (fios.eq.0) then ! read correctly as real
                get_data_type=2 ! so, real
                if ((il.ge.10).and.(idot.eq.0)) then ! check to see if "long integer"
             ! must be a long integer (fails integer test, passes float test,
             ! and no "." or "e")            
                   get_data_type=3 
                endif
             else ! character string: anything can be a character!
                get_data_type=3
             end if
          end if
       endif
    end if

  end function get_data_type
!!!!+!!!!1!!!!+!!!!2!!!!+!!!!3!!!!+!!!!4!!!!+!!!!5!!!!+!!!!6!!!!+!!!!7!!!!+!!!!8
end module envi_header
