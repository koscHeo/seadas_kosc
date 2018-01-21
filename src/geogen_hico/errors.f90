module errors_mod

implicit none
private
public:: fatal_error, open_error

contains
!+++*++++1++++*++++2++++*++++3++++*++++4++++*++++5++++*++++6++++*++++7++++*++++8
  subroutine fatal_error(file_name,module_name, subroutine_name, error_message)
    ! 2000-06-17 MJM: Slight change in output text.
    ! 2000-04-17 MJM Creation date. 
    ! Purpose: provide a generic "stop point". Especially necessary on windows
    ! machines. For example, when spawned by ENVI, if an error is encountered, 
    ! the program - and the DOS window- disappear without one being able to 
    ! read the error message. The intent is that this program will allow the 
    ! process to not die until input is received.
    implicit none
    character (len=*), intent(in)::file_name
    character (len=*), intent(in):: module_name, subroutine_name, error_message
    character (len=1)::ans

    write(*,*)'A FATAL ERROR HAS BEEN DETECTED'
    write(*,*)' Source code file: ',file_name
    write(*,*)' Module in file where error is found: ',module_name
    write(*,*)' Subroutine in module where error is found: ',subroutine_name
    write(*,*)' Text of error message: BEGIN'
    write(*,*)error_message
    write(*,*)' END'
    write(*,*)'STOPPING'
    write(*,*)''
    write(*,*)'Please copy the above messages'
    write(*,*)'Depending on the implementation and how this program was called,'
    write(*,*)'the window may disappear following your response to the'//&
         &' question below.'
    write(*,*)''
    write(*,*)'To stop, ENTER any character'
    read(*,*)ans
    stop

  end subroutine fatal_error
!+++*++++1++++*++++2++++*++++3++++*++++4++++*++++5++++*++++6++++*++++7++++*++++8
  subroutine open_error(file_name,module_name,subroutine_name,unit_num,&
       infile_name,ios)
    implicit none
    character (len=*), intent(in)::file_name
    character (len=*), intent(in):: module_name, subroutine_name, infile_name
    integer, intent(in)::ios,unit_num
    character (len=1)::ans
    
    write(*,*)'A FILE COULD NOT BE OPENED'
    write(*,*)' Source code file: ',file_name
    write(*,*)' Module in file where error is found: ',module_name
    write(*,*)' Subroutine in module where error is found: ',subroutine_name
    write(*,*)' Attempted to open unit number = ',unit_num
    write(*,*)' Attempted to open a file named = ',infile_name
    write(*,*)' Value of IOSTAT = ',ios
    write(*,*)' Please ensure the above file exists. '
    write(*,*)'STOPPING'
    write(*,*)''
    write(*,*)'Please copy the above messages'
    write(*,*)'Depending on the implementation and how this program was called,'
    write(*,*)'the window may disappear following your response to the'//&
         &' question below.'
    write(*,*)''
    write(*,*)'To stop, ENTER any character'
    read(*,*)ans
    stop
   

  end subroutine open_error
!+++*++++1++++*++++2++++*++++3++++*++++4++++*++++5++++*++++6++++*++++7++++*++++8
end module errors_mod
