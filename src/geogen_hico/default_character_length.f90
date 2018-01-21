module default_character_length
  !2006 June 7 MJM Define default length characters for filenames, etc. 
  !                This is temporary. Fortran 2003 has the ability to use 
  !                any length character strings. Since I am not sure how many
  !                compilers support it, I'll just use this module to define a
  !                large default character length.
  
  implicit none
  private
  
  integer, parameter, public:: charlen= 1024
  
  
end module default_character_length
