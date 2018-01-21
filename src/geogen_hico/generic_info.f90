module generic_info
!2008 August, Marcos J. Montes
! Naval Research Laboratory 
! 
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

   use default_character_length,only:charlen
   implicit none
   private
   public::get_username,get_username_s,get_hostname,get_hostname_s
   
! Fixed parameters for "hostname" routines
! One or more are commonly used in windows and various unix shells
   integer, parameter::h_tries=3
   character(len=charlen),dimension(h_tries),parameter::host=['COMPUTERNAME', &
        'HOST        ','HOSTNAME    ']

! Fixed parameters for "username" routines
! One or more are commonly used in windows and various unix shells
   integer, parameter::u_tries=3
   character(len=charlen),dimension(u_tries),parameter::user=['USERNAME','USER    ','LOGNAME ']


contains
! Note: If a character function, then need to set length at top.
! If character subroutine, then I would not need to set the length:
! The length would be whatever the length of the char that is passed in.
   subroutine get_username_s(username,len_out,stat_out)
!2008 August, Marcos J. Montes
! Naval Research Laboratory 
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
      implicit none
      character (len=*),intent(out)::username ! length from var in caller
      integer, intent(out),optional::len_out,stat_out
      integer::len_o,i,stat,lu

      lu=len(username)
   
      do i=1,u_tries
         username=''
         len_o=0
         call get_environment_variable(user(i),value=username,length=len_o,&
              trim_name=.true.,status=stat)
         if (stat.le.0) exit
      enddo
      if (len_o.eq.0) then
         username=''
         if (lu.ge.25) username='<USERNAME_IS_UNAVAILABLE>'
      endif

      if (present(len_out)) len_out=len_o
      if (present(stat_out)) stat_out=stat

   end subroutine get_username_s

   function get_username() result(username)
! 2008 MJM Naval Research Laboratory, Washington, DC
! A function to get the username whether Windows or Linux 
! May depend on shell
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
      implicit none
      character (len=charlen)::username
      integer::len,i,stat
   
      do i=1,u_tries
         username=''
         len=0
         call get_environment_variable(user(i),value=username,length=len,&
              trim_name=.true.,status=stat)
         if (stat.le.0) exit
      enddo
      if (len.eq.0) then
         username=''
         username='<USERNAME_IS_UNAVAILABLE>'
      endif

            
   end function get_username

   subroutine get_hostname_s(hostname,len_out,stat_out)
!2008 August, Marcos J. Montes
! Naval Research Laboratory 
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
      implicit none
      character (len=*),intent(out)::hostname
      integer, intent(out),optional::len_out,stat_out
      integer::len_o,i,stat,lu

      lu=len(hostname)
   
      do i=1,h_tries
         hostname=''
         len_o=0
         call get_environment_variable(host(i),value=hostname,length=len_o,&
              trim_name=.true.,status=stat)
         if (stat.le.0) exit
      enddo
      if (len_o.eq.0) then
         hostname=''
         if (lu.ge.25) hostname='<USERNAME_IS_UNAVAILABLE>'
      endif

      if (present(len_out)) len_out=len_o
      if (present(stat_out)) stat_out=stat

   end subroutine get_hostname_s

   function get_hostname() result(hostname)
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
! 2008 MJM Naval Research Laboratory, Washington, DC
! A function to get the username whether Windows or Linux 
! May depend on shell
      implicit none
      character (len=charlen)::hostname
      integer::len,i,stat

      do i=1,h_tries
         hostname=''
         len=0
         call get_environment_variable(host(i),value=hostname,length=len,&
              trim_name=.true.,status=stat)
         if (stat.le.0) exit
      enddo 
      if (len.eq.0) then
         hostname=''
         hostname='<HOSTNAME_IS_UNAVAILABLE>'
      endif

                  
   end function get_hostname


end module generic_info
