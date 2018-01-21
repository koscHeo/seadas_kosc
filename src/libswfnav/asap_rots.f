        subroutine asap_rots(iyinit,idinit,tsap,asap,nstp,vecs)
c
c  Purpose:  This subroutine rotates orbit vectors from inertial to ECEF
c               coordinates.
c
c  Calling Arguments:
c
c  Name         Type    I/O     Description
c  --------     ----    ---     -----------
c  iyinit       I*4      I      Time tag year
c  idinit       I*4      I      Time tag day
c  tsap(nstp)   R*8      I      Time tag seconds-of-day
c  asap(6,nstp) R*8      I      Inertial orbit vectors
c  nstp         I*4      I      Number of input vectors
c  vecs(6,nstp) R*8      I      ECEF orbit vectors
c
c
c
c  By: Frederick S. Patt, GSC, December 21, 1993
c
c  Notes:  
c
c  Modification History:
c

      implicit none

#include "nav_cnst.fin"

      integer*4 nstp, iyinit, idinit, i
      real*8 asap(6,nstp), vecs(6,nstp), tsap(nstp)
      real*8 gha, cogha, sigha, td

      do i=1,nstp
        td = idinit + tsap(i)/864.d2
        call gha2000(iyinit,td,gha)
        cogha = cos(gha/radeg)
        sigha = sin(gha/radeg)
        vecs(1,i) = asap(1,i)*cogha + asap(2,i)*sigha
        vecs(2,i) = asap(2,i)*cogha - asap(1,i)*sigha
        vecs(3,i) = asap(3,i)
        vecs(4,i) = asap(4,i)*cogha + asap(5,i)*sigha + omegae*vecs(2,i)
        vecs(5,i) = asap(5,i)*cogha - asap(4,i)*sigha - omegae*vecs(1,i)
        vecs(6,i) = asap(6,i)

      end do
      return
      end
