      subroutine fd_sss(tht,sst,win,surtb, sss_coef, dummy1, dummy2, sss)
      implicit none

      real(4),    intent(in )  :: tht,sst,win,surtb(2)
      real(4),    intent(in )  :: sss_coef(4,251,451)
      integer(2), intent(in )  :: dummy1, dummy2
      real(4),    intent(out)  :: sss

      real(4) brief,a1,a2,b1,b2
      real(4) sss_coef1,sss_coef2,sss_coef3,sss_coef4
      integer(4) i1,i2,j1,j2

      ! JMG  10/12/11
      if(tht.lt.25 .or. tht.gt.49.999) then
         write(*,*)'pgm stopped, tht oob in sss_algo' 
         return
      endif

!      if(sst.lt.-5 .or. sst.gt.39.999) stop 'pgm stopped, sst oob in sss_algo'

      brief=10*(tht-25.)
      i1=int(1+brief)                                                        
      i2=i1+1 
      a1=i1-brief                                                       
      a2=1.-a1

        
      brief=10*(sst+5.)         
      j1=int(1+brief)                                                        
      j2=j1+1 
      b1=j1-brief                                                       
      b2=1.-b1

      sss_coef1=a1*b1*sss_coef(1,i1,j1)       + a1*b2*sss_coef(1,i1,j2) + a2*b1*sss_coef(1,i2,j1) + a2*b2*sss_coef(1,i2,j2)
      sss_coef2=a1*b1*sss_coef(2,i1,j1)       + a1*b2*sss_coef(2,i1,j2) + a2*b1*sss_coef(2,i2,j1) + a2*b2*sss_coef(2,i2,j2)
      sss_coef3=a1*b1*sss_coef(3,i1,j1)       + a1*b2*sss_coef(3,i1,j2) + a2*b1*sss_coef(3,i2,j1) + a2*b2*sss_coef(3,i2,j2)
      sss_coef4=a1*b1*sss_coef(4,i1,j1)       + a1*b2*sss_coef(4,i1,j2) + a2*b1*sss_coef(4,i2,j1) + a2*b2*sss_coef(4,i2,j2)

      sss=sss_coef1 + sss_coef2*surtb(1) + sss_coef3*surtb(2) + sss_coef4*win

      return
      end
