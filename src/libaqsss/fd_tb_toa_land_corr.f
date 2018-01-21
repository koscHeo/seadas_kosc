      subroutine fd_tb_toa_land_corr(irad, xmon, zang, sclon, dtb, tb_landcorr_vpol, tb_landcorr_hpol)

      implicit none

      integer(4), intent(in) :: irad
      real(8),    intent(in) :: xmon
      real(8),    intent(in) :: zang,sclon
      integer(2), intent(in) :: dtb(1440,1440,12,2,3) 

      real(4),    intent(out) :: tb_landcorr_vpol, tb_landcorr_hpol

      integer(4) i1,i2,j1,j2,k1,k2
      real(4) a1,a2,b1,b2,c1,c2,brief,dtb1,dtb2

!      brief=int(xmon-0.5)
      brief=xmon-0.5 ! Fixed per F.Wentz email 06/27/11
      i1=1+brief
      i2=i1+1
      a1=i1-brief
      a2=1.-a1
      if(i1.eq. 0) i1=12
      if(i2.eq.13) i2= 1

      brief=4*(zang-0.125)
      j1=int(1+brief)
      j2=j1+1
      b1=j1-brief
      b2=1.-b1
      if(j1.eq.   0) j1=1440
      if(j2.eq.1441) j2=   1
 
      brief=4*(sclon-0.125)
      k1=int(1+brief)
      k2=k1+1
      c1=k1-brief
      c2=1.-c1
      if(k1.eq.  0) k1=1440
      if(k2.eq.1441) k2=  1 

      dtb1=b1*(c1*dtb(k1,j1,i1,1,irad)+c2*dtb(k2,j1,i1,1,irad)) + b2*(c1*dtb(k1,j2,i1,1,irad)+c2*dtb(k2,j2,i1,1,irad))
      dtb2=b1*(c1*dtb(k1,j1,i2,1,irad)+c2*dtb(k2,j1,i2,1,irad)) + b2*(c1*dtb(k1,j2,i2,1,irad)+c2*dtb(k2,j2,i2,1,irad))
      tb_landcorr_vpol = 0.01*(a1*dtb1 + a2*dtb2) !land comtamination

      dtb1=b1*(c1*dtb(k1,j1,i1,2,irad)+c2*dtb(k2,j1,i1,2,irad)) + b2*(c1*dtb(k1,j2,i1,2,irad)+c2*dtb(k2,j2,i1,2,irad))
      dtb2=b1*(c1*dtb(k1,j1,i2,2,irad)+c2*dtb(k2,j1,i2,2,irad)) + b2*(c1*dtb(k1,j2,i2,2,irad)+c2*dtb(k2,j2,i2,2,irad))
      tb_landcorr_hpol = 0.01*(a1*dtb1 + a2*dtb2) !land comtamination]

      return
      end

