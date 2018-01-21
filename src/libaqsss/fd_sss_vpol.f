      subroutine fd_sss_vpol(tht,sst,win,surtb, 
     1                       dtb_dwin_tab, ireflecv, ireflech, sss)
      implicit none

      integer(4), parameter :: n_rad=3  

      real(4),    intent(in )  :: tht,sst,win,surtb(2)
      real(4),    intent(in) :: dtb_dwin_tab(26,2,n_rad) 
      integer(2), intent(in) :: ireflecv(226,371,901), ireflech(226,371,901)
      real(4),    intent(out)  :: sss

      real(4),    parameter :: missing_val=-9999. ! indicator for missing value

      integer(4) iter
      real(4) surtep,dems(2),sss0,sss1,dsss,refl(2)
      real(4) a,y(2),y0(2),y1(2)
        
      data dsss/0.5/ !sss increment
        
      surtep=273.15+sst
        
      call wentz_dems(tht,sst,win, dems, dtb_dwin_tab) 
      y=surtb - surtep*dems !this is estimate of specular emission
        
c     solve for salinity using newtons method

      sss0=35
      sss1=sss0 + dsss
        
      do iter=1,50

!     0. is wind, -999 means ignore phir        
         call fd_water_refl(tht,sss0,sst,0.,-999., 
     1        ireflecv, ireflech, dtb_dwin_tab, refl) 
         y0=surtep*(1-refl)
         call fd_water_refl(tht,sss1,sst,0.,-999., 
     1        ireflecv, ireflech, dtb_dwin_tab, refl) 
         y1=surtep*(1-refl)
        
         a=(y1(1)-y0(1))/dsss   !slope

         ! JMG  02/24/12
         if (abs(a).lt.1.0E-8) then
            sss=missing_val
            return
         endif

         sss=sss0 + (y(1)-y0(1))/a
         if(sss.lt.0) sss=0
        
         if(abs(sss-sss0).lt.0.0005) exit
        
         sss0=sss
         sss1=sss0 + dsss  
      enddo  !iter
        
      return
      end

