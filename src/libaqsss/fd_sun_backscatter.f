      subroutine fd_sun_backscatter(irad,solar_flux,wind,sun_zenith, tasun_bak_tab, tasun_bak)
      implicit none

      real(4), parameter :: solar_flux_bak= 264. !nominal solar flux for solar backscatter tables

      integer(4), intent(in)  :: irad
      real(4),    intent(in)  :: solar_flux,wind
      real(8),    intent(in)  :: sun_zenith
      real(4),    intent(in)  :: tasun_bak_tab(161,26,3,3)

      real(4),    intent(out) :: tasun_bak(3) 

      real(4) brief,b1,b2,c1,c2
      integer(4) j1,j2,k1,k2,ipol

      if(sun_zenith.ge.90) then
         tasun_bak=0
         return
      endif
      
      brief=wind
      if(brief.lt. 0.00) brief= 0.00
      if(brief.gt.24.99) brief=24.99
      j1=1+brief                                                        
      j2=j1+1                                                           
      b1=j1-brief                                                       
      b2=1-b1

      brief=5*(sun_zenith-58.) 
      if(brief.lt.  0.00) brief=  0.00
      if(brief.gt.159.99) brief=159.99
      k1=1+brief                                                       
      k2=k1+1                                                           
      c1=k1-brief                                                       
      c2=1-c1                                                          
        
      do ipol=1,3
         tasun_bak(ipol)=b1*(c1*tasun_bak_tab(k1,j1,ipol,irad)+c2*tasun_bak_tab(k2,j1,ipol,irad))+
     &        b2*(c1*tasun_bak_tab(k1,j2,ipol,irad)+c2*tasun_bak_tab(k2,j2,ipol,irad))   
         tasun_bak(ipol)=(solar_flux/solar_flux_bak)*tasun_bak(ipol)
      enddo     
                
      return
        end
