      subroutine fd_climate_land_refl(irad,idayjl,isecdy,xlat,xlon, reflv_land_tab, reflh_land_tab, refl)                            
        implicit none

        integer(4), intent(in)  :: irad,idayjl,isecdy
        real(4),    intent(in)  :: xlat,xlon
        real(4),    intent(in)  :: reflv_land_tab(360,181,12,3),reflh_land_tab(360,181,12,3)  

        real(4),    intent(out) :: refl(2)
                                                                      
        integer(4) i1,i2,j1,j2,k1,k2

        real(4) a1,a2,b1,b2,c1,c2,brief
        real(8) timeyr
        
c     check inputs                                            
          
      if(idayjl.lt.0 .or. idayjl.gt.366)   stop 'idayjl oob in fd_climate_land_refl, pgm stopped'                                                        
      if(isecdy.lt.0 .or. isecdy.gt.86400) stop 'idayjl oob in fd_climate_land_refl, pgm stopped'                                                        
      if(abs(xlat).gt.90.)                 stop 'xlat   oob in fd_climate_land_refl, pgm stopped'
      if(xlon.lt.0. .or. xlon.gt.360.)     stop 'xlon   oob in fd_climate_land_refl, pgm stopped'

c     do time,lat,lon interpolation,  2629800 is the average num sec in a month, 86400*365.25/12        

      timeyr=86400.d0*(idayjl-1) + isecdy
      brief=(timeyr-1314900.d0)/2629800.d0                                  
      i1=1+brief                                                        
      i2=i1+1                                                           
      a1=i1-brief                                                       
      a2=1-a1                                                          
      if(i1.eq. 0) i1=12
      if(i2.eq.13) i2= 1

!     following is ncep lat/lon grid

      brief=90.-xlat
        if(brief.gt.179.999) brief=179.999
      j1=int(1+brief)
      j2=j1+1
      b1=j1-brief
      b2=1.-b1

      brief=xlon
      if(brief.gt.359.999) brief=0.
      k1=int(1+brief)
      k2=k1+1
      if(k2.eq.361) k2=1
      c1=k1-brief
      c2=1.-c1

      refl(1) = a1*b1*(c1*reflv_land_tab(k1,j1,i1,irad)+c2*reflv_land_tab(k2,j1,i1,irad)) +
     &          a1*b2*(c1*reflv_land_tab(k1,j2,i1,irad)+c2*reflv_land_tab(k2,j2,i1,irad)) +                 
     &          a2*b1*(c1*reflv_land_tab(k1,j1,i2,irad)+c2*reflv_land_tab(k2,j1,i2,irad)) +
     &          a2*b2*(c1*reflv_land_tab(k1,j2,i2,irad)+c2*reflv_land_tab(k2,j2,i2,irad))                  
                                                                       
      refl(2) = a1*b1*(c1*reflh_land_tab(k1,j1,i1,irad)+c2*reflh_land_tab(k2,j1,i1,irad)) +
     &          a1*b2*(c1*reflh_land_tab(k1,j2,i1,irad)+c2*reflh_land_tab(k2,j2,i1,irad)) +                 
     &          a2*b1*(c1*reflh_land_tab(k1,j1,i2,irad)+c2*reflh_land_tab(k2,j1,i2,irad)) +
     &          a2*b2*(c1*reflh_land_tab(k1,j2,i2,irad)+c2*reflh_land_tab(k2,j2,i2,irad))  
     
      return                                                            
      end

