      subroutine fd_climate_sss(idayjl,isecdy,xlat,xlon, sss_tab, sss)
        implicit none

        integer(4), intent(in)  :: idayjl,isecdy
        real(4),    intent(in)  :: xlat,xlon
        real(4),    intent(in)  :: sss_tab(360,180,12)

        real(4),    intent(out) :: sss
                                                                      
        integer(4) i1,i2,j1,j2,k1,k2
        real(4) a1,a2,b1,b2,c1,c2,brief
        real(8) timeyr
        
c     check inputs                                            
          
      if(idayjl.lt.0 .or. idayjl.gt.366)   stop 'idayjl oob in fd_climate_sss, pgm stopped'                                                      
      if(isecdy.lt.0 .or. isecdy.gt.86400) stop 'idayjl oob in fd_climate_sss, pgm stopped'                                                      
      if(abs(xlat).gt.90.)                 stop 'xlat   oob in fd_climate_sss, pgm stopped'
      if(xlon.lt.0. .or. xlon.gt.360.)     stop 'xlon   oob in fd_climate_sss, pgm stopped'

c     do time,lat,lon interpolation,  2629800 is the average num sec in a month, 86400*365.25/12        

      timeyr=86400.d0*(idayjl-1) + isecdy
      brief=(timeyr-1314900.d0)/2629800.d0                                  
      i1=1+brief                                                        
      i2=i1+1                                                           
      a1=i1-brief                                                       
      a2=1-a1                                                          
      if(i1.eq. 0) i1=12
      if(i2.eq.13) i2= 1

      brief=xlat+89.5
      j1=1+brief                                                        
      j2=j1+1                                                           
      b1=j1-brief                                                       
      b2=1-b1
      if(j1.eq.  0) j1=  1
      if(j2.eq.181) j2=180

      brief=xlon-0.5
      k1=1+brief                                                       
      k2=k1+1                                                           
      c1=k1-brief                                                       
      c2=1-c1                                                          
      if(k1.eq.  0) k1=360
      if(k2.eq.361) k2=  1 
                                                                       
      sss= a1*b1*(c1*sss_tab(k1,j1,i1)+c2*sss_tab(k2,j1,i1)) + a1*b2*(c1*sss_tab(k1,j2,i1)+c2*sss_tab(k2,j2,i1)) +                 
     &     a2*b1*(c1*sss_tab(k1,j1,i2)+c2*sss_tab(k2,j1,i2)) + a2*b2*(c1*sss_tab(k1,j2,i2)+c2*sss_tab(k2,j2,i2))                  
                                                                       

      return                                                            
      end

