subroutine fd_sss_clm(itime,xlat,xlon, staticdata, sss) 

  use static_data_module

  implicit none                                 
      
  integer(4), intent(in) :: itime ! seconds of year
  real(4), intent(in)    :: xlon, xlat
  type(static_data_struct), intent(in)    :: staticdata
  real(4), intent(out)   :: sss

  real(4) xlon0
  real(4)          :: brief, a1, a2, b1, b2, c1, c2
  integer(4)       :: i1, i2, j1, j2, k1, k2

!                                                                       
!     check inputs                                                      
!                                                                       

  if(itime.lt.0 .or. itime.gt.31622400) stop 'error2 in fdreynold, pgm stopped'  

  xlon0 = xlon
  if(xlon.lt.0. .or. xlon.gt.360.) then
     xlon0 = 360 + xlon
  endif

  if(abs(xlat).gt.90.)  stop 'error4 in fdreynold, pgm stopped'                         
!                                                                       
!     do time,lat,lon interpolation                                    
!                                                                       
!     2629800 is the average num sec in a month, 86400*365.25/12   

  brief=(itime-1314900)/2629800.d0                                  
  if(brief.lt.0) then                                               
     i1=12                                                             
     i2=1                                                              
     a1=-brief                                                         
     a2=1.-a1                                                          
  else                                                              
     i1=1+brief                                                        
     i2=i1+1                                                           
     if(i2.eq.13) i2=1                                                 
     a1=i1-brief                                                       
     a2=1.-a1                                                          
  endif
      
  brief=xlat+89.5
  j1=1+brief                                                        
  j2=j1+1                                                           
  b1=j1-brief                                                       
  b2=1-b1
  if(j1.eq.  0) j1=  1
  if(j2.eq.181) j2=180
  !
  brief=xlon0-0.5 
  k1=1+brief                                                       
  k2=k1+1                                                           
  c1=k1-brief                                                       
  c2=1-c1                                                          
  if(k1.eq.  0) k1=360
  if(k2.eq.361) k2=  1 
!                                                                       
  sss=                                         &
       a1*b1*(c1*staticdata%climate_sal(k1,j1,i1)+c2*staticdata%climate_sal(k2,j1,i1))+ & 
       a1*b2*(c1*staticdata%climate_sal(k1,j2,i1)+c2*staticdata%climate_sal(k2,j2,i1))+ &
       a2*b1*(c1*staticdata%climate_sal(k1,j1,i2)+c2*staticdata%climate_sal(k2,j1,i2))+ &
       a2*b2*(c1*staticdata%climate_sal(k1,j2,i2)+c2*staticdata%climate_sal(k2,j2,i2)) 
  
  sss=0.002*sss   !convert to parts/thousand
                                                                 
  return                                                            
end subroutine fd_sss_clm
