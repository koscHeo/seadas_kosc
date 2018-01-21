!     this collection of routines is the same as contained in the file 'aquarius_geolocation_v5.f' we sent previously,
!     with the following exceptions
!     1.  use module name change from 'sim_orbit5_module' to 'constants_module'
!     2.  new way to find boresight vector b(3) for icase.eq.1
!     3.  subroutine fd_sunr added at end of this file
!     4.  moon glitter is added
!     5.  all trigometric functions use generic sind, cosd, etc rather than real(8) dsind, dcosd, etc
!        (you may have to change this for your compiler)
!     this is the rss v7 routine dated june 6 2009


!     x, y, z  s/c vectors
!     x vector points in general direction of s/c velocity
!     z vector points is s/c nadir and points towards earth for normal attitude
!     y vector is given by z cross x                                         

        subroutine fd_cell_parameters(icase,rotangle,r,ru,x,y,z,beta,alpha,sun,moon,  b,
     &                              ref,xlat,xlon_rot,range,thtinc,azim,thtinc_sun,azim_sun,sunang,moonang)
!       use constants_module
        implicit none

        real(8), parameter :: re=6378.137d3  !earth equatorial radius (meters)
        real(8), parameter :: rp=6356.752d3 !earth      polar radius (meters)
        real(8), parameter :: ffac=(rp/re)*(rp/re)

c       earth dimensions match wgs84          

        integer(4), intent(in) :: icase
        real(8), intent(in)    :: rotangle,r,ru(3),x(3),y(3),z(3),beta,alpha,sun(3),moon(3)
        real(8), intent(inout) :: b(3)
        real(8), intent(out)   :: ref(3),xlat,xlon_rot,range,thtinc,azim,thtinc_sun,azim_sun,sunang,moonang

        integer(4) ierror                                                                                                 
        real(8) rearth,sinlat,coslat,costht,xlon
        real(8) sinbeta,cosbeta
        real(8) bx,by,bz,cel1,cel2,cel3
        real(8) coslon,sinlon
        real(8) geoid(3),cosazim,cosang
        real(8) deltau

        real(8) dcosd, dsind, datand, datan2d, dacosd

        if(icase.eq.1) then
           sinlat=ru(3)
           coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))
           rearth=rp/sqrt(ffac*coslat**2+sinlat**2)  
           deltau=((r-rearth)/re)*sqrt(1 - (1-ffac)*sinlat**2/((ffac*coslat)**2+sinlat**2))
           b(1)=-ru(1)
           b(2)=-ru(2)
           b(3)=-ru(3)*(1+deltau)/(ffac+deltau)
           b=b/sqrt(dot_product(b,b))      
        endif

        if(icase.eq.2) then
           cosbeta=dcosd(beta)                                                                      
           sinbeta=dsind(beta)                                                                      
           bx=sinbeta*dcosd(90-alpha) !input alpha defined relative to the y axes,  alpha=0 denotes bx=0, by=1
           by=sinbeta*dsind(90-alpha) !input alpha defined relative to the y axes , alpha=0 denotes bx=0, by=1   
           bz=cosbeta                                 
           b=bx*x+by*y+bz*z                            
        endif
        

!     for icase.eq.3, boresight vector b is an input 

        call loccel_range(r,ru(1),ru(2),ru(3),b(1),b(2),b(3),re,rp, cel1,cel2,cel3,rearth,range,ierror)

! change to warning JMG 08/07/11

        if(ierror.ne.0 .or. range.le.0) write(*,*) 'obs off of earth in fd_cell_parameters_2, pmg stopped'

! if OOB then return -999.0  JMG  10/12/11
        if(ierror.eq.0 .and. range.gt.0) then
           coslat=sqrt(cel1*cel1+cel2*cel2) !cosine of geocentric lat                                 
           xlat=datand(cel3/(ffac*coslat)) !geodetic cell latitude
           xlon=datan2d(cel2,cel1)                                              
           if(xlon.lt.  0) xlon=xlon+360                                  
           if(xlon.ge.360) xlon=xlon-360
        else
           xlat = -999.0
           xlon_rot = -999.0
           thtinc = -999.0
           return
        endif

        xlon_rot=xlon - rotangle !convert from inertial to earth
        if(xlon_rot.lt.  0) xlon_rot=xlon_rot+360                                  
        if(xlon_rot.ge.360) xlon_rot=xlon_rot-360

        coslon=dcosd(xlon)
        sinlon=dsind(xlon)
        coslat=dcosd(xlat)      !cosine of geodetic lat
        geoid(1)=coslat*coslon                                
        geoid(2)=coslat*sinlon                                
        geoid(3)=dsind(xlat)

        costht=-dot_product(b,geoid)
        if(costht.lt.-1.) costht=-1.
        if(costht.gt. 1.) costht= 1.
        thtinc=dacosd(costht)
      
        if(abs(costht).gt.0.9999999999d0) then !to close to nadir to get azim
           azim=0 
        else                         
           cosazim=(-sinlon*b(1)+coslon*b(2))/sqrt(1.-costht*costht)
           if(cosazim.lt.-1.) cosazim=-1.
           if(cosazim.gt. 1.) cosazim= 1.
           azim=dacosd(cosazim)
           if(b(3)+costht*geoid(3).lt.0.) azim=-azim
           azim=90-azim         ! convert to relative to clockwise from north
           if(azim.lt.  0) azim=azim+360
           if(azim.gt.360) azim=azim-360
        endif
        
        
        ref=b + 2*costht*geoid  !specular reflected boresight ray

!     compute sun glitter angle
        cosang=dot_product(ref,sun)
        if(cosang.gt. 1) cosang= 1
        if(cosang.lt.-1) cosang=-1
        sunang=dacosd(cosang)

!     compute moon glitter angle
        cosang=dot_product(ref,moon)
        if(cosang.gt. 1) cosang= 1
        if(cosang.lt.-1) cosang=-1
        moonang=dacosd(cosang)

!     compute sun tht and azimuth angles
        costht=dot_product(sun,geoid) 
        if(costht.lt.-1.) costht=-1.
        if(costht.gt. 1.) costht= 1.
        thtinc_sun=dacosd(costht)

        if(abs(costht).gt.0.9999999999d0) then !to close to nadir to get azim_sun
           azim_sun=0 
        else                         
           cosazim=(-sinlon*sun(1)+coslon*sun(2))/sqrt(1.-costht*costht)
           if(cosazim.lt.-1.) cosazim=-1.
           if(cosazim.gt. 1.) cosazim= 1.
           azim_sun=dacosd(cosazim)
           if(sun(3)-costht*geoid(3).lt.0.) azim_sun=-azim_sun !minus sign used here becuase sun is pointing away from surface
           azim_sun=90-azim_sun ! convert to relative to clockwise from north
           if(azim_sun.lt.  0) azim_sun=azim_sun+360
           if(azim_sun.gt.360) azim_sun=azim_sun-360
        endif

        return
        end









      subroutine loccel_range(r,r1,r2,r3,s1,s2,s3,re,rp, cel1,cel2,cel3,rearth,range,ierror)     
        implicit none

        real(8), intent(in)  :: r,r1,r2,r3,s1,s2,s3
        real(8), intent(out) :: cel1,cel2,cel3,rearth,range
        real(8) delta,a,cosx,b,c,celmag,arg
        real(8) re,rp
        integer(4) ierror

        
        delta=(re/rp)*(re/rp)-1.                                               
        a=1.+delta*s3*s3                                                  
        cosx=-r1*s1-r2*s2-r3*s3
        b=cosx-delta*r3*s3
        c=1.-(re/r)*(re/r)+delta*r3*r3

        arg=b*b-a*c
        if(arg.lt.0) then !vector does not intercept earth
        ierror=1
        return
        else
        ierror=0
        endif

      range=(b-sqrt(arg))/a                                         
      cel1=r1+range*s1                                                  
      cel2=r2+range*s2                                                  
      cel3=r3+range*s3                                                  
      celmag=sqrt(cel1*cel1+cel2*cel2+cel3*cel3)                              
      cel1=cel1/celmag                                                  
      cel2=cel2/celmag                                                  
      cel3=cel3/celmag                                                  
      rearth=r*celmag                                                   
      range=range*r
      return                                                            
      end                                                               





        subroutine find_geoid(scpos, geoid,rearth)      !geoid points away from earth like the surface normal vector
!       use constants_module
        implicit none

        real(8), parameter :: re=6378.137d3  !earth equatorial radius (meters)
        real(8), parameter :: rp=6356.752d3 !earth      polar radius (meters)
        real(8), parameter :: ffac=(rp/re)*(rp/re)

        real(8), intent(in)  :: scpos(3)
        real(8), intent(out) :: geoid(3),rearth

        real(8) r,ru(3),sinlat,coslat
        real(8) sindel,cosdel,xsinlat,xcoslat,ratio

        r=sqrt(dot_product(scpos,scpos))
        ru=scpos/r

        sinlat=ru(3)                                                        
        coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))    
        rearth=re*rp/sqrt((rp*coslat)*(rp*coslat)+(re*sinlat)*(re*sinlat))                  
        sindel=sin((rearth/r)*(1.003*(1-ffac)*coslat*sinlat))                                    
        cosdel=sqrt(1.-sindel*sindel)                                         
        xsinlat=sinlat*cosdel+coslat*sindel                               
        xcoslat=coslat*cosdel-sinlat*sindel                               
        ratio=xcoslat/coslat                                              
        geoid(1)=ratio*ru(1)                                      
        geoid(2)=ratio*ru(2)                                                    
        geoid(3)=xsinlat
        return
        end





        subroutine find_sc_axes_from_geoid(geoid,scvel, x,y,z)
        implicit none
        real(8), intent(in):: geoid(3), scvel(3)
        real(8), intent(out)::  x(3), y(3), z(3)
        z=-geoid
        call cross_norm(z,scvel, y)
        call cross_norm(y,z, x)
        return
        end subroutine find_sc_axes_from_geoid




        subroutine apply_rpy_tot_sc_axes(rpy, x,y,z)      

        implicit none

        !arguments
        real(8), intent(in):: rpy(3)
        real(8), intent(inout)::  x(3), y(3), z(3)
        
        !local variables
        real(8) roll,pitch,yaw
        real(8) att_mat(3,3),a0(3,3), a(3,3)
        real(8) cos_roll,sin_roll,cos_pitch,sin_pitch,cos_yaw,sin_yaw

        real(8) dsind, dcosd

        roll= rpy(1)
        pitch=rpy(2)
        yaw=  rpy(3)

        a0(1,:)=x
        a0(2,:)=y
        a0(3,:)=z

        cos_roll= dcosd(roll)
        sin_roll= dsind(roll)  
        cos_pitch=dcosd(pitch)
        sin_pitch=dsind(pitch)
        cos_yaw=  dcosd(yaw)
        sin_yaw=  dsind(yaw)

!     following is the 'rpy' or  '123' convention

        att_mat(1,1)= cos_pitch*cos_yaw       
        att_mat(1,2)= sin_roll*sin_pitch*cos_yaw + cos_roll*sin_yaw           
        att_mat(1,3)=-cos_roll*sin_pitch*cos_yaw + sin_roll*sin_yaw           
        att_mat(2,1)=-cos_pitch*sin_yaw       
        att_mat(2,2)=-sin_roll*sin_pitch*sin_yaw + cos_roll*cos_yaw           
        att_mat(2,3)= cos_roll*sin_pitch*sin_yaw + sin_roll*cos_yaw           
        att_mat(3,1)= sin_pitch       
        att_mat(3,2)=-sin_roll*cos_pitch              
        att_mat(3,3)= cos_roll*cos_pitch

        a=matmul(att_mat,a0)

        x=a(1,:)
        y=a(2,:)
        z=a(3,:)

        return
        end





      subroutine cross_norm(x, y, z) 
        real(8) x(3),y(3),z(3),xmag
         
      z(1)=x(2)*y(3)-x(3)*y(2)                                                    
      z(2)=x(3)*y(1)-x(1)*y(3)                                                    
      z(3)=x(1)*y(2)-x(2)*y(1)
                                                            
      xmag=sqrt(dot_product(z,z))
        z=z/xmag                                         
      return                                                            
      end subroutine cross_norm                                                              
   
      subroutine fd_sunr(scpos,scalt,sun, iflag_sun,sunr)
      implicit none
 
      integer(4), parameter :: nstep=9000
      real(8),    parameter ::  step=90./nstep
 
      real(8) scpos(3),scalt,sun(3),geoid(3),rearth_geoid,sunr(3)
      real(8) rs,rg,gs,aa,bb,fov_ang
      real(8) cos_alpha,alpha,beta_glint,ratio,z,y,y1,y2,e
      real(8) sinz(0:nstep),cosz(0:nstep),sin2z(0:nstep),sina,cosa
      integer(4) iflag_sun,istart,i,i_sv,n
 
      data istart/1/
 
      real(8) dsind, dcosd, dasind, dacosd

      if(istart.eq.1) then
         istart=0
         do i=0,nstep
            z=step*i
            sinz(i) =dsind(z)
            cosz(i) =dcosd(z)
            sin2z(i)=dsind(2*z)
         enddo
      endif

      call find_geoid(scpos, geoid,rearth_geoid) !geoid point away from earth like surface normal vector
 
      ratio=(scalt+rearth_geoid)/rearth_geoid
 
      fov_ang=dasind(1./ratio)
      cos_alpha=dot_product(sun,geoid)
      alpha=dacosd(cos_alpha)
 
      if(alpha.gt.180-fov_ang-0.1) then !0.1 is a tolerance to handle uncertainties at limb
         iflag_sun=0
         return
      endif
 
      iflag_sun=1
 
      if(alpha.le.1.0) then
         e=ratio-1
         beta_glint=alpha*(1+e)/(1+2*e)
         goto 100
      endif
 
      sina=dsind(0.5*alpha)
      cosa=dcosd(0.5*alpha)
        
      if(alpha.le.90) then
         n=int((0.5*alpha)/step) + 2 !n is the max limit for z when alpha<90
      else
         n=int((90-0.5*alpha)/step) + 2 !when alpha>90, max limit set by beta_glint must be 90 or less
      endif
 
      i_sv=0
      do i=1,n
         y=(sinz(i)*cosa + cosz(i)*sina)/sin2z(i)
         if(y.lt.ratio) then
            i_sv=i
            y2=y
            exit
         endif
         y1=y
      enddo

      if(i_sv.le.2) stop 'error1 in fd_sun'
 
      z=(i_sv-1 + (ratio-y1)/(y2-y1))*step
 
      beta_glint=z + 0.5*alpha
        
 100  continue
        
      rs=-dcosd(2*beta_glint)
      rg=-dcosd(2*beta_glint-alpha)
      gs= cos_alpha
        
      aa=(1-rg**2)/(rs-gs*rg)
      bb=rg-aa*gs
        
      sunr=aa*sun + bb*geoid          
        
      return
      end


      subroutine fd_boresight_pols(time_ut1_2000, cellat, cellon, bore_sight,
     1                             ve,he)
      implicit none

c     ve is unit vector in vpol direction
c     he is unit vector in hpol direction

      real(8),    intent(in)  ::  time_ut1_2000
      real(4),    intent(in)  ::  cellat,cellon
      real(8),    intent(in)  ::  bore_sight(3)
      real(8),    intent(out) ::  ve(3),he(3)

      real(8) rotangle,days,eci_lon, geoid(3),k(3)

      real(4) sind, cosd
      real(8) dsind, dcosd

      days=time_ut1_2000/86400.d0 -0.5d0
      call get_gm_angle(time_ut1_2000,days, rotangle)

      eci_lon=cellon + rotangle !convert from earth to inertial

      geoid(1)=cosd(cellat)*dcosd(eci_lon) 
      geoid(2)=cosd(cellat)*dsind(eci_lon) 
      geoid(3)=sind(cellat)

      k=-bore_sight !k is propagation vector

      call cross_norm(k,geoid, he)
      call cross_norm(he,k, ve)

      return
      end
  
 
