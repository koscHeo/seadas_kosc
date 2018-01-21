        subroutine geolocation( rpy_adj, time_utc_2000,scpos_j2k,scvel_j2k, scrpy,
     1                          cellat, cellon, celtht, celphi, suntht, sunphi, sunglt, moonglt,
     2                          glxlat, glxlon, zang, sun_zenith, sclon, sclat, scalt, cellat_corners, cellon_corners,
     3                          sund, sunr, moond, moonr, bore_sight)

        implicit none

        integer(4), parameter :: n_rad=3 ! number of radiometers (i.e., horns)  (inner, middle, outer)  

!     ==================================================================================================
!     ==================================== geolocation constants =======================================
!     ==================================================================================================
        real(8), parameter ::  global_ang=33.d0 !rotation angle of global system relative to s/c system
        real(8), parameter ::   beta0(n_rad)=(/ 25.828d0,  33.818d0, 40.367d0/) !boresight nadir   angle, Adams email 5/28/2009
        real(8), parameter ::  alpha0(n_rad)=(/  9.848d0, -15.291d0,  6.547d0/) !boresight azimuth angle, Adams email 5/28/2009
        real(8), parameter ::  half_beamwidth(n_rad)=(/3.07d0, 3.17d0, 3.24d0/) !half the full 3db beam widthd deg

        real(8), parameter :: re=6378.137d3  !earth equatorial radius (meters)
        real(8), parameter :: rp=6356.752d3 !earth      polar radius (meters)
        real(8), parameter :: ffac=(rp/re)*(rp/re)


        real(8) time_utc_2000, time_ut1_2000
        real(8),    dimension(3) :: scpos_j2k, scvel_j2k, scrpy

!     ==========================================================================================================
!     =================================== geolocation arrays ===================================================
!     ==========================================================================================================

!     iflag_sun is 1 if sun is visible from spacecraft (usual case), otherwise it = 0
!     iflag_moon is 1 if moon is visible from spacecraft (usual case), otherwise it = 0
!     zang is the intra-orbit angle. it is zero at the south  pole, 90 deg at the equator ascending node, etc.
!     sclat is the geodetic latitude of the spacecraft earth nadir point (deg)
!     sclon is the east longitude    of the spacecraft earth nadir point (deg)
!     scalt is the spacecraft altitude, nadir point to spacecraft  (meters)
!     sund  is the earth-to-sun unit vector in eci coordinates 
!     sunr  is the unit vector from spacecraft to sun reflection point on earth in eci coordinates
!     moond is the earth-to-moon unit vector in eci coordinates 
!     moonr is the unit vector from spacecraft to moon reflection point on earth
!     bore_sight is the unit vector from the spacecraft to the center of the observation cell

!     cellat is geodetic latitude (deg)
!       cellon is    east longitude (deg)
!     celtht is boresight earth incidence angle (deg)
!     celphi is boresight  earth-projection relative to clockwise from north (deg), boresight points towards earth
!     suntht is sun vector earth incidence angle (deg)
!     sunphi is sun vector earth_projection relative to clockwise from north (deg), sun vector points away from earth
!     sunglt   is the sun  glint angle: angle between the specular reflected boresight vector and the  sun vector (deg)
!     mooonglt is the moon glint angle: angle between the specular reflected boresight vector and the moon vector (deg)
!     glxlat is j2k declination latitude from where the specular galactic radiation originated (deg)
!     glxlon is j2k ascension from where the specular galactic radiation originated (deg)
!     cellat_corners is the geodetic latitudes of the four corners of the 3-db footprint
!     cellon_corners is the    east longitudes of the four corners of the 3-db footprint

        integer(4) iflag_sun, iflag_moon
        real(8)    zang,sclat,sclon,scalt,sun_zenith
        real(8),    dimension(      3) ::  sund,sunr,moond,moonr,rpy_adj
        real(8),    dimension(3,n_rad) ::  bore_sight
        real(4),    dimension(  n_rad) ::  cellat,cellon,celtht,celphi,suntht,sunphi,sunglt,moonglt,glxlat,glxlon
        real(4),    dimension(4,n_rad) ::  cellat_corners,cellon_corners

        real(8) secyr,secdy
        real(8) scpos(3),scvel(3)                                                                                                                        
        real(8) sun_j2k(3), sundis_km, moon_j2k(3), moondis_km
        real(8) ru(3),vu(3)
        real(8) rotangle,days,np(3,3),np_inv(3,3),det 
        real(8) r,v
        real(8) x(3), y(3), z(3),xg(3),yg(3),zg(3)
        real(8) b0(3),ref(3)
        real(8) xlat,xlon,range,azim,thtinc,sunang,moonang,thtinc_sun,azim_sun
        real(8) dummy
        real(8) cos_global,sin_global
        real(8) d_alpha(n_rad),alphax,betax,b0x(3)

        integer(4) irad,icorner,istart
        integer(4) lyear,idayjl,imon,idaymo,isecdy

        data dummy/0./
        data istart/1/

        real(8) dsind, dcosd, dasind, dacosd, datan2d

!        real(8), dimension(n_rad), parameter  :: rpy_adj = (/-0.51, + 0.16, 0.00/)

        if(istart.eq.1) then
           istart=0
           cos_global=dcosd(global_ang)
           sin_global=dsind(global_ang)
           do irad=1,3
              d_alpha(irad)=dacosd( ( dcosd(half_beamwidth(irad))-(dcosd(beta0(irad)))**2 ) / (dsind(beta0(irad)))**2 )
           enddo
        endif

c     ======================================================================================================================
c     ============================================= do time computations ===================================================
c     ======================================================================================================================
        call fd_date_2000(time_utc_2000,   secyr,lyear,idayjl,imon,idaymo,secdy)
        isecdy=int(secdy)       !truncate to make sure it does not = 86400 for the wind and sst programs
        time_ut1_2000 = time_utc_2000
        days=time_ut1_2000/86400.d0 -0.5d0
        call get_gm_angle(time_ut1_2000,days, rotangle)
        call fd_precession(days, np) !i do not bother doing the nutation.  its effect is way to small for aquarius
        call invert_3by3(np, np_inv,det)
        if(abs(det-1).gt.1.e-6) stop 'error in np_inv,pgm stopped'

        scpos=matmul(np,scpos_j2k(:)) !convert to eci 
        scvel=matmul(np,scvel_j2k(:)) !convert to eci 

c     ======================================================================================================================
c     ============================= do coorindate transformations and s/c location==========================================
c     ======================================================================================================================

        call  sun_vector(days,  sun_j2k, sundis_km) !returns unit earth-to-sun vector in j2K system, dis is km
        sund(:)=matmul(np,  sun_j2k) !convert to eci

        call moon_vector(days, moon_j2k,moondis_km) !returns unit earth-to-moon vector in j2k, dis is km
        moond(:)=matmul(np,moon_j2k) !convert to eci

        r=sqrt(dot_product(scpos,scpos))
        v=sqrt(dot_product(scvel,scvel))
 
        ru=scpos/r
        vu=scvel/v
 
        sun_zenith=dacosd(dot_product(ru,sund(:)))

        zang=datan2d(ru(3),vu(3)) + 90.
        if(zang.ge.360.) zang=zang-360.
        if(zang.lt.  0.) zang=zang+360.

c     find s/c x,y,z axes and global xg,yg,zg axes

        z=-ru 
        call cross_norm(z,scvel, y)
        call cross_norm(y,z, x)
        call apply_rpy_tot_sc_axes(scrpy(:), x,y,z)
         
c       pointing adjustment (Kyle)
        call apply_rpy_tot_sc_axes(rpy_adj, x,y,z)     
   
        xg=y*cos_global - z*sin_global
        yg=-x
        zg=y*sin_global + z*cos_global

        call fd_cell_parameters(1,rotangle,r,ru,x,y,z,dummy,dummy,sund(:),moond(:),
     &                          b0,ref,sclat,sclon,scalt,thtinc,azim,thtinc_sun,azim_sun,sunang,moonang)         !s/c lat, lon, and alt
      
c     ======================================================================================================================
c     ========= compute some general parameters that are independent of horn, including the sun glint vector ===============
c     ======================================================================================================================

c     sunr(3,icyc) is from unit vector from spacecraft to sun reflection point on earth
c     iflag_sun(icyc)=1 if sun is visible from spacecraft, otherwise it = 0
        call fd_sunr(scpos,scalt,sund(:), iflag_sun,sunr(:)) 

c     moonr(3,icyc) is from unit vector from spacecraft to moon reflection point on earth
c     iflag_moon(icyc)=1 if moon is visible from spacecraft, otherwise it = 0
        call fd_sunr(scpos,scalt,moond(:), iflag_moon,moonr(:)) 
c     ======================================================================================================================
c     ================================== compute earth ta and sun ta for each horn =========================================
c     ======================================================================================================================

        do irad=1,n_rad                                                                 

           call fd_cell_parameters(2,rotangle,r,ru,x,y,z,beta0(irad),alpha0(irad),sund(:),moond(:),
     &                          b0,ref,xlat,xlon,range,thtinc,azim,thtinc_sun,azim_sun,sunang,moonang) !find boresight and cell location

           bore_sight(:,irad)=b0

           cellat(irad)=real(xlat)
           cellon(irad)=real(xlon)
           celtht(irad)=real(thtinc)                  
           celphi(irad)=real(azim)
           suntht(irad)=real(thtinc_sun)
           sunphi(irad)=real(azim_sun)
           sunglt(irad)=real(sunang)
           moonglt(irad)=real(moonang)

           ref=matmul(np_inv,ref) !convert reflection vector back to j2000
           glxlat(irad)=real(dasind(ref(3))) !geocentric latitude
           glxlon(irad)=real(datan2d(ref(2),ref(1)))
           if(glxlon(irad).lt.  0) glxlon(irad)=glxlon(irad)+360
           if(glxlon(irad).ge.360) glxlon(irad)=glxlon(irad)-360

           do icorner=1,4       !get footprint location
        
              if(icorner.eq.1) then
                 betax = beta0(irad) + half_beamwidth(irad)
                 alphax=alpha0(irad)
              endif
        
              if(icorner.eq.2) then
                 betax = beta0(irad) 
                 alphax=alpha0(irad) + d_alpha(irad)
              endif
         
              if(icorner.eq.3) then
                 betax = beta0(irad) - half_beamwidth(irad)
                 alphax=alpha0(irad) 
              endif

              if(icorner.eq.4) then
                 betax = beta0(irad) 
                 alphax=alpha0(irad) - d_alpha(irad)
              endif

              call fd_cell_parameters(2,rotangle,r,ru,x,y,z,betax,alphax,sund(:),moond(:),
     &                          b0x,ref,xlat,xlon,range,thtinc,azim,thtinc_sun,azim_sun,sunang,moonang) !find boresight and cell location

              cellat_corners(icorner,irad)=real(xlat)
              cellon_corners(icorner,irad)=real(xlon)
           enddo                !icorner

        enddo                   !irad

        return
        end


