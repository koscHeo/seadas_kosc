subroutine tb2vsm(theta,tbright,tsoil,sand,clay,por,vc,b_factor,tau,omega,h,mv)

        ! This subroutine performs soil moisture retrieval based on input
        ! brightness temperature, air/surface temperature, incidence angle
        ! and soil texture information and outputs volumetric soil moisture.
        !
        ! Routines called
        ! ---------------
        ! wang_sch - converts dielectric constant to volumetric soil moisture
        !            using the empirical equations developed by Wang and
        !            Schmugge
        !
        ! Variables used in this subroutine:
        !
        ! Input
        ! tex_index - soil texture index ranging from 1 to 7
        ! ndvi      - average ndvi values based on 2001 to 2011 MODIS data
        ! tbright   - Aquarius brightness temperature at 1.4 GHz H polarization
        ! tsoil      - surface temperature provided by NCEP
        !
        ! Output
        ! mv        - volumetric soil moisture

        !  byte sand,clay,por
        real theta,tbright,tsoil,sand,clay,por,ndvi,vc,b_factor,tau,omega,h,mv

        pi=acos(-1.0)
        ctheta = cos(theta*pi/180.)
        stheta = sin(theta*pi/180.)

        e_target=tbright/tsoil

        tau=b_factor*vc

        opt2=exp(-2.0*b_factor*vc/ctheta)
        opt=exp(-1.0*b_factor*vc/ctheta)
        ! opt=exp(-1.0*tau/ctheta)
        ! e_sur=1.0 - (1.0-e_target)/opt2
        e_sur=(e_target-1.0+opt2+omega-omega*opt2)/(opt2+omega*opt-omega*opt2)

        if ((e_sur.le.0.0).or.(e_sur.ge.1.0)) then !overcorrection on vegetation April,2007
                mv=por/100.0
                mv=-999.0
                return
        endif

        ! correct for surface roughness effect - Choudhury et al., 'Eeffect of
        ! surface roughness on microwave emission of soils', J. Geophys. Res., 84,
        ! pp. 5699-5706.

        e_soil=1.0 + (e_sur-1.0)*exp(h*ctheta**2)
        rs= 1.0 - e_soil
        if(rs .ge. 0) then
                x=sqrt(rs)
        a=(ctheta**2)*(1.0-x)**2
        b=-(1.0+x)**2
        c=(stheta**2)*(1.0+x)**2
        if ((b**2 - 4*a*c).ge.0.0) then
                eps_r = (-1.0*b + sqrt(b**2 - 4*a*c))/(2.0*a)
!       if (eps_r.lt.0.0) eps_r = (-1.0*b + sqrt(b**2 - 4*a*c))/(2.0*a)
                eps_r=stheta**2+(ctheta*(-1.0-x)/(x-1.0))**2
                call wang_sch(sand,clay,por,tsoil,eps_r,mv)
        else
        mv=-999.0
        endif
        else
                mv=-999.0
        endif
        return

end subroutine tb2vsm

subroutine wang_sch(sand,clay,por,temp,esp_r,sm)

        ! The method used in this program to convert real part of dielectric
        ! constant to volumetric soil moisture is based on J. Wang and T. J.
        ! Schmugge, "An empirical model for the complex dielectric permittivity
        ! of soil as a function of water content," IEEE Transactions on Geo-
        ! science and Remote Sensing, GE 18:288-295, 1980.
        !
        ! Routines called
        ! ---------------
        ! None
        !
        ! Variables used in this subroutine:
        !
        ! Input
        ! sand  - percent of sand
        ! clay  - percent of clay
        ! por   - percent of porosity
        ! temp  - air/surface temperature in degree k
        ! esp_r - real part of dielectric constant
        !
        ! Output
        ! sm    - volumetric soil moiture
        !
        !  byte sand,clay,por
        real*8 freq, ew0,rt, esp_water_r,esp_i
        real p_sand,p_clay,porosity

        data esp_rock_r /5.5/   !* real part of dielectric constant for rock
        data esp_ice_r /3.2/    !* real part of dielectric constant for ice
        data freq/1.41e+09/     !* freq of Aquarius data

        p_sand=sand
        p_clay=clay
        porosity=por/100.0

!       esp_r=esp_r+1.0
!       esp_i=(esp_r-esp_ice_r)/4.0
!       esp_r=(esp_r**2.0+esp_i**2.0)**0.5

        wp = 0.06774 - 0.00064*p_sand + 0.00478*p_clay
        gamma = 0.481 - 0.57*wp
        wt = 0.165 + 0.49*wp
        tt = temp - 273.15

        ! This section computes water dielectric constant based on
        ! L. A. Klein and C. T. Swift, "An improved model for the
        ! dielectric constant of sea water at microwave frequencies,"
        ! IEEE Trans. Antenna and Prop., AP-25:104,111, 1977.

        ew0 = 88.045 - 0.4147*tt + 6.295e-04*tt**2 + &
                1.075e-05*tt**3
        rt = 1.1109e-10 - 3.824e-12*tt + 6.938e-14*tt**2 - &
                5.096e-16*tt**3
        rt = rt * freq
        esp_water_r=4.9 + (ew0-4.9)/(1+rt**2)
        a = (esp_water_r - esp_ice_r) * gamma / wt
        b = esp_ice_r - 1.0
        c = porosity + (1.0 - porosity)*esp_rock_r - esp_r
        if ((b**2 - 4*a*c).ge.0.0) then
                ws1 = (-1.0*b + sqrt(b**2 - 4*a*c))/(2.0*a)
                ex = esp_ice_r + (esp_water_r - esp_ice_r)*gamma
                ws2 = (esp_r - wt*ex + wt*esp_water_r - porosity - &
                        (1.0-porosity)*esp_rock_r)/(esp_water_r-1.0)

                if (ws1 .gt. wt) then
                        sm = ws2
                else
                        sm = ws1
                end if
        else
                sm = -10.0/a
!       print*,b,a
        endif
!       if (esp_r.le.2.8) then
!               sm=0.01
!       end if
        if (esp_r.ge.esp_water_r) then
                sm=1.0
        end if
        if (sm.le.0.02) sm=0.1+sm
        if (sm.le.0.02) sm=0.02
        return

end subroutine wang_sch
