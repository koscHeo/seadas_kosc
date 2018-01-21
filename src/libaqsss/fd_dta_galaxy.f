c     empirical galaxy correction  (I)
c     time is observation time (seconds since jan 1 0z 2000)
c     zang is intra-orbit angle in degrees

        subroutine fd_dta_galaxy(irad,time,zang, time_galaxy,
     1     dtagal_ref_tab, dtagal_ref)
        implicit none

        integer(4), parameter :: n_rad=3
        integer(4), parameter :: nomega=1441, nzang=1441
        real(8), parameter :: sidereal_year= 365.25636d0 ! astronomical almanac (days)

        integer(4), intent(in)        ::  irad
        real(8),    intent(in)        ::  time,zang
        real(8),    intent(in)        ::  time_galaxy(nzang)
        real(8),    intent(in)        ::  dtagal_ref_tab(n_rad,nomega,nzang)
        real(4),    intent(out)       ::  dtagal_ref(3)

        real(8)                       ::  deltime,days
        real(4)                       ::  a1,a2,b1,b2,brief
        real(4)                       ::  dtagal_ref1(2)

        integer(4)                    ::  izang1,izang,i,j1,j2


        if(zang.lt.0 .or. zang.gt.360) stop 'zang oob in fd_ta_galaxy, pgm stopped'

        brief=4*zang
        if(brief.gt.1439.99) brief=1439.99
        izang1=1+brief
        a1=izang1-brief
        a2=1-a1
        if(izang1.lt.1 .or. izang1.gt.1440) stop 'i oob in fd_dta_galaxy, pgm stopped'

        do i=1,2
           izang=izang1+i-1

           deltime=(time-time_galaxy(izang))/86400.d0 !convert to days
           days=deltime - sidereal_year*floor(deltime/sidereal_year)
           if(days.lt.0 .or. days.gt.sidereal_year)  stop 'days oob in fd_dta_galaxy, pgm stopped' ! i may have to cut some slack here

           brief=(nomega-1)*days/sidereal_year
           if(brief.gt.1439.99) brief=1439.99
           j1=1+brief
           j2=j1+1
           b1=j1-brief
           b2=1-b1
           if(j1.lt.1 .or. j2.gt.nomega) stop 'j oob in fd_ta_galaxy, pgm stopped'

           dtagal_ref1(i)=b1*dtagal_ref_tab(irad,j1,izang) + b2*dtagal_ref_tab(irad,j2,izang) ! omega interpolation

        enddo                   !i

        dtagal_ref(1) =(a1*dtagal_ref1(1) + a2*dtagal_ref1(2)) ! zang interpolation
        dtagal_ref(2) = 0.d0
        dtagal_ref(3) = 0.d0
        

      return
      end subroutine fd_dta_galaxy    
