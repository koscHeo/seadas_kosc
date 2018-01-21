c     time is observation time (seconds since jan 1 0z 2000)
c     zang is intra-orbit angle in degrees

        subroutine fd_ta_galaxy(irad,time,zang,wind, eia_galaxy, tagal_dir_tab, tagal_ref_tab, time_galaxy,
     1                          eia_gal,tagal_dir,tagal_ref)
        implicit none

        integer(4), parameter :: nomega=1441, nzang=1441
        integer(4), parameter :: n_rad=3
        real(8), parameter :: sidereal_year= 365.25636d0 ! astronomical almanac (days)

        integer(4), intent(in)  :: irad
        real(8),    intent(in)  :: time,zang
        real(4),    intent(in)  :: wind
        real(4),    intent(in)  :: eia_galaxy(nzang,n_rad)
        real(4),    intent(in)  :: tagal_dir_tab(nomega,nzang,3,n_rad)
        real(4),    intent(in)  :: tagal_ref_tab(nomega,nzang,3,n_rad,5)
        real(8),    intent(in)  :: time_galaxy(nzang)
        real(4),    intent(out) :: eia_gal,tagal_dir(3),tagal_ref(3)




      real(8) deltime,days
        real(4) a1,a2,b1,b2,c1,c2,brief
      real(4) tagal_dir0(3,2),tagal_ref1(3,2),tagal_ref2(3,2)

        integer(4) izang1,izang,i,j1,j2,iwin1,iwin2

      if(wind.lt.0)                  stop 'wind oob in fd_ta_galaxy, pgm stopped'
        if(zang.lt.0 .or. zang.gt.360) stop 'zang oob in fd_ta_galaxy, pgm stopped'

      brief=wind/5.
        if(brief.gt.3.999) brief=3.999
        iwin1=1+brief
        iwin2=iwin1+1
        c1=iwin1-brief
        c2=1-c1


      brief=4*zang
        if(brief.gt.1439.99) brief=1439.99
        izang1=1+brief
        a1=izang1-brief
        a2=1-a1
        if(izang1.lt.1 .or. izang1.gt.1440) stop 'i oob in fd_ta_galaxy, pgm stopped'

        eia_gal=a1*eia_galaxy(izang1,irad) + a2*eia_galaxy(izang1+1,irad)

      do i=1,2
        izang=izang1+i-1

        deltime=(time-time_galaxy(izang))/86400.d0  !convert to days
        days=deltime - sidereal_year*floor(deltime/sidereal_year)
        if(days.lt.0 .or. days.gt.sidereal_year)  stop 'days oob in fd_ta_galaxy, pgm stopped' ! i may have to cut some slack here

        brief=(nomega-1)*days/sidereal_year
        if(brief.gt.1439.99) brief=1439.99
        j1=1+brief
        j2=j1+1
        b1=j1-brief
        b2=1-b1
        if(j1.lt.1 .or. j2.gt.nomega) stop 'j oob in fd_ta_galaxy, pgm stopped'

      tagal_dir0(:,i)=b1*tagal_dir_tab(j1,izang,:,irad)       + b2*tagal_dir_tab(j2,izang,:,irad)
      tagal_ref1(:,i)=b1*tagal_ref_tab(j1,izang,:,irad,iwin1) + b2*tagal_ref_tab(j2,izang,:,irad,iwin1)
      tagal_ref2(:,i)=b1*tagal_ref_tab(j1,izang,:,irad,iwin2) + b2*tagal_ref_tab(j2,izang,:,irad,iwin2)

        enddo !i

        tagal_dir=    a1*tagal_dir0(:,1) + a2*tagal_dir0(:,2)
        tagal_ref=c1*(a1*tagal_ref1(:,1) + a2*tagal_ref1(:,2)) +  c2*(a1*tagal_ref2(:,1) + a2*tagal_ref2(:,2))

      return
        end     
