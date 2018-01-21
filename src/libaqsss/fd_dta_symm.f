c     TM September 16, 2013  
c     correction to symmetrize ascending/descending swaths
c     can be regarded as an effective galactic correction
c     
c     it updates the empirical galaxy correction from March 2013 
 
c     the table is produced by the IDL routine dta_table_analysis.pro
c
c     time is observation time (seconds since jan 1 0z 2000)
c     zzang is intra-orbit angle in degrees

        subroutine fd_dta_symm(irad,time,zzang, time_galaxy, dtagal_symm_tab, dtagal_symm)
        implicit none

        integer(4), parameter :: n_rad=3
        integer(4), parameter :: nomega=1441, nzang=1441
        real(8),    parameter :: sidereal_year= 365.25636d0 ! astronomical almanac (days)

        integer(4), intent(in)        ::  irad
        real(8),    intent(in)        ::  time,zzang
        real(8),    intent(in)        ::  time_galaxy(nzang)
        real(8),    intent(in)        ::  dtagal_symm_tab(2,n_rad,nomega,nzang)
        real(4),    intent(out)       ::  dtagal_symm(3)

        real(8)                       ::  deltime,days
        real(4)                       ::  a1,a2,b1,b2,brief
        real(4)                       ::  dtagal_symm1(2,2)

        integer(4)                    ::  izzang1,izzang,i,j1,j2
        
        

        if(zzang.lt.0 .or. zzang.gt.360) stop 'zzang oob in fd_ta_symm, pgm stopped'

      brief=4*zzang
        if(brief.gt.1439.99) brief=1439.99
        izzang1=1+brief
        a1=izzang1-brief
        a2=1-a1
        if(izzang1.lt.1 .or. izzang1.gt.1440) stop 'i oob in fd_dta_symm, pgm stopped'

      do i=1,2
        izzang=izzang1+i-1

        deltime=(time-time_galaxy(izzang))/86400.d0  !convert to days
        days=deltime - sidereal_year*floor(deltime/sidereal_year)
        if(days.lt.0 .or. days.gt.sidereal_year)  stop 'days oob in fd_dta_galaxy, pgm stopped' ! i may have to cut some slack here

        brief=(nomega-1)*days/sidereal_year
        if(brief.gt.1439.99) brief=1439.99
        j1=1+brief
        j2=j1+1
        b1=j1-brief
        b2=1-b1
        if(j1.lt.1 .or. j2.gt.nomega) stop 'j oob in fd_ta_galaxy, pgm stopped'

        dtagal_symm1(1:2,i)=b1*dtagal_symm_tab(1:2,irad,j1,izzang) + b2*dtagal_symm_tab(1:2,irad,j2,izzang) ! omega interpolation

        enddo !i

        dtagal_symm(1:2) =(a1*dtagal_symm1(1:2,1) + a2*dtagal_symm1(1:2,2)) ! zzang interpolation
        dtagal_symm(3)   = 0.0
        

      return
      end subroutine fd_dta_symm
