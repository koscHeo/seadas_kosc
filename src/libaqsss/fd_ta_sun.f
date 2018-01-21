c     time is observation time (seconds since jan 1 0z 2000)
c     zang is intra-orbit angle in degrees

        subroutine fd_ta_sun(time,zang, time_sun, tasun_dir_tab,tasun_ref_tab, tasun_dir,tasun_ref)
        implicit none

        integer(4), parameter :: n_rad=3 ! number of radiometers (i.e., horns)  (inner, middle, outer)
        integer(4), parameter :: nomega=1441, nzang=1441
        real(8), parameter :: tropical_year= 365.24219d0 ! wikipedia (days)

        real(4),    parameter :: transq(n_rad)=(/ 0.98294, 0.98105, 0.97864/) !avg values for 2007 run, see test_retreival\memo10
!        real(4),    parameter :: refl_cor(2*n_rad)=(/0.98861,0.97190,0.98863,0.98177,0.98910,0.98571/) !refl 7.5 m/s correction, see '\taspace\memo3.txt'
        real(4),    parameter :: refl_cor(2*n_rad)=(/0.99023,0.97198,0.98926,0.96685,0.98740,0.96534/) ! V2.0  08/10/12

        real(8),    intent(in)  :: time,zang
        real(4), intent(in) ::  tasun_dir_tab(nomega,nzang,3,n_rad),tasun_ref_tab(nomega,nzang,3,n_rad)
        real(8), intent(in) ::  time_sun(nzang)
        real(4),    intent(out) :: tasun_dir(3,n_rad),tasun_ref(3,n_rad)

        real(8) deltime,days
        real(4) a1,a2,b1,b2,brief
        real(4) tasun_dirx(3,n_rad,2),tasun_refx(3,n_rad,2)

        integer(4) irad,izang1,izang,i,j1,j2

        if(zang.lt.0 .or. zang.gt.360) stop 'zang oob in fd_ta_sun, pgm stopped'

        brief=4*zang
        if(brief.gt.1439.99) brief=1439.99
        izang1=1+brief
        a1=izang1-brief
        a2=1-a1
        if(izang1.lt.1 .or. izang1.gt.1440) stop 'i oob in fd_ta_sun, pgm stopped'

        do i=1,2
           izang=izang1+i-1

           deltime=(time-time_sun(izang))/86400.d0 !convert to days
           days=deltime - tropical_year*floor(deltime/tropical_year)
           if(days.lt.0 .or. days.gt.tropical_year)  stop 'days oob in fd_ta_sun, pgm stopped' ! i may have to cut some slack here

           brief=(nomega-1)*days/tropical_year
           if(brief.gt.1439.99) brief=1439.99
           j1=1+brief
           j2=j1+1
           b1=j1-brief
           b2=1-b1
           if(j1.lt.1 .or. j2.gt.nomega) stop 'j oob in fd_ta_sun, pgm stopped'
           if(j2.eq.nomega) j2=1 !use wrap around rather than element 1441, but the dif is extremely small see \ta_space\memo3

           do irad=1,n_rad
              tasun_dirx(:,irad,i)=b1*tasun_dir_tab(j1,izang,:,irad) + b2*tasun_dir_tab(j2,izang,:,irad)
              tasun_refx(:,irad,i)=b1*tasun_ref_tab(j1,izang,:,irad) + b2*tasun_ref_tab(j2,izang,:,irad)
           enddo

        enddo                   !i

        do irad=1,n_rad
           tasun_dir(:,irad)=a1*tasun_dirx(:,irad,1) + a2*tasun_dirx(:,irad,2)
           tasun_ref(:,irad)=a1*tasun_refx(:,irad,1) + a2*tasun_refx(:,irad,2)
           tasun_ref(1,irad)=refl_cor(1+(irad-1)*2)*tasun_ref(1,irad)
           tasun_ref(2,irad)=refl_cor(2+(irad-1)*2)*tasun_ref(2,irad)
           tasun_ref( : ,irad)=transq(irad)*tasun_ref(:,irad)
        enddo

        return
        end     
