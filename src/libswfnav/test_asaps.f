        real*8 tvec(1000),xvec(6,1000)
        real*8 orb(6),cdrg/2.02/,sec,xmin
        data li/30/,mi/15/,iplt/0/,iyr/1992/,nstp/720/

        print *,'Enter start day and minutes of day'
        accept *,iday,xmin
        sec = xmin*60.d0
        print *,'Enter orbital elements'
        accept *,orb
      call asaps(li,mi,iplt,orb,iyr,iday,sec,nstp,cdrg,tvec,xvec)
        print *,nstp,orb
        do i=1,nstp
          write(9,1100) tvec(i),(xvec(j,i),j=1,6)
        end do
 1100   format (f14.3,3f11.4,3f11.7)
        stop 
        end
