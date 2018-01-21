        real*8 vecs(6,721),driv(6,3,721),gpsec(721),orb(6),secin, t
        ngps = 721
        igyr = 1992


        print *, 'Enter element epoch'
        accept *, iyinit, idinit, secin
        print *, 'Enter orbital elements'
        accept *, orb

        open(11,status='old',readonly)
        call cdata

        do i=1,721
          read (11,1100) t, (vecs(j,i),j=1,6)
          if (i.eq.1) igday = t
          gpsec(i) = (t-igday)*864.d2
        end do
        call pderiv(ngps,igyr,igday,gpsec,vecs,orb,
     *    iyinit,idinit,secin,driv)
        do m=1,3
          do i=1,721
            write (12,1101) gpsec(i),(driv(j,m,i),j=1,6)
          end do
        end do

        stop
 1100   format(f14.8,3f11.4,3f11.7)
 1101   format(f14.3,f11.4,f11.2,2f11.4,f11.6,f11.4)
        end
