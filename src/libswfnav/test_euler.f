        real*4 ang(3), mpry(3,3)

        radeg=180.d0/acos(-1.d0)
        print *, 'Enter yaw, roll, pitch'
        accept *, ang

        call euler(ang, mpry)

        print *, mpry

        ang(1) = radeg*atan2( mpry(2,3), mpry(3,3) )
        ang(2) = radeg*asin( mpry(1,3) )
        ang(3) = radeg*atan2( -mpry(1,2), mpry(1,1) )
        
        print *, ang

        stop 
        end
