      subroutine fd_ta_moon(irad, bore_sight, moonr,refl,transq,  tamon_ref)
        implicit none

        integer(4), intent(in)  :: irad
        real(8),    intent(in)  :: bore_sight(3),moonr(3)
        real(4),    intent(in)  :: refl(2),transq
        real(4),    intent(out) :: tamon_ref(3)

      real(4) moon_ang,tbv,tbh,tb_moon_toa(2),gain_fac
        real(4) gain(2,2,3), half_beamwidth(3), tb_mon_solid_angle(3)
        real(8) dp, das

      data gain/  74.41968,  -0.47552,  -0.42802,  74.36780,  70.84354,   0.13980,   0.12591,  70.82498,
     &            65.86123,  -0.58609,  -0.52414,  65.78670/

      data half_beamwidth/3.04e0, 3.17e0, 3.24e0/

        data tb_mon_solid_angle/ 0.010813, 0.010430,  0.009982/ ! 'O:\aquarius\moon\memo1.txt'

        real(8) dasind

!        moon_ang=2*dasind(sqrt(dot_product(moonr-bore_sight,moonr-bore_sight))/2.)
        ! Check for asin() argument greater than 1  JMG  02/10/12
        dp = sqrt(dot_product(moonr-bore_sight,moonr-bore_sight)) / 2.
        if ( dp .gt. 1) dp = 1.0
        moon_ang = 2*dasind(dp)

        tbv=transq*refl(1)*tb_mon_solid_angle(irad)
        tbh=transq*refl(2)*tb_mon_solid_angle(irad)

        tb_moon_toa(1)=tbv+tbh
        tb_moon_toa(2)=tbv-tbh

        if(moon_ang.gt.15) then
        tamon_ref=0
        return
        endif

        gain_fac=10.**(-0.3*(moon_ang/half_beamwidth(irad))**2)

        tamon_ref(1)=(gain(1,1,irad)*tb_moon_toa(1) + gain(2,1,irad)*tb_moon_toa(2))*gain_fac
        tamon_ref(2)=(gain(1,2,irad)*tb_moon_toa(1) + gain(2,2,irad)*tb_moon_toa(2))*gain_fac
        tamon_ref(3)=0
        return
        end


