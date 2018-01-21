        subroutine fd_land(irad, zang,sclon, fpt_lnd, frc_lnd, fland,gland)
        implicit none

        integer(4), parameter :: n_rad=3  
        integer(4), parameter :: nlon_lnd=2160, nzang_lnd=2881

        integer(4), intent(in)  :: irad
        real(8),    intent(in)  :: zang,sclon
        real(4),    intent(in)  :: fpt_lnd(nlon_lnd,nzang_lnd,n_rad)
        real(4),    intent(in)  :: frc_lnd(nlon_lnd,nzang_lnd,n_rad)
        real(4),    intent(out) :: fland,gland

        integer(4) i1,i2,j1,j2
        real(4) a1,a2,b1,b2,brief


        if(sclon.lt.0 .or. sclon.gt.360) stop 'sclon oob in fd_land, pgm stopped'
        if(zang .lt.0 .or.  zang.gt.360) stop 'zang  oob in fd_land, pgm stopped'

      brief=8*zang
        if(brief.gt.2879.99) brief=2879.99
        i1=1+brief
        i2=i1+1
        a1=i1-brief
        a2=1-a1
        if(i1.lt.1 .or. i1.gt.nzang_lnd-1) stop 'i oob in fd_land, pgm stopped'

        brief=6.*sclon
        if(brief.gt.2159.99) brief=2159.99
        j1=1+brief
        j2=j1+1
        b1=j1-brief
        b2=1-b1
        if(j1.lt.1 .or. j1.gt.nlon_lnd) stop 'j oob in fd_land, pgm stopped'
        if(j2.eq.nlon_lnd+1) j2=1  !wrap around

        fland=a1*(b1*fpt_lnd(j1,i1,irad) + b2*fpt_lnd(j2,i1,irad)) + a2*(b1*fpt_lnd(j1,i2,irad) + b2*fpt_lnd(j2,i2,irad)) 
        gland=a1*(b1*frc_lnd(j1,i1,irad) + b2*frc_lnd(j2,i1,irad)) + a2*(b1*frc_lnd(j1,i2,irad) + b2*frc_lnd(j2,i2,irad)) 

      return
        end     
