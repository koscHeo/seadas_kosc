        subroutine moon_vector(days, moonvec,moondis)    !returns vector in j2000 system, dis is km
        implicit none

        real(8), intent(in) :: days
        real(8), intent(out):: moonvec(3),moondis  
        real(8) t,l_0, l, l_p, f, d, l_m, b_m, dl, twiced, twicel, l_2d, l_p_2d,l_plus_2d, twicel_twiced
        real(8)  cosl_m, sinl_m, cosb_m, sinb_m, eps,coseps, sineps 

        real(8) dcosd, dsind

      t=days/36525.d0

        eps=23.4393d0 - 0.01300d0*t
        coseps=dcosd(eps)
        sineps=dsind(eps)

        l_0 = 218.31617d0 + 481266.48368d0*t     !degrees
        l   = 134.96292d0 + 477198.86753d0*t
        l_p = 357.52543d0 +  35999.04944d0*t
        f   =  93.27283d0 + 483202.01873d0*t
        d   = 297.85027d0 + 445267.11135d0*t

        twiced = 2*d
        twicel = 2*l
        l_2d = l - twiced
        l_p_2d = l_p - twiced
        l_plus_2d = l + twiced
        twicel_twiced = twicel - twiced
        
        dl = 22640*dsind(l)        + 769*dsind(twicel)   - 4586*dsind(l_2d)     + 2370*dsind(twiced)
     &      -668*dsind(l_p)      - 412*dsind(2*f)       - 212*dsind(twicel_twiced) 
     &      -206*dsind(l+l_p_2d) + 192*dsind(l_plus_2d) - 165*dsind(l_p_2d) 
     &      +148*dsind(l-l_p)    - 125*dsind(d)         - 110*dsind(l+l_p) 
     &       -55*dsind(2*f-twiced)   !arcseconds
                
        l_m = l_0 +  dl/3600.d0

        b_m = 18520*dsind(f + (dl +  412*dsind(2*f) +     541*dsind(l_p))/3600.d0)
     &       -526*dsind(f-twiced) + 44*dsind(l_2d+f )  - 31*dsind(f-l_plus_2d)
     &        -25*dsind(f-twicel) - 23*dsind(l_p_2d+f) + 21*dsind(f-l) 
     &       + 11*dsind(f-l_p-twiced)  !arcseconds

        b_m = b_m/3600.d0        !degrees

        moondis= 385000 - 20905*dcosd(l)         - 3699*dcosd(l_2d)         - 2956*dcosd(twiced)
     &                   -570*dcosd(twicel)    +  246*dcosd(twicel_twiced) - 205*dcosd(l_p_2d) 
     &                   -171*dcosd(l_plus_2d) -  152*dcosd(l + l_p_2d)          !km
        
        cosl_m = dcosd(l_m)
        sinl_m = dsind(l_m)
        cosb_m = dcosd(b_m)
        sinb_m = dsind(b_m)
                
        moonvec(1) = cosl_m * cosb_m
        moonvec(2) = coseps * sinl_m * cosb_m - sineps * sinb_m
        moonvec(3) = sineps * sinl_m * cosb_m + coseps * sinb_m     
        
        end subroutine moon_vector              
