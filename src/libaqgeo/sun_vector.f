C     This routine is based on:
C     1.  Supplement to Astron. Almanac. 1992 
C     2.  Satellite Orbits: Models, Methods, and Applications by Montenbruck and Gill.
C     See MEMO2

      SUBROUTINE sun_vector(days, sunVec,sunDis)         !returns vector in J2000 system, dis is km
        IMPLICIT NONE

        REAL(8), INTENT(IN):: days
        REAL(8), INTENT(OUT)::sunVec(3),sunDis
        REAL(8) T,XM,L_M,cosL_M,sinL_M,eps,coseps,sineps
        real(8) dsind, dcosd

      T=days/36525.D0

        eps=23.4393D0 - 0.01300D0*T
        coseps=dcosd(eps)
        sineps=dsind(eps)

        XM=357.528D0 + 35999.050D0*T
        L_M= 280.460d0 + 35999.3728D0*T + 1.915D0*dSIND(XM) + 0.020D0*dSIND(2*XM) 

        sunDis  =(149.619D0 - 2.499D0*dCOSD(XM) - 0.021D0*dCOSD(2*XM))*1.D6  !km

        cosL_M = dcosd(L_M)
        sinL_M = dsind(L_M)

        sunVec(1) = cosL_M
        sunVec(2) = coseps * sinL_M  
        sunVec(3) = sineps * sinL_M  
        RETURN
        END


