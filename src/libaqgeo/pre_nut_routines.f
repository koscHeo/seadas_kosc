C     June 15, 2005, 01:30:17 PM version change on July 2, 2005.  Only comments were changed regarding utc vs ut1.
C     No functional changes were made.
C     
        SUBROUTINE FD_PRECESSION(DAYS, PRECESSION)
        IMPLICIT NONE
        REAL(8) DAYS,T,T2,T3,CHI,ZZZ,THT,COSC,COSZ,COST,SINC,SINZ,SINT,PRECESSION(3,3)
        real(8) dsind, dcosd
C     Table 3.211.1 with capital T equal 0

C       t=(UTCSec93/86400.d0 + Dif_EpOCH)/daysPerCentury 
        t=DAYS/36525.D0
        t2=t*t
        t3=t*t2
        CHI=(2306.2181D0*t + 0.30188D0*t2 + 0.017998D0*t3)/3600.D0      
        ZZZ=(2306.2181D0*t + 1.09468D0*t2 + 0.018203D0*t3)/3600.D0
        THT=(2004.3109D0*t - 0.42665D0*t2 - 0.041833D0*t3)/3600.D0

        COSC=dCOSD(CHI)
        COSZ=dCOSD(ZZZ)
        COST=dCOSD(THT)
        SINC=dSIND(CHI)
        SINZ=dSIND(ZZZ)
        SINT=dSIND(THT)

C       Eq.  3.21-8
        PRECESSION(1,1)= COSZ*COST*COSC - SINZ*SINC
        PRECESSION(1,2)=-COSZ*COST*SINC - SINZ*COSC
        PRECESSION(1,3)=-COSZ*SINT
        PRECESSION(2,1)= SINZ*COST*COSC + COSZ*SINC
        PRECESSION(2,2)=-SINZ*COST*SINC + COSZ*COSC
        PRECESSION(2,3)=-SINZ*SINT
        PRECESSION(3,1)= SINT*COSC
        PRECESSION(3,2)=-SINT*SINC
        PRECESSION(3,3)= COST


        RETURN
        END


C     Note that this routine does one additioal rotation to put the x axis at the mean equinox
        SUBROUTINE FD_NUTATION(DAYS, NUTATION)
        IMPLICIT NONE
        REAL(8) DAYS,T,T2,T3,E0,E,DPSI,DEPS,COSP,COSE,COS0,SINP,SINE,SIN0,NUTATION(3,3)
        REAL(8) R_ANGLE,COSR,SINR 
        REAL(8) N11,N12,N13,N21,N22,N23,N31,N32,N33
        real(8) dsind, dcosd

        t=DAYS/36525.D0
        t2=t*t
        t3=t*t2

c     Eq 3.222-1
        E0=(84381.448D0 - 46.8150D0*t -0.00059D0*t2 + 0.001813D0*t3)/3600.D0 !DEG

C     Eq. 3.225-4
        DPSI=-0.0048D0*dSIND(125.D0-0.05295D0*DAYS) - 0.0004D0*dSIND(200.9D0 + 1.97129D0*DAYS) !DEG
        DEPS= 0.0026D0*dCOSD(125.D0-0.05295D0*DAYS) + 0.0002D0*dCOSD(200.9D0 + 1.97129D0*DAYS)   !DEG

c      CALL iau_NUT00B(2451545.D0,DAYS, DPSI, DEPS)
c      DPSI=DPSI*180.D0/3.141592654D0
c      DEPS=DEPS*180.D0/3.141592654D0

        E=E0 + DEPS

        COS0=dCOSD(E0)
        COSE=dCOSD(E)
        COSP=dCOSD(DPSI)
        SIN0=dSIND(E0)
        SINE=dSIND(E)
        SINP=dSIND(DPSI)

C       Eq.  3.222-4
        N11= COSP 
        N12=-SINP*COS0
        N13=-SINP*SIN0
        N21= SINP*COSE
        N22= COSP*COSE*COS0 + SINE*SIN0
        N23= COSP*COSE*SIN0 - SINE*COS0
        N31= SINP*SINE
        N32= COSP*SINE*COS0 - COSE*SIN0
        N33= COSP*SINE*SIN0 + COSE*COS0

        R_ANGLE=DPSI*COSE

        COSR=dCOSD(R_ANGLE)
        SINR=dSIND(R_ANGLE)

        NUTATION(1,1)= COSR*N11 + SINR*N21
        NUTATION(1,2)= COSR*N12 + SINR*N22
        NUTATION(1,3)= COSR*N13 + SINR*N23

        NUTATION(2,1)=-SINR*N11 + COSR*N21
        NUTATION(2,2)=-SINR*N12 + COSR*N22
        NUTATION(2,3)=-SINR*N13 + COSR*N23
  
        NUTATION(3,1)= N31
        NUTATION(3,2)= N32
        NUTATION(3,3)= N33
        RETURN
        END


        SUBROUTINE get_gm_angle(UT1Sec93,DAYS, angle)

        IMPLICIT NONE

        REAL(8), INTENT(IN)::UT1Sec93,DAYS
        REAL(8), INTENT(OUT)::angle

        REAL(8), PARAMETER:: daysPerCentury = 36525.d0
        REAL(8), PARAMETER:: C0 = 24110.54841d0, C1 = 8640184.812866d0, C2 = 0.093104d0, C3 = -6.2d-6
        REAL(8)  T,T2,T3,UT1

        t=DAYS/36525.D0
        t2=t*t
        t3=t*t2

        UT1=UT1Sec93 - 86400.D0*FLOOR(UT1Sec93/86400.D0)
        angle = C0 + C1*t  + C2*t2 + C3*t3 +  UT1 !seconds
        angle = angle/240.d0       !convert seconds to degrees
        angle=  angle-360.d0*floor(angle/360.d0)
        RETURN
        END SUBROUTINE get_gm_angle 




