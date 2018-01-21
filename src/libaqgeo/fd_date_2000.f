      SUBROUTINE FD_DATE_2000(TIME_2000, SECYR,LYEAR,IDAYJL,IMON,IDAYMO,SECDY)
        IMPLICIT NONE
C     TIME_2000= SEC FROM BEGIN OF 2000
C     ISECYR= SEC FROM BEGIN OF LYEAR
C     LYEAR= YEAR, FULL FOR DIGITS, I.E., 1987
C     IDAYJL= JULIAN DAT OF YEAR
C     IMON= MONTH (1 TO 12)                   
C     IDAYMO= DAY OF MONTH (1 TO 31)
C     ISECDY= SECONDS OF DAY (0 TO 86399)
      
        REAL(8) TIME_2000,SECYR,SECDY,XTIME 
        INTEGER(4) LYEAR,IDAYJL,IMON,IDAYMO

        INTEGER(4) IDAYTOT,IDAYBG,ILEAP,JMON,IDAYFX(12,2) 
                                                   
      DATA IDAYFX/1,32,60,91,121,152,182,213,244,274,305,335,           
     1            1,32,61,92,122,153,183,214,245,275,306,336/  
     
      XTIME=TIME_2000 + 410227200.D0 !CONVERT TO TIME87 
        
      IF(XTIME.LT.0 .OR. XTIME.GT.3534451200.D0)  then
         write(*,*) 'TIME OOB IN FD_DATE_2000: ', xtime !3534451200 IS JAN 1 2099       
         call exit(1)
      endif

      IDAYTOT=1 + INT(XTIME/86400.D0)
      LYEAR=1987 + INT((IDAYTOT-1)/365)                        !KYEAR MAY BE 1 YEAR TOO BIG
      IDAYBG= 1 + 365*(LYEAR-1987) + INT((LYEAR-1985)/4)       !BEGIN DAY OF KYEAR
      IF(IDAYTOT.LT.IDAYBG) LYEAR=LYEAR-1
  
      SECYR=XTIME-31536000.D0*(LYEAR-1987)-86400.D0*INT((LYEAR-1985)/4)
                                                  
      ILEAP=1                                                           
      IF(LYEAR.EQ.4*INT(LYEAR/4)) ILEAP=2                               
      IDAYJL=1+INT(SECYR/86400.D0)                                           
      DO JMON=2,12                                                  
      IMON=JMON-1                                                       
      IF(IDAYFX(JMON,ILEAP).GT.IDAYJL) GOTO 210  
        ENDDO                        
      IMON=12                                                         
  210 CONTINUE                                                          
      IDAYMO=IDAYJL-IDAYFX(IMON,ILEAP)+1                                  
      SECDY=SECYR-(IDAYJL-1)*86400.D0                                         
      RETURN                                                            
      END                                                               
