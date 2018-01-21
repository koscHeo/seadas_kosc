C ***********************************************************************

      SUBROUTINE GET_SOA_PROD(idProd,pix,val)

      IMPLICIT NONE

#include "l12_parms.h"
#include "l2prod.h"

      INTEGER*4 idProd,pix
      REAL*4    val

      INTEGER*4 NSCANSZ
      PARAMETER (NSCANSZ=MAXPIX)
      REAL*4    chl_out(NSCANSZ),acdm_out(NSCANSZ),bbp_out(NSCANSZ),
     &          chl4s_out(NSCANSZ),pcentcdm_out(NSCANSZ)
      COMMON /GS97_OUT/ chl_out,acdm_out,bbp_out,chl4s_out,
     &                  pcentcdm_out
      REAL*4    w0_out(NSCANSZ)
      COMMON /GENERIC_OUT/ w0_out
      REAL*4    v_out(NSCANSZ)
      COMMON /SOA_OUT/ v_out

      INTEGER j,isDebug

      isDebug=0

      IF (idProd.EQ.CAT_chl_soa) THEN
         val=chl_out(pix)
      ELSE IF (idProd.EQ.CAT_adg_soa) THEN
         val=acdm_out(pix)
      ELSE IF (idProd.EQ.CAT_bbp_soa) THEN
         val=bbp_out(pix)
      ELSE IF (idProd.EQ.CAT_pcentcdm_soa) THEN
         val=pcentcdm_out(pix)
      ELSE IF (idProd.EQ.CAT_w0_soa) THEN
         val=w0_out(pix)
      ELSE IF (idProd.EQ.CAT_v_soa) THEN
         val=v_out(pix)
      ELSE
         WRITE(*,*) 'Unknown SOA product ID'
         CALL EXIT(1)
      ENDIF

      IF(isDebug.EQ.1) THEN
        DO j=1,NSCANSZ
           PRINT*,bbp_out(j), w0_out(j)
        ENDDO
      ENDIF

      RETURN
      END

C ***********************************************************************

C ***********************************************************************

      SUBROUTINE GET_SMA_PROD(idProd,pix,val)

      IMPLICIT NONE

#include "l12_parms.h"
#include "l2prod.h"

      INTEGER*4 idProd,pix
      REAL*4    val

      INTEGER*4 NSCANSZ
      PARAMETER (NSCANSZ=MAXPIX)
      REAL*4    chl_out(NSCANSZ),acdm_out(NSCANSZ),bbp_out(NSCANSZ),
     &        chl4s_out(NSCANSZ),pcentcdm_out(NSCANSZ)
      COMMON /GS97_OUT/ chl_out,acdm_out,bbp_out,chl4s_out,
     &                  pcentcdm_out
      REAL*4    w0_out(NSCANSZ)
      COMMON /GENERIC_OUT/ w0_out
      INTEGER*4 dom_out(NSCANSZ,4)
      COMMON /SMA_OUT/ dom_out

      INTEGER j,isDebug

      isDebug=0

      IF (idProd.EQ.CAT_chl_sma) THEN
         val=chl_out(pix)
      ELSE IF (idProd.EQ.CAT_adg_sma) THEN
         val=acdm_out(pix)
      ELSE IF (idProd.EQ.CAT_bbp_sma) THEN
         val=bbp_out(pix)
      ELSE IF (idProd.EQ.CAT_w0_sma) THEN
         val=w0_out(pix)
      ELSE IF (idProd.EQ.CAT_dom_sma) THEN
         val=real(dom_out(pix,1))
      ELSE
         WRITE(*,*) 'Unknown SMA product ID'
         CALL EXIT(1)
      ENDIF

      IF(isDebug.EQ.1) THEN
        DO j=1,NSCANSZ
           PRINT*,bbp_out(j), w0_out(j)
        ENDDO
      ENDIF

      RETURN
      END

C ***********************************************************************

C ***********************************************************************

      SUBROUTINE CLEAR_SOA_PROD

      IMPLICIT NONE


      INTEGER*4 NSCANSZ
      PARAMETER (NSCANSZ=MAXPIX)
      REAL*4    chl_out(NSCANSZ),acdm_out(NSCANSZ),bbp_out(NSCANSZ),
     &          chl4s_out(NSCANSZ),domtot_out(NSCANSZ)
      COMMON /GS97_OUT/ chl_out,acdm_out,bbp_out,chl4s_out,domtot_out
      REAL*4    w0_out(NSCANSZ)
      COMMON /GENERIC_OUT/ w0_out
      REAL*4    v_out(NSCANSZ)
      COMMON /SOA_OUT/ v_out
      INTEGER*4 dom_out(NSCANSZ,4)
      COMMON /SMA_OUT/ dom_out


      INTEGER j,jj

      DO j=1,NSCANSZ
         chl_out(j)     = BAD_FLT
         acdm_out(j)    = BAD_FLT
         bbp_out(j)     = BAD_FLT
         chl4s_out(j)   = BAD_FLT
         domtot_out(j)  = BAD_FLT
         w0_out(j)      = BAD_FLT
         v_out(j)       = BAD_FLT
      ENDDO


      RETURN
      END

C ***********************************************************************

C ***********************************************************************

      SUBROUTINE CLEAR_SMA_PROD

      IMPLICIT NONE

      INTEGER*4 NSCANSZ
      PARAMETER (NSCANSZ=MAXPIX)
      REAL*4    chl_out(NSCANSZ),acdm_out(NSCANSZ),bbp_out(NSCANSZ),
     &        chl4s_out(NSCANSZ),pcentcdm_out(NSCANSZ)
      COMMON /GS97_OUT/ chl_out,acdm_out,bbp_out,chl4s_out,
     &                  pcentcdm_out
      REAL*4    w0_out(NSCANSZ)
      COMMON /GENERIC_OUT/ w0_out
      INTEGER*4 dom_out(NSCANSZ,4)
      COMMON /SMA_OUT/ dom_out

      INTEGER j,jj

      DO j=1,NSCANSZ
         chl_out(j)     = BAD_FLT
         acdm_out(j)    = BAD_FLT
         bbp_out(j)     = BAD_FLT
         chl4s_out(j)   = BAD_FLT
         w0_out(j)      = BAD_FLT
         DO jj=1,4
            dom_out(j,jj) = BAD_INT
         ENDDO
      ENDDO


      RETURN
      END

C ***********************************************************************
