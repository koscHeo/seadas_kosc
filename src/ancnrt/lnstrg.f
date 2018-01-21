      integer function LNSTRG(STRING)
C**********************************************************************
C
C     This function returns with its value equal to the position in
C     STRING of the last non-blank and non-null character.  If STRING
C     is all blanks or null, LNSTRG will equal zero.
C
C     Created by Michael Darzi, GSC, 10/86.
C               To UNIX: BD Schieber, SAIC/GSC, 93
C**********************************************************************
      implicit none
      character STRING*(*), BLNK*(*), NULL
      byte BNULL, BTMP
      parameter (BLNK=' ', BNULL='0'o)
      data BTMP/BNULL/
      equivalence (NULL,BTMP)

      LNSTRG = LEN(STRING)
   10 continue
      if (LNSTRG .lt. 1) return
      if ((STRING(LNSTRG:LNSTRG) .ne. BLNK) .and.
     .    (STRING(LNSTRG:LNSTRG) .ne. NULL)) return
      LNSTRG = LNSTRG - 1
      go to 10
      end
