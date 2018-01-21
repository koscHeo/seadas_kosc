        SUBROUTINE INVERT_3BY3(A, AINV,DET)
        IMPLICIT NONE
        REAL(8) A(3,3),AINV(3,3)
        REAL(8) Q1,Q2,Q3,DET
 
      Q1=A(2,2)*A(3,3)-A(3,2)*A(2,3)
      Q2=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      Q3=A(1,2)*A(2,3)-A(2,2)*A(1,3)
      DET=A(1,1)*Q1 + A(2,1)*Q2 + A(3,1)*Q3
      IF(ABS(DET).LT.1.E-30) RETURN

        AINV(1,1)=Q1
        AINV(1,2)=Q2
        AINV(1,3)=Q3

        AINV(2,1)=A(3,1)*A(2,3)-A(2,1)*A(3,3)
        AINV(2,2)=A(1,1)*A(3,3)-A(3,1)*A(1,3)
        AINV(2,3)=A(2,1)*A(1,3)-A(1,1)*A(2,3)

        AINV(3,1)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
        AINV(3,2)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
        AINV(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)

      AINV=AINV/DET

        RETURN
        END





