      SUBROUTINE WAGAUS(A,B,AGS,WGS,NGAUSS)
C   ******************************************************************
C   *                                                                *
C   *    Obtain abscissas and weight factors for Gauss-Legendre      *
C   *    n-point integration in the range [A,B]. To see abscissas    *
C   *    and weigths in the reduced range [-1,1] set IDBG = 1.       *
C   *                                                                *
C   *    A     : Lower integration limit                             *
C   *    B     : Upper integration limit                             *
C   *    NGAUSS: Number of points in the range [A,B]                 *
C   *    AGS   : Abscissas in the range [A,B]                        *
C   *    WGS   : Corresponding weight factors                        *
C   *    EPS   : Floating precission in the root-finding             *
C   *                                                                *
C   *    HLS: 30-JAN-87                                              *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ONE=1,EPS=1.D-15,IDBG=0)
      DIMENSION AGS(NGAUSS),WGS(NGAUSS)
  101 FORMAT('0',10X,'Mesh and weights found in WAGAUSS',//,3X,'I',9X,
     1       'X',17X,'WEIGHT',/)
  102 FORMAT(' ',I3,2F20.15)
  103 FORMAT('0WAGAUS: EPS too small for machine precission',/,8X,
     1       'EPS = ',E15.8)
C
      PI=ACOS(-ONE)
      NHALF=(NGAUSS+1)/2
C
C     Guess at the I'th root of the Legendre polynomial and refine
C     by Newton to obtain abscissas and weights in the normalized
C     range [-1,1]
C
      DO 11 I=1,NHALF
      Z=COS(PI*(I-.25D0)/(NGAUSS+0.5D0))
      NEWTON=0
   10    NEWTON=NEWTON+1
         IF(NEWTON.GT.20) THEN
            WRITE(6,103) EPS
            STOP
         ENDIF
         P1=1.D0
         P2=0.D0
         DO 12 J=1,NGAUSS
         P3=P2
         P2=P1
         P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
   12    CONTINUE
C
C        Newton root-finding using the derivative PP
C
         PP=NGAUSS*(Z*P1-P2)/(Z*Z-1.D0)
         Z1=Z
         Z=Z1-P1/PP
      IF(ABS(Z-Z1).GT.EPS) GO TO 10
      AGS(I)=-Z
      AGS(NGAUSS+1-I)=Z
      WGS(I)=2.D0/((1.D0-Z*Z)*PP*PP)
      WGS(NGAUSS+1-I)=WGS(I)
   11 CONTINUE
      IF(IDBG.EQ.1) THEN
         WRITE(6,101)
         DO 16 I=1,NGAUSS
   16    WRITE(6,102) I,AGS(I),WGS(I)
      ENDIF
C
C     Transform abscissas and weights to the range [A,B]
C
      C=(B-A)/2.0D0
      D=(B+A)/2.0D0
      DO 20 I=1,NGAUSS
      AGS(I)=C*AGS(I)+D
   20 WGS(I)=C*WGS(I)
      RETURN
      END
