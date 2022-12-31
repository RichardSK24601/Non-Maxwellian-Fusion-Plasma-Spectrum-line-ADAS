CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/dnaq.for,v 1.4 2007/07/20 09:46:35 allan Exp $ Date $Date: 2007/07/20 09:46:35 $
CX
      SUBROUTINE DNAQ(A0,A,B0,B,Q,NMAX,JSWICH)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C      PURPOSE: GIVEN A0 AND ITS FIRST NMAX DERIVATIVES IN ARRAY A,
C      AND GIVEN Q AND NMAX, CALCULATES B0 AND ARRAY B, BEING THE
C      VALUE AND FIRST NMAX DERIVATIVES OF (A0)**Q
C
C      FOR PERHAPS GREATER SPEED, YOU MAY SET JSWICH TO 2 IF Q IS -2.0
C      JSWICH TO 3 IF Q IS -1.0
C      JSWICH TO 4 IF Q IS -0.5
C      JSWICH TO 5 IF Q IS -0.25
C      JSWICH TO 6 IF Q IS 0.25
C      JSWICH TO 7 IF Q IS 0.5
C      JSWICH TO 8 IF Q IS 2.0
C      OTHERWISE SET JSWICH TO 1
C-----------------------------------------------------------------------
C UNIX-IDL PORT:
C
C AUTHOR:  WILLIAM OSBORN (TESSELLA SUPPORT SERVICES PLC)
C
C DATE:    4TH JULY 1996
C
C VERSION: 1.1                          DATE: 04-07-96
C MODIFIED: WILLIAM OSBORN
C               - FIRST VERSION.
C
C VERSION: 1.2                          DATE: 19-12-01
C MODIFIED: Martin O'Mullane
C               - Removed junk from > column 72.
C
C VERSION: 1.3                          DATE: 16-05-07
C MODIFIED: Allan Whiteford
C               - Modified comments as part of subroutine documentation
C                 procedure.
C
C VERSION: 1.4                          DATE: 20-07-07
C MODIFIED: Allan Whiteford
C               - Further modification to comments as part of
C                 subroutine documentation procedure.
C
C-----------------------------------------------------------------------
       DIMENSION A(20),B(20),W(20)
       GO TO (1,2,3,4,5,6,7,8),JSWICH
    1  B0=A0**Q
       GO TO 20
    2  B0=1.0/(A0*A0)
       GO TO 20
    3  B0=1.0/A0
       GO TO 20
    4  B0=1.0/DSQRT(A0)
       GO TO 20
    5  B0=1.0/DSQRT(A0)
       B0=DSQRT(B0)
       GO TO 20
    6  B0=DSQRT(A0)
       B0=DSQRT(B0)
       GO TO 20
    7  B0=DSQRT(A0)
       GO TO 20
    8  B0=A0*A0
   20  CONTINUE
       W(1)=Q
      C=1.0/A0
      B(1)=Q*A(1)*B0*C
       DO 22 N=2,NMAX
      B(N)=Q*B0*A(N)
      W(N)=Q
      U=-1.0
      J1=N-1
       DO 21 J=1,J1
      V=W(J)
      W(J)=U+W(J)
      U=V
      J2=N-J
      B(N)=B(N)+W(J)*A(J)*B(J2)
   21  CONTINUE
      B(N)=C*B(N)
   22  CONTINUE
      RETURN
      END
