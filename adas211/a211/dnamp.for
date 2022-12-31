CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/dnamp.for,v 1.2 2004/07/06 13:32:39 whitefor Exp $ Date $Date: 2004/07/06 13:32:39 $
CX
       SUBROUTINE DNAMP(A0,A,E,ELL1,Z,X,NMAX,JMAX)
       IMPLICIT REAL*8(A-H,O-Z)
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
C MODIFIED: Martin O'MULLANE
C               - Removed junk from > column 72.
C
C-----------------------------------------------------------------------
       DIMENSION A(20),B(20),W(20),C(20)
       X1=1.0/X
       F=-2.0*Z*X1
       G=-ELL1*X1*X1
       W0=E+F+G
       S=1.0
       IF(W0-1.0D-50)5,6,6
    5  W0=-W0
       F=-F
       G=-G
       S=-S
    6  CONTINUE
       NMAX2=NMAX+2
       DO 1 N=1,NMAX2
       EN=N
       F=-EN*F*X1
       G=-(EN+1.0)*G*X1
       W(N)=F+G
    1  CONTINUE
       CALL DNAQ(W0,W,A0,A,-0.25D0,NMAX2,5)
       DO 4 J=1,JMAX
       CALL DNAQ(A0,A,B0,B,-1.0D0,NMAX,3)
       A0=A(2)
       DO 2 N=1,NMAX
       A(N)=A(N+2)
    2  CONTINUE
       CALL DNPROD(A0,A,B0,B,C0,C,NMAX)
       B0=W0+C0*S
       DO 3 N=1,NMAX
       B(N)=W(N)+C(N)*S
    3  CONTINUE
       CALL DNAQ(B0,B,A0,A,-0.25D0,NMAX,5)
    4  CONTINUE
       RETURN
      END
