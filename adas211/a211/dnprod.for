CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/dnprod.for,v 1.2 2004/07/06 13:32:51 whitefor Exp $ Date $Date: 2004/07/06 13:32:51 $
CX
      SUBROUTINE DNPROD(A0,A,B0,B,C0,C,NMAX)
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
      DIMENSION A(20),B(20),C(20),W(20)
      N=0
      C0=A0*B0
1     N=N+1
      C(N)=A0*B(N)+B0*A(N)
      IF (N-NMAX) 1,2,2
2     IF (NMAX-1) 6,6,3
3     W(1)=1.0
      U=1.0
      DO 5 N=2,NMAX
      W(N)=1.0
      JMAX=N-1
      DO 4 J=1,JMAX
      V=W(J)
      W(J)=U+W(J)
      U=V
      J1=N-J
      C(N)=C(N)+W(J)*A(J1)*B(J)
4     CONTINUE
5     CONTINUE
6     RETURN
      END
