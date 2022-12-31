CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/phase.for,v 1.3 2007/07/20 09:46:35 allan Exp $ Date $Date: 2007/07/20 09:46:35 $
CX
      REAL*8 FUNCTION PHASE(E,EL,Z,X)
C-----------------------------------------------------------------------
C
C PURPOSE: CALCULATES ASYMPTOTIC PHASE OF FREE COULOMB
C          REAL FUNCTION FCF4
C
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
C VERSION  : 1.2                          
C DATE     : 19-12-2001
C MODIFIED : Martin O'Mullane
C               - Changed function defintion to a more standard form.
C               - Removed mainframe listing information beyond column 72.
C
C VERSION  : 1.3
C DATE     : 20-07-2007
C MODIFIED : Allan Whiteford
C               - Small modification to comments to allow for automatic
C                 documentation preparation.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      L=EL+0.5
      PI=3.141592654
      ZZ=Z*Z
      CK=DSQRT(E)
       XK=X*CK
      XZ=X*Z
      C=EL*(EL+1.0)
      CHI=DSQRT(XK*XK-XZ-XZ-C)
      C1=1.0/CHI
      IF(C) 1,1,2
    1 THETA=0.25/(CHI+XK)
      GO TO 8
    2 T=CK*C*CHI+ZZ*X+C*Z
      T1=DABS (T)
      IF(T1-1.0D-15)3,3,4
    3 THETA =0.5*PI
      GO TO 7
    4 ARG=DSQRT(C*(ZZ+E*C)*(CK*C*CHI+E*C*X+2.*ZZ*X+C*Z)/(CK*CHI+E*X-
     1Z))/T
      IF(T)5,5,6
    5 THETA =PI-DATAN (-ARG)
      GOTO 7
    6 THETA=DATAN (ARG)
    7  THETA=THETA*(C+0.125)/DSQRT(C)
    8 IF(E-1.0D-50) 9,9,10
    9 PHASE=CHI+CHI+THETA-(EL+0.25)*PI-0.166666667*C1*(1.0+1.25*(XZ+C)*C
     11*C1)
      GO TO 11
   10 A=Z/CK
      B=E*C+ZZ
      D=3.0*(CHI+XK)*B+ZZ*(CHI-XK)-CK*C*Z
      D=(D/(24.0*B*(CHI+XK))+0.208333333*(XZ+C)*C1*C1)*C1
      AG=ARGAM(L,A)
      PHASE=CHI-A*DLOG(CHI+XK-A)+AG+A-0.5*PI*EL+THETA-D
   11 RETURN
      END
