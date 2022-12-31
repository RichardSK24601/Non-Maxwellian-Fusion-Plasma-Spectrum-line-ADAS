CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/wig6j.for,v 1.4 2004/07/06 15:28:52 whitefor Exp $ Date $Date: 2004/07/06 15:28:52 $
CX
       REAL*8 FUNCTION WIG6J(A1,A2,A3,B1,B2,B3)
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
C VERSION: 1.2					DATE: 24-11-2000
C MODIFIED: MARTIN O'MULLANE & RICHARD MARTIN
C		CHANGED DIMENSIONS OF GAM & JGAM TO BE CONSISITENT WITH
C		DELTA2.FOR. NAMED COMMON BLOCK.
C
C VERSION  : 1.3                          
C DATE     : 19-12-2001
C MODIFIED : Martin O'Mullane
C               - Changed function defintion to a more standard form.
C               - Removed mainframe listing information beyond column 72.
C
C VERSION  : 1.4                         
C DATE     : 18-03-2003
C MODIFIED : Richard Martin
C               - Increased GAM, JGAM dimensions to 500.
C
C-----------------------------------------------------------------------
       DIMENSION GAM(500),JGAM(500),K(4),L(3)
       COMMON /FAC/GAM,JGAM
       CALL DELTA2(A1,A2,A3,D1,J1)
       CALL DELTA2(A1,B2,B3,D2,J2)
       CALL DELTA2(A2,B1,B3,D3,J3)
       CALL DELTA2(A3,B1,B2,D4,J4)
       D=D1*D2*D3*D4
       IF(D-1.0D-20)1,1,2
    1  WIG6J=0.0
       RETURN
    2  J=J1+J2+J3+J4
       D=DSQRT(D)
       K(1)=A1+A2+A3+0.001
       K(2)=A1+B2+B3+0.001
       K(3)=A2+B1+B3+0.001
       K(4)=A3+B1+B2+0.001
       L(1)=A1+B1+A2+B2+0.001
       L(2)=A2+B2+A3+B3+0.001
       L(3)=A1+B1+A3+B3+0.001
       NMIN=K(1)
       DO 4 I=2,4
       IF(K(I)-NMIN)4,4,3
    3  NMIN=K(I)
    4  CONTINUE
       NMAX=L(1)
       DO 6 I=2,3
       IF(L(I)-NMAX)5,6,6
    5  NMAX=L(I)
    6  CONTINUE
       D=D*(-1.0)**NMIN
       N1=NMIN-K(1)
       N2=NMIN-K(2)
       N3=NMIN-K(3)
       N4=NMIN-K(4)
       N5=L(1)-NMIN
       N6=L(2)-NMIN
       N7=L(3)-NMIN
       D=D*GAM(NMIN+2)/(GAM(N1+1)*GAM(N2+1)*GAM(N3+1)*GAM(N4+1)
     1 *GAM(N5+1)*GAM(N6+1)*GAM(N7+1))
       M1=JGAM(NMIN+2)-JGAM(N1+1)-JGAM(N2+1)-JGAM(N3+1)-JGAM(N4+1)
     1  -JGAM(N5+1)-JGAM(N6+1)-JGAM(N7+1)
       D=D*8.0**(J+M1+M1)
       IF(NMAX-NMIN)7,7,8
    7  WIG6J=D
       RETURN
    8  DSUM=1.0
       DT=1.0
       DP1=N1
       DP2=N2
       DP3=N3
       DP4=N4
       DP5=N5+1
       DP6=N6+1
       DP7=N7+1
       DP8=NMIN+1
       NMAX1=NMAX-NMIN
       DO 9 N=1,NMAX1
       DP1=DP1+1.0
       DP2=DP2+1.0
       DP3=DP3+1.0
       DP4=DP4+1.0
       DP5=DP5-1.0
       DP6=DP6-1.0
       DP7=DP7-1.0
       DP8=DP8+1.0
       DT=-DT*DP5*DP6*DP7*DP8/(DP1*DP2*DP3*DP4)
       DSUM=DSUM+DT
    9  CONTINUE
       SUM=DSUM
       D=D*SUM
       GO TO 7
      END
