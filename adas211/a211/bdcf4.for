CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/bdcf4.for,v 1.2 2004/07/06 11:37:29 whitefor Exp $ Date $Date: 2004/07/06 11:37:29 $
CX
       REAL*8 FUNCTION BDCF4(E,N,L,Z,X)
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
C-----------------------------------------------------------------------
       IMPLICIT REAL*8(A-H,O-Z)
       EL=L
       GNU=-Z/DSQRT(-E)
       T1=GNU+EL+1.0
       T2=GNU-EL
       IF(T2-0.5)4,4,5
    4  T2=1.0
    5  T3=GAMA6(T1)*GAMA6(T2)
       B=DSQRT(-Z/T3)/GNU
       B=B*(-1.0)**(N-L+1)
       R=-Z*X
       T=2.0*R/GNU
       T3=1.0/T
       S=1.0
       A=1.0
       K0=GNU+GNU+1.5+T
       TH=K0
       TH=GNU+GNU+1.0+T-TH
       DO 1 K=1,K0
       TK=K
       A=A*(TK-GNU+EL)*(GNU-TK+EL+1.0)*T3/TK
       T1=DABS(A)
       T2=DABS(S)
       IF(T2-T1*1.0D7)1,1,2
    1  S=S+A
       T1=(GNU-EL)*(GNU+EL+1.0)
       T2=-2.0*GNU*T1
       C=-0.5+0.125*T3*((2.0*TH-1.0)+(TH*TH-1.5*TH+0.25-2.0*T1)*T3)
       S=S+C*A
    2  BDCF4=(T**GNU)*S*DEXP(-(0.5*T))*B
       RETURN
      END
