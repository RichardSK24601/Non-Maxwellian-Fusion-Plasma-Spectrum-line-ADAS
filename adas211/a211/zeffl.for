CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/zeffl.for,v 1.2 2004/07/06 15:41:19 whitefor Exp $ Date $Date: 2004/07/06 15:41:19 $
CX
      REAL*8 FUNCTION ZEFFL(I,Z0,NSHELL,N,NUMEL,ALFA,R,ZL,H,R1,
     &                      Z1,Z2,Z3)
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
       DIMENSION ZL(1000),N(10),NUMEL(10),ALFA(10)
       IF(R-R1)2,2,1
    1  T=1.0D0/(R*R)
       ZEFFL=Z1+T*(Z2+T*Z3)+ZEFF(I,Z0,NSHELL,N,NUMEL,ALFA,R)
       RETURN
    2  J=R/H
       T=J
       T=T*H
       J1=R1/H+0.5D0
       IF(J-J1+1)4,3,1
    3  J=J-1
       T=T-H
    4  IF(J-1)5,6,7
    5  J=1
       T=H
    6  T1=Z0
       GO TO 8
    7  T1=ZL(J-1)
    8  T2=ZL(J)
       T3=ZL(J+1)
       T4=ZL(J+2)
       P=(R-T)/H
       ZEFFL=P*(P-1.0D0)*((P+1.0D0)*T4-(P-2.0D0)*T1)/6.0D0
     1 +(P+1.0D0)*(P-2.0D0)*((P-1.0D0)*T2-P*T3)*0.5D0
       RETURN
      END
