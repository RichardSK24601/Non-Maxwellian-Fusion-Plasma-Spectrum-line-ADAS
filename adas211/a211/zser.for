CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/zser.for,v 1.2 2004/07/06 15:41:49 whitefor Exp $ Date $Date: 2004/07/06 15:41:49 $
CX
      SUBROUTINE ZSER(I,Z0,NSHELL,N,NUMEL,ALFA,ZS)
C      SERIES EXPANSION, ABOUT R=0, OF ZEFF.
C      INPUT AS IN ZEFF
C      OUTPUT IS ZS(J), J=1,100,
C      WHERE ZEFF=ZS(1)+ZS(2)*R+ZS(3)*R**2+...
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
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION N(10),NUMEL(10),ALFA(10),ZS(100),B(100),A(100)
       ZS(1)=Z0
       DO 1 J=2,100
       ZS(J)=0.0D0
    1  CONTINUE
       IF(NSHELL)14,14,15
   14  RETURN
   15  IF(I)16,22,22
   16  T2=(-Z0)**(1.0D0/3.0D0)
       DO 17 J=1,NSHELL
       T1=NUMEL(J)
       T=0.0D0
       A(1)=-1.0D60*T1
       B(1)=1.0D0
       DO 18 I1=2,60
       T=T+1.0D0
       A(I1)=-0.2075D0*A(I1-1)/T
       B(I1)=-1.19D0*B(I1-1)
   18  CONTINUE
       DO 19 I1=61,100
       A(I1)=0.0D0
       B(I1)=-1.19D0*B(I1-1)
   19  CONTINUE
       T3=T2*ALFA(J)
       T=1.0D0
       S=1.0D0
       DO 21 I1=2,100
       IF(S-1.0D50)24,23,23
   23  ZS(I1)=-ZS(I1-1)
       GO TO 21
   24  S=0.0D0
       DO 20 J2=1,I1
       J1=I1-J2+1
       S=S+A(J2)*B(J1)
   20  CONTINUE
       S=S*1.0D-60
       T=T*T3
       S=S*T
       ZS(I1)=ZS(I1)+S
   21  CONTINUE
   17  CONTINUE
       RETURN
   22  Z=-Z0
       DO 10 J=1,NSHELL
       T1=NUMEL(J)-1
       Z=Z-0.5D0*T1
       EN=N(J)
       RHO=2.0D0*Z*ALFA(J)/EN
       IMAX=N(J)+N(J)-1
       T2=EN+EN
       T3=1.0D0
       T4=1.0D0/T2
       B(1)=1.0D0
       DO 2 I1=1,IMAX
       T5=I1
       T2=T2-1.0D0
       T4=T4*RHO/T5
       B(I1+1)=T4*T2
    2  CONTINUE
       IF(I-J)3,4,3
    3  T1=NUMEL(J)
    4  A(1)=1.0D0
       DO 5 I1=1,99
       T=I1
       A(I1+1)=-A(I1)*RHO/T
       IF(DABS(A(I1+1))-1.0D-60)11,5,5
    5  CONTINUE
       GO TO 13
   11  I2=I1+1
       DO 12 I1=I2,100
       A(I1)=0.0D0
   12  CONTINUE
   13  DO 9 I1=3,99
       T=0.0D0
       JMAX=I1+1
       IF(I1-IMAX)7,7,6
    6  JMAX=IMAX+1
    7  CONTINUE
       DO 8 J1=1,JMAX
       J2=I1+2-J1
       T=T+B(J1)*A(J2)
    8  CONTINUE
       ZS(I1+1)=ZS(I1+1)-T1*T
    9  CONTINUE
       ZS(2)=ZS(2)+T1*RHO/(EN+EN)
       T1=NUMEL(J)+1
       Z=Z-0.5D0*T1
   10  CONTINUE
       RETURN
      END
