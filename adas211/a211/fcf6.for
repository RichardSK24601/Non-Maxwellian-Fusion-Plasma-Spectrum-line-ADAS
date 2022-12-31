CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/fcf6.for,v 1.2 2004/07/06 13:52:03 whitefor Exp $ Date $Date: 2004/07/06 13:52:03 $
CX
       SUBROUTINE FCF6(F,C,A2,AMP,PHI,DEL,X0,N,L,E,JSN,Z0,NSHELL,NC,
     1 NUMEL,ALFA,ZL,Z1,Z2,Z3,ZS,X1,H,X2,H2)
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
C
C-----------------------------------------------------------------------
C
       DIMENSION F(1000),A(100),ZL(1000),ZS(100),AMP(20),AMP3(20)
       DIMENSION R(1000),NC(10),NUMEL(10),ALFA(10)
       XV=DABS(ZS(2)/ZS(1))
       ZA=Z0+Z1
       IF(NSHELL)108,108,106
  106  DO 107 J=1,NSHELL
       T=NUMEL(J)
       IF(J-JSN)110,109,110
  109  T=NUMEL(J)-1
  110  ZA=ZA+T
  107  CONTINUE
  108  CONTINUE
       NMAX=10
       JMAX=5
       TMAX=0.015625
       I1=0.5+X1/H
       X1=I1
       X1=X1*H
       EL=L
       L1=L+1
       S=ZS(2)
       ZS(2)=S-0.5D0*E
       W2=-EL*(EL+1.0)
       W3=0.25-W2
       Z=Z0
       M=0
       M1=1
       X0=W3/(DSQRT(Z*Z+W3*E)-Z)
       IF(XV*X0-0.35D0)2,2,1
    1  X0=0.35D0/XV
    2  CONTINUE
       H0=-1.0
       IF(X1-X0)3,4,4
    3  X0=X1
    4  IF(X0-2.0*H)5,13,13
    5  IF(X0-H)6,12,12
    6  X00=H
    7  X00=0.5*X00
       IF(X0-X00)7,8,8
    8  H0=X00
       T1=E-2.0D0*ZEFFL(JSN,Z0,NSHELL,NC,NUMEL,ALFA,X00,ZL,H,X1,Z1,Z2,
     2 Z3)/X00
       T1=DABS(T1+W2/(X00*X00))
    9  IF(T1*H0*H0-TMAX)11,11,10
   10  H0=0.5*H0
       GO TO 9
   11  X0=X00+H0
       JH=(H-X0)/H0+0.5
       I=0
       GO TO 14
   12  X0=H+H
   13  I0=0.5+X0/H
       X0=I0
       X0=X0*H
   14  A(1)=1.0
       J=0
   15  J=J+1
       IF(J-99)17,17,16
   16  X0=0.49*X0
       GO TO 2
   17  CONTINUE
       T=0.0
       DO 18 K=1,J
       K1=J+1-K
       T=X0*T+A(K)*ZS(K1)
   18  CONTINUE
       TJ=J
       A(J+1)=2.0D0*X0*T/((EL+EL+TJ+1.0)*TJ)
       IF(J-5)15,15,19
   19  C2=DABS(A(J+1))+DABS(A(J))
       IF(C2-1.0D6)20,16,16
   20  IF(C2-1.0D-20)21,21,15
   21  ZS(2)=S
       J0=J+1
       IF(H0)30,22,22
   22  X=X0
       J=J0
       T1=A(J)
   23  J=J-1
       T1=A(J)+T1
       IF(J-1)24,24,23
   24  F2=T1*X**L1
       IF(F2)25,26,26
   25  M2=-1
       GO TO 27
   26  M2=1
   27  IF(M1+M2)29,28,29
   28  M=M+1
       M1=M2
   29  CONTINUE
       V2=2.0D0*ZEFFL(JSN,Z0,NSHELL,NC,NUMEL,ALFA,X,ZL,H,X1,Z1,Z2,Z3)/X
       GO TO 49
   30  X=0.0
       DO 38 I=1,I0
       X=X+H
       J=J0
       T=X/X0
       T1=A(J)
   31  J=J-1
       T1=A(J)+T*T1
       IF(J-1)32,32,31
   32  F(I)=T1*X**L1
       IF(T1)33,34,34
   33  M2=-1
       GO TO 35
   34  M2=1
   35  IF(M1+M2)37,36,37
   36  M=M+1
       M1=M2
   37  CONTINUE
   38  CONTINUE
       IF(I0-I1)39,69,69
   39  I=I0
       X=X0
       T=DABS(E-2.0D0*ZL(I-1)/(X-H)+W2/((X-H)*(X-H)))*H*H
       IF(T-TMAX)40,40,47
   40  HH=H*H
       H1=0.0833333333*HH
       W0=E-2.0D0*ZL(I)/X+W2/(X*X)
       W1=E-2.0D0*ZL(I-1)/(X-H)+W2/((X-H)*(X-H))
       C0=F(I)*(1.0+(H1-HH)*W0)-F(I-1)*(1.0+H1*W1)
       C1=F(I)*(1.0+H1*W0)
   41  I=I+1
       X=X+H
       C1=C1+C0
       C2=E-2.0D0*ZL(I)/X+W2/(X*X)
       F(I)=C1/(1.0+H1*C2)
       IF(C1)42,43,43
   42  M2=-1
       GO TO 44
   43  M2=1
   44  IF(M1+M2)46,45,46
   45  M=M+1
       M1=M2
   46  CONTINUE
       C0=C0-HH*C2*F(I)
       IF(I-I1)41,70,70
   47  H0=H
       F2=F(I)
       V2=2.0D0*ZL(I)/X
       JH=1
   48  H0=0.5*H0
       JH=JH+JH
       T=DABS(E-2.0D0*ZL(I-1)/(X-H)+W2/((X-H)*(X-H)))*H0*H0
       IF(T-TMAX)49,48,48
   49  J=J0
       T1=A(J)
       T=(X-H0)/X0
   50  J=J-1
       T1=A(J)+T*T1
       IF(J-1)51,51,50
   51  F1=T1*(X-H0)**L1
   52  HH=H0*H0
       H1=0.0833333333*HH
       W0=E-V2+W2/(X*X)
       T=X-H0
       T1=1.0/T
       V1=2.0D0*ZEFFL(JSN,Z0,NSHELL,NC,NUMEL,ALFA,T,ZL,H,X1,Z1,Z2,Z3)/T
       W1=E-V1+W2*T1*T1
       C0=F2*(1.0+(H1-HH)*W0)-F1*(1.0+H1*W1)
       C1=F2*(1.0+H1*W0)
   53  J=0
       T=X
   54  J=J+1
       T=T+H0
       V0=V1
       V1=V2
       C1=C1+C0
       T1=1.0/T
       V2=2.0D0*ZEFFL(JSN,Z0,NSHELL,NC,NUMEL,ALFA,T,ZL,H,X1,Z1,Z2,Z3)/T
       C2=E-V2+W2*T1*T1
       F0=F1
       F1=F2
       F2=C1/(1.0+H1*C2)
       IF(F2)55,56,56
   55  M2=-1
       GO TO 57
   56  M2=1
   57  IF(M1+M2)59,58,59
   58  M=M+1
       M1=M2
   59  CONTINUE
       C0=C0-HH*C2*F2
       IF(J-JH)54,60,60
   60  I=I+1
       X=T
       F(I)=F2
       IF(I-1)61,61,62
   61  JH=0.5+H/H0
   62  CONTINUE
       IF(I-I1)63,71,71
   63  H0=H0+H0
       T=DABS(E-V0+W2/((X-H0)*(X-H0)))*H0*H0
       IF(T-TMAX)65,64,64
   64  H0=0.5*H0
       GO TO 53
   65  JH=0.5+H/H0
       IF(JH-1)66,66,68
   66  IF(I-1)67,67,40
   67  H0=0.5*H0
       JH=2
       GO TO 53
   68  F1=F0
       GO TO 52
   69  HH=H*H
       H1=0.0833333333*HH
       W0=E-2.0D0*ZL(I1)/X1+W2/(X1*X1)
       W1=E-2.0D0*ZL(I1-1)/(X1-H)+W2/((X1-H)*(X1-H))
       C0=F(I1)*(1.0+(H1-HH)*W0)-F(I1-1)*(1.0+H1*W1)
       C1=F(I1)*(1.0+H1*W0)
   70  H2=H
       GO TO 72
   71  H2=H0
   72  F2=F(I1)
       X2=X1
       IF(E+1.0D-40)73,95,95
   73  X12=0.2/(DSQRT(-E))
       I12=X12/H2
       IF(I12-2)74,75,75
   74  I12=2
   75  T=I12
       X2=X1+T*H2
       X=X1
       XT=X1
       R(1)=E-2.0D0*ZL(I1)/X+W2/(X*X)
       I=1
       I12=I12+1
       FI1=F2
   76  X=X+H2
       I=I+1
       T=1.0/X
       T3=2.0D0*ZEFFL(JSN,Z0,NSHELL,NC,NUMEL,ALFA,X,ZL,H,X1,Z1,Z2,Z3)/X
       T1=2.0*ZA*T-T3
       T2=-E+T3-W2*T*T
       C2=-T2
       IF(I-1000)77,77,78
   77  R(I)=C2
   78  CONTINUE
       IF(C2)80,79,79
   79  XT=X
   80  CONTINUE
       IF(I-I12)81,81,82
   81  C1=C1+C0
       F2=C1/(1.0+H1*C2)
       C0=C0-HH*C2*F2
       GO TO 76
   82  IF(DABS(T2)-1.0D4*DABS(T1))76,83,83
   83  FI2=F2
       X3=X
       X=X3+H2+H2
       T=1.0/X
       W1=E-2.0*ZA*T+W2*T*T
       R1=BDCF4(E,N,L,ZA,X)
       X=X3+H2
       T=1.0/X
       W0=E-2.0*ZA*T+W2*T*T
       R2=BDCF4(E,N,L,ZA,X)
       C0=R2*(1.0+(H1-HH)*W0)-R1*(1.0+H1*W1)
       C1=R2*(1.0+H1*W0)
       I=I+1
   84  X=X-H2
       I=I-1
       C1=C1+C0
       IF(I-1000)85,85,86
   85  C2=R(I)
       GO TO 87
   86  T=1.0/X
       T3=2.0D0*ZEFFL(JSN,Z0,NSHELL,NC,NUMEL,ALFA,X,ZL,H,X1,Z1,Z2,Z3)/X
       C2=E-T3+W2*T*T
   87  CONTINUE
       R2=C1/(1.0+H1*C2)
       C0=C0-HH*C2*R2
       IF(I-I12)89,88,89
   88  RI2=R2
   89  IF(I-1)90,90,84
   90  RI1=R2
       IF(N-L-1-M)91,92,93
   91  DEL=-1.0D0
       GO TO 94
   92  C1=RI1*RI2*(FI2*RI1-FI1*RI2)
       C2=FI1*RI1-FI2*RI2
       DEL=0.636619772D0*DATAN2(C1,C2)
       GO TO 94
   93  DEL=2.0D0
   94  CONTINUE
       A2=0.0
       PHI=M
       X1=XT
       X2=X3
       C=1.0
       RETURN
   95  X2=X2+H2
       C1=C1+C0
       T=1.0/X2
       T1=T*T
       V2=2.0D0*ZEFFL(JSN,Z0,NSHELL,NC,NUMEL,ALFA,X2,ZL,H,X1,Z1,Z2,Z3)/
     3 X2
       T2=V2-2.0D0*ZA*T
       C2=E-V2+W2*T1
       F2=C1/(1.0+H1*C2)
       IF(F2)96,97,97
   96  M2=-1
       GO TO 98
   97  M2=1
   98  IF(M1+M2)100,99,100
   99  M=M+1
       M1=M2
  100  CONTINUE
       C0=C0-HH*C2*F2
       CHI2=C2*X2*X2+W2+W2
       IF(CHI2-60.0D0)95,101,101
  101  IF(C2-1.0D4*DABS(T2))95,102,102
  102  ELL1=-W2
       CALL DNAMP(A2,AMP,E,ELL1,ZA,X2,NMAX,JMAX)
       X23=1.5708*A2*A2
       X3=X2+X23
       X=X2
  103  X=X+H2
       C1=C1+C0
       T=1.0/X
       C2=E+T*(W2*T-2.0*ZA)
       F3=C1/(1.0+H1*C2)
       C0=C0-HH*C2*F3
       IF(X-X3)103,104,104
  104  X3=X
       CALL DNAMP(A3,AMP3,E,ELL1,ZA,X3,NMAX,JMAX)
       PHI2=PHASE(E,EL,ZA,X2)
       PHI3=PHASE(E,EL,ZA,X3)
       PHI=PHI2
       F2=F2/A2
       F3=F3/A3
       S2=DSIN(PHI2)
       S3=DSIN(PHI3)
       C2=DCOS(PHI2)
       C3=DCOS(PHI3)
       S23=DSIN(PHI3-PHI2)
       C23=DCOS(PHI3-PHI2)
       C=DSQRT(S23*S23/(F2*F2+F3*F3-2.0*F2*F3*C23))
       T=C/S23
       SD=(F2*S3-F3*S2)*T
       CD=(F3*C2-F2*C3)*T
       DEL=DATAN2(SD,CD)
       M1=(PHI+DEL)/3.1415926536
       EM1=M-M1
       DEL=EM1+DEL/3.1415926536
       DO 105 I=1,I1
       F(I)=C*F(I)
  105  CONTINUE
       RETURN
      END
