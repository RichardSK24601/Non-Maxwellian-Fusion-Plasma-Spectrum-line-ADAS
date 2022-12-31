CX UNIX PORT - SCCS info: Module @(#)aphotdw.for	1.2 Date 03/15/02
CX
       SUBROUTINE APHOTDW(B,B1,V,N,L,L1,LP,ISP,LT,LT1,IS,PREC,PION,
     &PSTIM,IRESOL,ndgnt,gaunt,energy,ie)
       IMPLICIT REAL*8(A-H,O-Z)
C ---------------------------------------------------------------------
C
C  VERSION OF PHOTDW FOR USE BY ADASRRC WHICH USES AGIIDW
C  ******************* H.P.SUMMERS, JET   30 JUNE 1992 ***************
C
C  PURPOSE: CALCULATE PHOTO INTEGRALS USING GIIDW BOUND-FREE
C           GAUNT-FACTORS
C
C  SAME AS RECOM.FORT(PHOTO5) BUT CALLS GIIDW
C  ******************* H.P.SUMMERS, JET   19 AUG. 1984****************
C  INPUT
C      B=1.5789D5*Z**2/(V**2*TE)
C      B1=1.5789D5*Z**2/(V**2*TR)
C               WHERE
C               TE=ELECTRON TEMPERATURE (K)
C               TR=RADIATION TEMPERATURE (K)
C               Z=BOUND STATE ION CHARGE +1
C               (THUS Z**2/V***2 IS THE IONISATION POTENTIAL  (RYD))
C      V=EFFECTIVE PRINCIPAL QUANTUM NUMBER OF BOUND ELECTRON
C      N=PRINCIPAL QUANTUM NUMBER OF BOUND ELECTRON
C      L=ORBITAL QUANTUM NUMBER OF BOUND ELECTRON
C      L1=ORBITAL QUANTUM NUMBER OF FREE ELECTRON
C      ISP=2*SP+1 WHERE SP IS TOTAL SPIN OF PARENT STATE
C      LP=TOTAL ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER OF PARENT STATE
C      LT=TOTAL ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER OF BOUND SYSTEM
C      LT1=TOTAL ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER OF FREE SYSTEM
C      IS=2*S+1 WHERE S IS TOTAL SPIN OF SYSTEM
C      ndgnt = max number of Gaunt factors allowed
C  OUTPUT
C      PREC=RADIATIVE RECOMBINATION INTEGRAL
C      PION=PHOTOIONISATION INTEGRAL
C      PSTIM=STIMULATED RECOMBINATION INTEGRAL
C          WHERE
C      IRESOL=1 FOR  ((LP,SP)N L LT S,(LP,SP)L1 LT1 S)
C            =2 FOR  ((LP,SP)N L LT S,(LP,SP)L1 S) =ABOVE LT1 SUM
C            =3 FOR  ((LP,SP)N L S,(LP,SP)L1 S)  = ABOVE LT SUM
C            =4 FOR  ((LP,SP)N L,(LP,SP)L1)  = ABOVE S SUM
C            =5 FOR  NO L RESOLUTION USING GBF
C      gaunt()  = Bound-free Gaunt factor at energy
C      energy() = v**2*e
C                 where e = (free electron energy)/z**2 (ryd)
C                       v = effective principal quantum number 
C                           of bound electron
C      ie       = number of Gaunt/energy pairs
C
C
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
C VERSION: 1.3                          DATE: 02-02-05
C MODIFIED: Paul Bryans
C               - Returns Gaunt factor, associated vve and number
C                 of Gaunt/vve pairs
C
C VERSION: 1.4                          DATE: 16-05-07
C MODIFIED: Allan Whiteford
C               - Modified comments as part of subroutine documentation
C                 procedure.
C
CC ----------------------------------------------------------------------
C
       DIMENSION X1(8),X2(8),W1(8),W2(8)
       dimension gaunt(ndgnt), energy(ndgnt)
C
C------------------------------------------------------------------
C     DATA FOR FOUR POINT QUADRATURE
C------------------------------------------------------------------
C      DATA X1/-0.8611363156D0, -0.3399810436D0, 0.3399810436D0,
C    &          0.8611363156D0,  4*0.0D0 /
C      DATA X2/ 0.3225476896D0,  1.7457611012D0, 4.5366202969D0,
C    &          9.3950709123D0,  4*0.0D0 /
C      DATA W1/ 0.3478548451D0,  0.6521451548D0, 0.6521451548D0,
C    &0.3478548451D0,4*0.0D0/
C      DATA W2/ 6.03154104340D-1, 3.5741869244D-1, 3.8887908515D-2,
C    &          5.3929470556D-4 , 4*0.0D0 /
C
C------------------------------------------------------------------
C     DATA FOR EIGHT POINT QUADRATURE
C------------------------------------------------------------------
       DATA X1/-0.9602898564D0, -0.7966664774D0, -0.5255324099D0,
     &         -0.1834346424D0,  0.1834346424D0,  0.5255324099D0,
     &          0.7966664774D0,  0.9602898564D0   /

       DATA X2/ 0.170279632305,  0.9037017768D0,  2.2510866299D0,
     &          4.2667001703D0,  7.0459054024D0, 1.07585160102D1,
     &         1.57406786413D1, 2.28631317369D1  /

       DATA W1/ 0.1012285362D0,  0.2223810344D0, 0.3137066458D0,
     &          0.3626837833D0,  0.3626837833D0, 0.3137066458D0,
     &          0.2223810344D0,  0.1012285362D0   /

       DATA W2/ 3.69188589342D-1, 4.18786780814D-1, 1.75794986637D-1,
     &          3.33434922612D-2, 2.79453623523D-3, 9.07650877336D-5,
     &          8.48574671627D-7, 1.04800117487D-9  /

C------------------------------------------------------------------
C
C  IB=1,B>=0.9; IB=2,B<0.9
C  IB1=1,B1>=0.9; IB1=2,B1<0.9
C  IT=1,TR>=TF; IT=2,TR<TE
C  FORTRAN STATEMENT FUNCTION FOLLOWS FOR CONVENIENCE
       GII(VVE)=AGIIDW(VVE,V,N,L,L1,LP,ISP,LT,LT1,IS,IRESOL)
       M=8
c       ie=0
       IF(B-0.9D0)20,10,10
   10  IB=1
C  INTEGRAL FORM A..........B>=0.9
       F2=1.0D0
       T2=GII(0.0D0)
         ie=ie+1
	 gaunt(ie)=t2
	 energy(ie)=0.0D0 
       U2=EEI(B)
       PREC=T2*U2
       GO TO 30
   20  IB=2
C  INTEGRAL FORM B..........B<0.9
       EX=DEXP(B)
       F1=-0.5D0*DLOG(B)*EX
       F2=3.678794117D-1*EX
       T2=GII(1.0D0/B-1.0D0)
         ie=ie+1
	 gaunt(ie)=t2
	 energy(ie)=1.0D0/B-1.0D0
       PREC=F2*T2*0.596347D0
   30  IF(B1-0.9D0)50,40,40
   40  IB1=1
C  INTEGRAL FORM C..........B1>=0.9
       F4=1.0D0
       T4=GII(0.0D0)
         ie=ie+1
	 gaunt(ie)=t4
	 energy(ie)=0.0D0
       PION=T4*EEI(B1)
       IF(B1-B)42,42,46
   42  IT=1
C  INTEGRAL FORM E..........B>=0.9,B1>=0.9,B1<=B
   43  F6=1.0D0
       T6=T4
       U6=1.0D0+B/B1
       PSTIM=T6*EEI(B+B1)
       GO TO 60
   46  IT=2
       GO TO (49,47),IB
C  INTEGRAL FORM G..........B<0.9 ,B1>=0.9 ,B1>B
   47  GO TO 43
C  INTEGRAL FORM I..........B>=0.9 , B1>=0.9 ,B1>B
   49  GO TO 43
   50  IB1=2
C  INTEGRAL FORM D..........B1<0.9
       EX1=DEXP(B1)
       F3=(1.0D0-B1)*EX1
       F4=3.678794117D-1*EX1
       T4=GII(1.0D0/B1-1.0D0)
         ie=ie+1
	 gaunt(ie)=t4
	 energy(ie)=1.0D0/B1-1.0D0
       PION=F4*T4*0.596347D0
       IF(B1-B)52,52,56
   52  IT=1
       GO TO (54,53),IB
C INTEGRAL FORM F..........B<0.9 ,B1<0.9,B1<=B
   53  F5=(1.0D0-B)*EX*EX1
       F6=F2*EX1
       T6=GII(1.0D0/B-1.0D0)
         ie=ie+1
	 gaunt(ie)=t6
	 energy(ie)=1.0D0/B-1.0D0
       PSTIM=F6*B*T6*0.403653D0/B1
       GO TO 60
C INTEGRAL FORM J..........B>=0.9 ,B1<0.9 ,B1<=B
   54  F6=EX1
       T6=GII(0.0D0)
         ie=ie+1
	 gaunt(ie)=t6
	 energy(ie)=0.0D0
       PSTIM=F6*T6*(1.0D0-B*U2)/B1
       GO TO 60
   56  IT=2
C  INTEGRAL FORM H..........B<0.9 ,B1<0.9 ,B1>B
       F5=F3*EX
       U6=1.0D0+B/B1
       F6=DEXP(-U6)*EX*EX1
       T6=T2
       PSTIM=F6*T6*EEI(U6)
   60  CONTINUE
C  60  WRITE(I4UNIT(-1),100)PREC,PION,PSTIM,B,B1
       DO 80 I=1,M
       GO TO (62,64),IB
C  INTEGRAL FORM A
   62  Z2=X2(I)+B
       t1=GII(Z2/B-1.0D0)
       PREC=PREC+F2*W2(I)*(t1-T2)/Z2
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z2/B-1.0D0
       GO TO (66,70),IB1
C  INTEGRAL FORM B
   64  Z1=B**(0.5D0*(X1(I)+1.0D0))
       Z2=X2(I)+1.0D0
       t1=GII(Z1/B-1.0D0)
       t3=GII(Z2/B-1.0D0)
       PREC=PREC+F1*W1(I)*DEXP(-Z1)*t1+F2*W2(I)*(
     &t3-T2)/Z2
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z1/B-1.0D0
	 ie=ie+1
	 gaunt(ie)=t3
	 energy(ie)=Z2/B-1.0D0
       GO TO (66,70),IB1
C  INTEGRAL FORM C
   66  Z2=X2(I)+B1
       T=0.0D0
       IF(Z2.LE.193.0D0)T=DEXP(-Z2)
       t1=GII(Z2/B1-1.0D0)
       PION=PION+F4*W2(I)*(t1/(1.0D0-T)-T4)/Z2
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z2/B1-1.0D0
C  INTEGRAL FORM E ,I ,G
   67  Z2=(X2(I)+B+B1)/U6
       T=0.0D0
       IF(Z2.LE.193.0D0)T=DEXP(-Z2)
       t1=GII(Z2/B1-1.0D0)
       PSTIM=PSTIM+F6*W2(I)*(t1/(1.0D0-T)-T6)/(Z2*U6)
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z2/B1-1.0D0
       GO TO 80
C  INTEGRAL FORM D
   70  Z1=2.0D0*B1/((1.0D0-B1)*X1(I)+(1.0D0+B1))
       D1=2.0*B1
       IF(Z1-1.0D-6)82,82,81
   81  D1=D1*(DEXP(Z1)-1.0D0)/Z1
   82  Z2=X2(I)+1.0D0
       t1=GII(Z1/B1-1.0D0)
       t3=GII(Z2/B1-1.0D0)
       PION=PION+F3*W1(I)*t1/D1+F4*W2(I)*(t3
     &/(1.0D0-DEXP(-Z2))-T4)/Z2
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z1/B1-1.0D0
	 ie=ie+1
	 gaunt(ie)=t3
	 energy(ie)=Z2/B1-1.0D0
       GO TO (72,78),IT
   72  GO TO (76,74),IB
C  INTEGRAL FORM F
   74  Z1=2.0D0*B1/((1.0D0-B)*X1(I)+(1.0D0+B))
       D1=2.0D0*B1
       IF(Z1-1.0D-6)84,84,83
   83  D1=D1*(DEXP(Z1)-1.0D0)/Z1
   84  U6=B1/B
       Z2=U6*(X2(I)+1.0D0)
       IF(Z2-1.0D-6)86,86,85
   85  D2=(DEXP(Z2)-1.0D0)
       GO TO 87
   86  D2=Z2
   87  t1=GII(Z1/B1-1.0D0)
       t3=GII(Z2/B1-1.0D0)
       PSTIM=PSTIM+F5*W1(I)*DEXP(-B*Z1/B1)*t1/D1+F6*W2(I)*
     &U6*(t3/D2-T6/Z2)/Z2
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z1/B1-1.0D0
	 ie=ie+1
	 gaunt(ie)=t3
	 energy(ie)=Z2/B1-1.0D0
       GO TO 80
C  INTEGRAL FORM J
   76  U6=B1/B
       Z2=U6*(X2(I)+B)
       IF(Z2-1.0D-6)88,88,89
   88  D2=Z2
       GO TO 90
   89  D2=DEXP(Z2)-1.0D0
   90  t1=GII(Z2/B1-1.0D0)
       PSTIM=PSTIM+F6*W2(I)*U6*(t1/D2-T6/Z2)/Z2
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z2/B1-1.0D0
       GO TO 80
C  INTEGRAL FORM H
   78  Z1=2.0D0*B1/((1.0D0-B1)+X1(I)+(1.0D0+B1))
       D1=2.0D0*B1*DEXP(B*Z1/B1)
       IF(Z1-1.0D-6)92,92,91
   91  D1=D1*(DEXP(Z1)-1.0D0)/Z1
   92  Z2=X2(I)+U6
       t1=GII(Z1/B-1.0D0)
       t3=GII(Z2/(U6*B1)-1.0D0)
       PSTIM=PSTIM+F5*W1(I)*t1/D1+F6+W2(I)*(
     &t3/(1.0D0-DEXP(-Z2/U6))-T6)/Z2
         ie=ie+1
	 gaunt(ie)=t1
	 energy(ie)=Z1/B-1.0D0
	 ie=ie+1
	 gaunt(ie)=t3
	 energy(ie)=Z2/(U6*B1)-1.0D0
       CONTINUE
   80  CONTINUE

  100  FORMAT(1P5D12.4)
       
       RETURN
       END
