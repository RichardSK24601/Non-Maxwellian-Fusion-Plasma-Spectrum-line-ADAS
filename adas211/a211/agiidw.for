CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/agiidw.for,v 1.3 2007/05/16 10:16:49 allan Exp $ Date $Date: 2007/05/16 10:16:49 $
CX
       FUNCTION AGIIDW(VVE,V,N,L,L1,LP,ISP,LT,LT1,IS,IRESOL)
       IMPLICIT REAL*8(A-H,O-Z)
C ----------------------------------------------------------------------
C
C  VERSION OF GIIDW FOR USE BY APHOTDW CALLED BY ADASRRC
C  ******************** H.P.SUMMERS, JET   30 JUNE 1992 ***************
C
C  PURPOSE: CALCULATES BOUND-FREE G-FACTORS USING DISTORTED WAVES,
C  BURGESS-SEATON PEACH OR HYDROGENIC APPROXIMATIONS
C
C  ******************** H.P.SUMMERS, JET   19 AUG. 1984 ***************
C  MAY SELECT DISTORTED WAVE MATRIX ELEMENTS, FROM PEACH TABLES, FROM
C  ORIGINAL BURGESS-SEATON APPROXIMATION OR HYDROGENIC MATRIX ELEMENTS
C  USING SELECTOR IBSOPT IN THE /BSPARS/ COMMON BLOCK.
C  FOR COMPLETENESS, THE UNRESOLVED, BUNDLED N, GBF (BURGESS AND SUMMERS
C  ,1976) CAN ALSO BE OBTAINED.
C  THE DRIVING PROGRAM MUST SET COMMON BLOCKS /PCHGTB/ AND /PCHXTB/ FOR
C  USE BY FUNCTIONS PCHG AND PCHX AND SET IFIRST=IGONE=1 AT START UP.
C  /PCHGTB/ DATA IS REQUIRED FROM FILE PCHGTAB.DATA ON STREAM 13
C  /PCHXTB/ DATA IS REQUIRED FROM FILE PCHXTAB.DATA ON STREAM 14
C  INPUT
C      VVE=V**2*E WHERE E=(FREE ELECTRON ENERGY)/Z**2 (RYD)
C      V=EFFECTIVE PRINCIPAL QUANTUM NUMBER OF BOUND ELECTRON
C      N=PRINCIPAL QUANTUM NUMBER OF BOUND ELECTRON
C      L=ORBITAL QUANTUM NUMBER OF BOUND ELECTRON
C      L1=ORBITAL QUANTUM NUMBER OF FREE ELECTRON
C      ISP=2*SP+1 WHERE SP IS TOTAL SPIN OF PARENT STATE
C      LP=TOTAL ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER OF PARENT STATE
C      LT=TOTAL ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER OF BOUND SYSTEM
C      LT1=TOTAL ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER OF FREE SYSTEM
C      IS=2*S+1 WHERE S IS TOTAL SPIN OF SYSTEM
C      IRESOL=1 GIVES GII((LP,SP)N L LT S,(LP,SP)E L1 LT1 S)
C            =2 GIVES GII((LP,SP)N L LT S,(LP,SP)E L1 S) =ABOVE LT1 SUM
C            =3 GIVES GII((LP,SP)N L S,(LP,SP)E L1 S)   = ABOVE LT SUM
C            =4 GIVES GII((LP,SP)N L,(LP,SP)E L1)      = ABOVE S SUM
C            =5 GIVES GII(N,E) = GBF  (BURGESS AND SUMMERS)
C  OUTPUT
C      AGIIDW  THE BOUND-FREE GAUNT FACTOR
C-----------------------------------------------------------------------
C  UPDATE:  01/10/96  HP SUMMERS - BYPASS PEACH DATA INPUT IF IBSOPT=3
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
C VERSION: 1.2                          DATE: 14-10-96
C MODIFIED: WILLIAM OSBORN
C               - ADDED CHANGES DATED 01/10/96 ABOVE.
C
C VERSION: 1.3                          DATE: 16-05-07
C MODIFIED: Allan Whiteford
C               - Modified comments as part of subroutine documentation
C                 procedure.
C
C-----------------------------------------------------------------------
C ----------------------------------------------------------------------
       DIMENSION A(6),B(6),C(6),ALF(6),BET(6)
       DIMENSION G(72),GAM(72)
       DIMENSION G01(8),GAM01(8),X01(8),V01(8)
       DIMENSION RG10(11),G10(11),GAM10(11),X10(11),V10(11)
       DIMENSION RG12(11),G12(11),GAM12(11),X12(11),B12(11),V12(11)
       DIMENSION NLQS(10),NAA(2),LAA(2),EAA(2),QDAA(2),ALFAA(2,10)
       COMMON /DWPARS/Z0,EAA,QDAA,ALFAA,ACC,XMAX,H,NLQS,NSHELL,NAA,LAA,
     &JSN,JEALFA,IONCE
       COMMON /BSPARS/Z,U1,U2,U3,ZETA,IBSOPT,IWARN
C      Z=BOUND STATE ION CHARGE+1
C      U1,U2,U3=QUANTUM DEFECT EXPANSION PARAMETERS FOR FREE STATE
C      ZETA=ZETA PARAMETER FOR BOUND STATE
C      IBSOPT=1 TO USE FITTED PEACH PHASE IN MATRIX ELEMENT
C            =2 TO USE BURGESS-SEATON PHASE IN MATRIX ELEMENT
C            =3 TO USE HYDROGENIC MATRIX ELEMENT
C            =4 TO USE DISTORTED WAVE MATRIX ELEMENT
C      IWARN=0 NO SENSITIVITY (OR HYDROGENIC OPTION SELECTED)
C           =I1+2*I2+4*I3+8*I4   WHERE
C                I1=1 IMPLIES PEACH PHASE SENSITIVITY
C                I2=1 IMPLIES B&S PHASE SENSITIVITY
C                I3=1 IMPLIES IBSOPT CHOICE SENSITIVE
C                I4=1 IMPLIES B&S AND PEACH STRADDLE PI/2
       DATA A/-0.147D0,-0.216D0,-0.120D0,-0.247D0,-0.117D0,-0.362D0/
       DATA B/0.2515D0,-0.171D0,0.600D0,-0.272D0,1.170D0,0.599D0/
       DATA C/-0.078D0,0.0D0,0.0D0,0.0D0,0.0D0,-2.432D0/
       DATA ALF/0.310D0,0.0D0,0.362D0,-0.010D0,0.321D0,-0.390D0/
       DATA BET/0.0D0,0.0D0,0.0535D0,-0.019D0,0.106D0,0.050D0/
       DATA G/2.723D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &2.095D0,1.028D0,2.840D0,0.0D0,0.0D0,0.0D0,
     &1.856D0,1.117D0,2.264D0,0.669D0,3.000D0,0.000D0,
     &1.718D0,1.152D0,2.010D0,0.818D0,2.413D0,0.468D0,
     &1.623D0,1.168D0,1.856D0,0.899D0,2.139D0,0.599D0,
     &1.553D0,1.175D0,1.749D0,0.952D0,1.971D0,0.704D0,
     &1.498D0,1.177D0,1.666D0,0.988D0,1.854D0,0.793D0,
     &1.452D0,1.176D0,1.601D0,1.014D0,1.765D0,0.868D0,
     &1.414D0,1.173D0,1.546D0,1.033D0,1.694D0,0.933D0,
     &1.381D0,1.170D0,1.501D0,1.047D0,1.635D0,0.911D0,
     &1.352D0,1.165D0,1.461D0,1.058D0,1.585D0,1.041D0,
     &1.327D0,1.161D0,1.427D0,1.065D0,1.543D0,1.085D0/
       DATA GAM/1.754D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &1.605D0,1.667D0,1.574D0,0.0D0,0.0D0,0.0D0,
     &1.591D0,1.667D0,1.582D0,1.819D0,1.447D0,0.000D0,
     &1.590D0,1.667D0,1.579D0,1.771D0,1.535D0,1.850D0,
     &1.591D0,1.667D0,1.582D0,1.741D0,1.544D0,1.908D0,
     &1.594D0,1.667D0,1.587D0,1.722D0,1.549D0,1.918D0,
     &1.596D0,1.667D0,1.593D0,1.707D0,1.556D0,1.920D0,
     &1.599D0,1.667D0,1.598D0,1.697D0,1.564D0,1.921D0,
     &1.601D0,1.667D0,1.603D0,1.688D0,1.573D0,1.922D0,
     &1.603D0,1.667D0,1.608D0,1.682D0,1.581D0,1.924D0,
     &1.605D0,1.667D0,1.614D0,1.676D0,1.589D0,1.926D0,
     &1.607D0,1.667D0,1.618D0,1.672D0,1.596D0,1.928D0/
       DATA V01/0.6D0,0.8D0,1.0D0,1.2D0,1.4D0,1.6D0,1.8D0,2.0D0/
       DATA G01/3.259D0,2.976D0,2.739D0,2.527D0,2.360D0,2.244D0,
     &2.162D0,2.095D0/
       DATA GAM01/1.85D0,1.77D0,1.701D0,1.655D0,1.632D0,1.620D0,
     &1.612D0,1.604D0/
       DATA X01/0.143D0,0.085D0,0.043D0,0.011D0,-0.008D0,-0.020D0,
     &-0.031D0,-0.041D0/
       DATA V10/1.0D0,1.2D0,1.4D0,1.6D0,1.8D0,2.0D0,2.2D0,2.4D0,
     &2.6D0,2.8D0,3.0D0/
       DATA RG10/1.88D0,1.50D0,1.31D0,1.18D0,0.0D0,0.0D0,0.0D0,0.0D0,
     &0.0D0,0.0D0,0.0D0/
       DATA G10/0.0D0,0.670D0,0.826D0,0.911D0,0.962D0,0.999D0,
     &1.029D0,1.058D0,1.080D0,1.100D0,1.117D0/
       DATA GAM10/1.333D0,1.515D0,1.585D0,1.630D0,1.655D0,1.667D0,
     &1.667D0,1.667D0,1.667D0,1.667D0,1.667D0/
       DATA X10/-0.330D0,-0.321D0,-0.313D0,-0.306D0,-0.300D0,-0.295D0,
     &-0.290D0,-0.286D0,-0.281D0,-0.277D0,-0.273D0/
       DATA V12/1.0D0,1.2D0,1.4D0,1.6D0,1.8D0,2.0D0,2.2D0,2.4D0,
     &2.6D0,2.8D0,3.0D0/
       DATA RG12/5.69D0,5.57D0,5.02D0,4.25D0,0.0D0,0.0D0,0.0D0,
     &0.0D0,0.0D0,0.0D0,0.0D0/
       DATA G12/0.0D0,2.489D0,3.174D0,3.291D0,3.075D0,2.757D0,2.512D0,
     &2.415D0,2.386D0,2.340D0,2.251D0/
       DATA GAM12/2.340D0,1.911D0,1.703D0,1.625D0,1.624D0,1.658D0,
     &1.675D0,1.635D0,1.593D0,1.576D0,1.597D0/
       DATA X12/0.650D0,0.511D0,0.389D0,0.287D0,0.210D0,0.164D0,
     &0.1425D0,0.131D0,0.115D0,0.0945D0,0.073D0/
       DATA B12/0.079D0,0.069D0,0.054D0,0.038D0,0.029D0,0.035D0,
     &0.053D0,0.068D0,0.068D0,0.060D0,0.050D0/
       AGIIDW=0.0D0
       IWARN=0
       XN=N
       XL=L
       XL1=L1
       XLP=LP
       XSP=ISP
       IF(IRESOL.LE.0.OR.IRESOL.GT.5)RETURN
       GO TO (55,60,65,70,52),IRESOL
   52  AGIIDW=GBF(V,VVE)
       RETURN
   55  XLT1=LT1
   60  XLT=LT
   65  XS=IS
   70  CONTINUE
       GO TO (75,80,85,90),IRESOL
   75  W=WIG6J(XL,XL1,1.0D0,XLT1,XLT,XLP)
       ANG=(2.0*XLT+1.0)*(2.0*XLT1+1.0)*XS*DMAX1(XL,XL1)*W*W/(XSP*(2.0*X
     &LP+1.0))
       CL1=XSP*(2.0*XLP+1.0)*ANG/((2.0*XLT+1.0)*XS)
       GO TO 95
   80  ANG=(2.0*XLT+1.0)*XS*DMAX1(XL,XL1)/((2.0*XL+1.0)*XSP*(2.0*XLP+1.0
     &))
       CL1=XSP*(2.0*XLP+1.0)*ANG/((2.0*XLT+1.0)*XS)
       GO TO 95
   85  ANG=XS*DMAX1(XL,XL1)/XSP
       CL1=XSP*ANG/((2.0*XL+1.0)*XS)
       GO TO 95
   90  ANG=2.0*DMAX1(XL,XL1)
       CL1=ANG/(2.0*(2.0*XL+1.0))
   95  CONTINUE
       E=VVE/(V*V)
       RV=1.0D0/V
       RV2=RV*RV
C  SCALE CROSS-SECTION FOR E>ECUT USING HYDROGENIC FACTORS
C  ***  SET ECUT=3, ORIGINAL SETTING = 1  *******
       ECUT=3.0D0
       FAC=1.0
       EA=E
       IF(EA.GT.ECUT)GO TO 5
       GO TO 7
    5  E=ECUT
       E2=-E
       EA2=-EA
C      WRITE(I4UNIT(-1),*)'N=',N,'  L=',L,'  L1=',L1,'  E2=',E2,
C    & '  RD2FS=',RD2FS(N,L,L1,E2)
       FAC=(1.0D0+V*V*EA)*RD2FS(N,L,L1,EA2)/((1.0D0+V*V*ECUT)*
     &RD2FS(N,L,L1,E2))
    7  CONTINUE
       U=U1+E*(U2+E*U3)
       EAA(2)=E*Z*Z
       IF(IONCE.LE.0)QDAA(2)=U
       XL=L
       XL1=L1
       K=2*L+MAX0(0,L1-L)
       IF(IBSOPT.EQ.4.AND.K.LE.6)GO TO 23
C ******** FORCE HYDROGENIC OPTION IF K.EQ.6. AND V.LT.4.0D0 **********
C ******** FORCE HYDROGENIC OPTION IF K.EQ.5. AND V.LT.3.0D0 **********
C ******** FORCE HYDRODENIC OPTION IF K.EQ.4. AND V.LT.3.0D0 **********
C ******** ALTHOUGH ALLOWED BY PCHG AND PCHX                   ********
       IF(K.EQ.6.AND.V.LT.4.0D0)GO TO 40
       IF(K.EQ.5.AND.V.LT.3.0D0)GO TO 40
       IF(K.EQ.4.AND.V.LT.3.0D0)GO TO 40
C----------------------------------------------------------------------
C  ADD EXTRA BYPASS IF IBSOPT = 3                  01/10/96  HP SUMMERS
C----------------------------------------------------------------------
       IF(IBSOPT.EQ.3)GO TO 40
C
       IF(K-6)10,10,40
   10  IF(V-12.0D0)12,12,40
   12  IF(V.LT.XL+2.0D0.AND.L.LE.1)GO TO 25
   15  N0=V
       N0=N0-1
       IF(N0.LT.L+1)N0=L+1
       IF(N0.GE.11)N0=10
       XLL=A(K)+RV*B(K)+RV2*C(K)
       X=XLL+E*ALF(K)/(RV+E)+E*BET(K)/(RV2+E)
C      WRITE(I4UNIT(-1),101)X,XLL,ALF(K),BET(K)
       V0=N0
       T0=0.5D0*(V-V0-1.0D0)*(V-V0-2.0D0)
       T1=-(V-V0)*(V-V0-2.0D0)
       T2=0.5D0*(V-V0)*(V-V0-1.0D0)
       K1=6*N0+K
       GLL=T0*G(K1-6)+T1*G(K1)+T2*G(K1+6)
       GAMLL=T0*GAM(K1-6)+T1*GAM(K1)+T2*GAM(K1+6)
C      WRITE(I4UNIT(-1),100)N0,K,K1
C      WRITE(I4UNIT(-1),101)V,V0,T0,T1,T2,GAM(K1)
   16  PB=(-1.0D0)**(L+1)*GLL*(1.0D0+E*V*V)**(-GAMLL)
       PP=PCHG(V,L,L1,E,FACT)
       ISWIT=1
C  X IS BURGESS-SEATON PHASE, XX IS PEACH PHASE
       XX=PCHX(V,L,L1,E)
       BSARG=DMOD(V+U+X,1.0D0)
       PARG=DMOD(V+U+XX,1.0D0)
       GO TO (61,62,64,64),IBSOPT
   61  ARG=PARG
       GO TO 63
   62  ARG=BSARG
   63  IF(DABS(PARG-0.5D0).LT.0.1D0)IWARN=IWARN+1
       IF(DABS(BSARG-0.5D0).LT.0.1D0)IWARN=IWARN+2
       IF(DABS(ARG-0.5D0).LT.0.1D0)IWARN=IWARN+4
       IF((BSARG-0.5D0)*(PARG-0.5D0).LT.0.0D0)IWARN=IWARN+8
C      WRITE(I4UNIT(-1),105)V,E,L,L1,ZETA,U,PB,PP,X,XX,FAC,FACT
   64  GO TO (43,44,43,23),IBSOPT
   43  P=0.79788D0*PP*DSQRT(V)/(1.0D0+E*V*V)**2
       P=P*DCOS(3.14159D0*(V+U+XX))
       IF(FACT.LE.1.0D-3) GO TO 45
       P=FACT*P/DSQRT(ZETA)
       GO TO 45
   44  P=PB*DCOS(3.14159D0*(V+U+X))/DSQRT(ZETA)
C      WRITE(I4UNIT(-1),101)GLL,GAMLL,X,PB,PP,P,FAC,XX
   45  P2=P*P
   17  THET=(1.0D0+V*V*E)*P2*FAC
       Q=8.559D-19*CL1*THET*V*V/Z/Z
       GO TO (18,19),ISWIT
   18  E2=-E
       THETH=FAC*(1.0D0+V*V*E)*RD2FS(N,L,L1,E2)/(V**4)
       QH=8.559D-19*CL1*THETH*V*V/Z/Z
       ISWIT=1
       GO TO 20
   19  THETH=THET
       ISWIT=2
   20  E=EA
C      WRITE(I4UNIT(-1),103)THET,THETH,ANG,CL1
       GO TO (21,21,22),IBSOPT
   21  AGIIDW=5.412659D-2*(1.0D0+V*V*E)**3*ANG*THET/V
C      WRITE(I4UNIT(-1),104)IBSOPT,V,E,EA,ANG,THET,AGIIDW
       RETURN
   22  AGIIDW=5.412659D-2*(1.0D0+V*V*E)**3*ANG*THETH/V
C      WRITE(I4UNIT(-1),104)IBSOPT,V,E,EA,ANG,THETH,AGIIDW
       RETURN
   23  THETDW=FAC*(1.0D0+V*V*E)*ADWRD2(N,L,L1)*(Z/V)**4
C      WRITE(I4UNIT(-1),106)
  106  FORMAT(1H0,'ADWRD2 CALLED')
       E=EA
       AGIIDW=5.412659D-2*(1.0D0+V*V*E)**3*ANG*THETDW/V
C      WRITE(I4UNIT(-1),104)IBSOPT,V,E,EA,ANG,THETDW,AGIIDW
       RETURN
   25  I0=5.0D0*V
       GO TO (26,27,27),K
   26  I0=I0-2
       IF(I0.LT.1)I0=1
       IF(I0.GT.6)I0=6
       V0=I0+2
       GO TO 28
   27  I0=I0-4
       IF(I0.LT.1)I0=1
       IF(I0.GT.9)I0=9
       V0=I0+4
   28  V0=0.2D0*V0
       T0=12.5D0*(V-V0-0.2D0)*(V-V0-0.4D0)
       T1=-25.0D0*(V-V0)*(V-V0-0.4D0)
       T2=12.5D0*(V-V0)*(V-V0-0.2D0)
       GO TO (30,32,36),K
   30  XLL=T0*X01(I0)+T1*X01(I0+1)+T2*X01(I0+2)
       X=XLL+0.310D0*E/(RV+E)
       GLL=T0*G01(I0)+T1*G01(I0+1)+T2*G01(I0+2)
       GAMLL=T0*GAM01(I0)+T1*GAM01(I0+1)+T2*GAM01(I0+2)
C      WRITE(I4UNIT(-1),101)V,V0,T0,T1,T2,XLL,X,GLL,GAMLL
       GO TO 16
   32  XLL=T0*X10(I0)+T1*X10(I0+1)+T2*X10(I0+2)
       X=XLL
       GAMLL=T0*GAM10(I0)+T1*GAM10(I0+1)+T2*GAM10(I0+2)
       IF(I0-2)33,33,34
   33  RGLL=T0*RG10(I0)+T1*RG10(I0+1)+T2*RG10(I0+2)
       GLL=DSQRT(V-1.0D0)*RGLL
       GO TO 16
   34  GLL=T0*G10(I0)+T1*G10(I0+1)+T2*G10(I0+2)
       GO TO 16
   36  XLL=T0*X12(I0)+T1*X12(I0+1)+T2*X12(I0+2)
       BLL=T0*B12(I0)+T1*B12(I0+1)+T2*B12(I0+2)
       X=XLL+0.362D0*E/(RV+E)+BLL*E/(RV2+E)
       GAMLL=T0*GAM12(I0)+T1*GAM12(I0+1)+T2*GAM12(I0+2)
       IF(I0-2)37,37,38
   37  RGLL=T0*RG12(I0)+T1*RG12(I0+1)+T2*RG12(I0+2)
       GLL=DSQRT(V-1.0D0)*RGLL
       GO TO 16
   38  GLL=T0*G12(I0)+T1*G12(I0+1)+T2*G12(I0+2)
       GO TO 16
   40  E2=-E
       P2=RD2FS(N,L,L1,E2)/(V**4)
       ISWIT=2
       GO TO 17
  100  FORMAT(4I5,3F10.5)
  101  FORMAT(1P7D12.4)
  102  FORMAT(I5,3F10.5)
  103  FORMAT(1H ,' THET=',1PD12.4,' THETH=',1PD12.4,' ANG=',1PD12.4,' C
     &L1=',1PD12.4)
  104  FORMAT(1H ,'    IBSOPT=',I2,3X,'V=',F10.5,3X,'E=',F10.5,3X,'EA=',
     &F10.5,3X,'ANG=',1PE12.4,3X,'THET=',1PE12.4,3X,'GIIDW=',1PE12.4)
  105  FORMAT(1H ,'    V=',F10.5,3X,'E=',F10.5,3X,'L=',I2,3X,'L1=',I2,
     &3X,'ZETA=',1PD12.4,3X,'U=',1PD12.4/1H ,'    PB=',1PE12.4,3X,
     &'PP=',1PD12.4,3X,'X=',1PD12.4,3X,'XX=',1PD12.4,3X,'FAC=',1PD12.4,
     &3X,'FACT=',1PD12.4)
      END
