CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/pchx.for,v 1.3 2007/07/20 09:46:35 allan Exp $ Date $Date: 2007/07/20 09:46:35 $
CX
       FUNCTION PCHX(V,L,L1,E)
       IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C  PURPOSE: EVALUATES PEACH PHASE X(V,L,E,L1) BY TABLE INTERPOLATION
C
C  TABLES ARE READ IN FROM DISC FILE 'JETSHP.PCHXTAB.DATA' ON FIRST
C  CALL TO PCHX AND STORED IN LABELLED COMMON /PCHXTB/.
C  THE INPUT IS ON STREAM 14.  IFIRST MUST BE SET TO 1 IN MAIN ROUTINE
C  INITIALLY. IT IS RESET TO 0 AFTER TABLES ARE READ IN.
C  ********** H.P. SUMMERS, JET     26 JUNE 1985  ********************
C_______________________________________________________________________
C  INPUT
C      N=EFFECTIVE PRINCIPAL QUANTUM NUMBER OF BOUND STATE
C      L=ORBITAL ANGULAR MOMENTUM OF ELECTRON IN BOUND STATE
C      L1=ORBITAL ANGULAR MOMENTUM OF ELECTRON IN FREE STATE
C      E=ELECTRON ENERGY IN FREE STATE (REDUCED RYDBERG UNITS)
C  OUTPUT
C      PCHX=PEACH INTERPOLATED PHASE
C_______________________________________________________________________
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
C VERSION: 1.3                          DATE: 20-07-07
C MODIFIED: Allan Whiteford
C               - Small modification to comments to allow for
C                 automatic documentation preparation.
C
C-----------------------------------------------------------------------
       DIMENSION XA(17,23,6),REA(17),RVA(23),J0A(6)
       COMMON /PCHXTB/XA,REA,RVA,J0A,IFIRST
       IF(IFIRST.LT.1)GO TO 20
       READ(14,*)(REA(I),I=1,17)
       READ(14,*)(RVA(J),J=1,23)
       READ(14,*)(J0A(K),K=1,6)
       DO 10 K=1,6
       J0=J0A(K)
   10  READ(14,*)((XA(I,J,K),I=1,17),J=J0,23)
       IFIRST=0
   20  K=2*L+MAX0(0,L1-L)
       IF(K.GE.1.AND.K.LE.6)GO TO 30
   25  WRITE(I4UNIT(-1),100)
       PCHX=0.0D0
       RETURN
   30  IF(IABS(L1-L).NE.1)GO TO 25
       J0=J0A(K)
       RV=1.0D0/V
       IF(RV.GT.RVA(J0).OR.E.LT.0.0D0)GO TO 25
       RE=1.0D0/(E+0.02D0)
C  LOCATE RE IN VECTOR REA
       I0=1
   35  I0=I0+1
       IF(RE.LT.REA(I0))GO TO 35
       I0=I0-1
C  LOCATE RV IN VECTOR RVA
       J1=J0
   40  J1=J1+1
       IF(RV.LT.RVA(J1))GO TO 40
       J1=J1-1
       T1=(RE-REA(I0+1))/(REA(I0)-REA(I0+1))
       T2=(RE-REA(I0))/(REA(I0+1)-REA(I0))
       U1=(RV-RVA(J1+1))/(RVA(J1)-RVA(J1+1))
       U2=(RV-RVA(J1))/(RVA(J1+1)-RVA(J1))
       PCHX=T1*U1*XA(I0,J1,K)+T2*U1*XA(I0+1,J1,K)+T1*U2*XA(I0,J1+1,K)+
     &T2*U2*XA(I0+1,J1+1,K)
       RETURN
  100  FORMAT(1H0,'FAIL IN PCHX: V, L, L1 OR E OUT OF RANGE')
      END
