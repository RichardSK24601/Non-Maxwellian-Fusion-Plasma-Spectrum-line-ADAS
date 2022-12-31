CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/pchg.for,v 1.3 2007/07/20 09:46:35 allan Exp $ Date $Date: 2007/07/20 09:46:35 $
CX
       FUNCTION PCHG(V,L,L1,E,FACT)
       IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------
C
C  PURPOSE: EVALUATES PEACH AMPLITUDE G(V,L,E,L1) BY TABLE INTERPOLATION
C
C  TABLES ARE READ IN FROM DISC FILE 'JETSHP.PCHGTAB.DATA' ON FIRST
C  CALL TO PCHG AND STORED IN LABELLED COMMON /PCHGTB/.
C  THE INPUT IS ON STREAM 13.  IGONE MUST BE SET TO 1 IN MAIN ROUTINE
C  INITIALLY. IT IS RESET TO 0 AFTER TABLES ARE READ IN.
C  ********** H.P. SUMMERS, JET     26 JUNE 1985  ********************
C_______________________________________________________________________
C  INPUT
C      N=EFFECTIVE PRINCIPAL QUANTUM NUMBER OF BOUND STATE
C      L=ORBITAL ANGULAR MOMENTUM OF ELECTRON IN BOUND STATE
C      L1=ORBITAL ANGULAR MOMENTUM OF ELECTRON IN FREE STATE
C      E=ELECTRON ENERGY IN FREE STATE (REDUCED RYDBERG UNITS)
C  OUTPUT
C      PCHG=PEACH INTERPOLATED AMPLITUDE G OR G*
C      FACT=1                IF PEACH G  RETURNED
C          =DSQRT(DABS(V-L)) IF PEACH G* RETURNED
C                            NB OBTAIN G FROM G* BY MULTIPLYING BY FACT
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
       DIMENSION GA(16,22,6),EA(16),VA(22),J0A(6),JCA(6)
       COMMON /PCHGTB/GA,EA,VA,J0A,JCA,IGONE
       IF(IGONE.LT.1)GO TO 20
       READ(13,*)(EA(I),I=1,16)
       READ(13,*)(VA(J),J=1,22)
       READ(13,*)(J0A(K),K=1,6)
       READ(13,*)(JCA(K),K=1,6)
       DO 10 K=1,6
       J0=J0A(K)
   10  READ(13,*)((GA(I,J,K),I=1,16),J=J0,22)
       IGONE=0
   20  K=2*L+MAX0(0,L1-L)
       IF(K.GE.1.AND.K.LE.6)GO TO 30
   25  WRITE(I4UNIT(-1),100)K,V,L,L1,E
       PCHG=0.0D0
       RETURN
   30  IF(IABS(L1-L).NE.1)GO TO 25
       J0=J0A(K)
       JC=JCA(K)
       IF(V.LT.VA(J0).OR.V.GT.12.0D0.OR.E.LT.0.0D0)GO TO 25
       EL=L
       T=1.0D0
       FACT=T
C  LOCATE E IN VECTOR EA
       I0=1
   35  I0=I0+1
       IF(I0.GE.16)GO TO 37
       IF(E.GT.EA(I0))GO TO 35
       I0=I0-1
       GO TO 38
   37  I0=15
C  LOCATE V IN VECTOR VA
   38  J1=J0
   40  J1=J1+1
       IF(V.GT.VA(J1))GO TO 40
       J1=J1-1
       IF(J1-JC+1)44,42,46
   42  T=1.0D0/DSQRT(DABS(VA(JC)-EL))
   44  FACT=DSQRT(DABS(V-EL))
   46  VVJ=VA(J1)**2
       B=DLOG(GA(I0,J1,K)/GA(I0+1,J1,K))/DLOG((1.0D0+VVJ*EA(I0))/(1.0D0
     &+VVJ*EA(I0+1)))
       GP1=GA(I0,J1,K)*((1.0D0+VVJ*E)/(1.0D0+VVJ*EA(I0)))**B
       VVJ=VA(J1+1)**2
       B=DLOG(GA(I0,J1+1,K)/GA(I0+1,J1+1,K))/DLOG((1.0D0+VVJ*EA(I0))/
     &(1.0D0+VVJ*EA(I0+1)))
       GP2=GA(I0,J1+1,K)*T*((1.0D0+VVJ*E)/(1.0D0+VVJ*EA(I0)))**B
       U1=(V-VA(J1+1))/(VA(J1)-VA(J1+1))
       U2=(V-VA(J1))/(VA(J1+1)-VA(J1))
       PCHG=U1*GP1+U2*GP2
       RETURN
  100  FORMAT(1H0,'FAIL IN PCHG: V, L, L1 OR E OUT OF RANGE',2X,'K=',I2,
     &2X,'V=',F5.2,2X,'L=',I2,2X,'L1=',I2,2X,'E=',F8.2)
      END
