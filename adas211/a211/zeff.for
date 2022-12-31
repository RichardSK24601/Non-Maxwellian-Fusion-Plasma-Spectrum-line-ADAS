CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/zeff.for,v 1.2 2004/07/06 15:41:11 whitefor Exp $ Date $Date: 2004/07/06 15:41:11 $
CX
      REAL*8 FUNCTION ZEFF(I,Z0,NSHELL,N,NUMEL,ALFA,R)

C-----------------------------------------------------------------------
C      Z0=NUCLEAR CHARGE (IN ELECTRON CHARGE UNITS)
C      NSHELL=NUMBER OF ELECTRON SHELLS
C      N(J),J=1, NSHELL, IS THE PRINCIPAL QUANTUM NUMBER OF SHELL J
C      NUMEL(J) IS THE NUMBER OF ELECTRONS IN SHELL J
C      R IS THE RADIAL DISTANCE COORDINATE IN UNITS OF A0 AND ALFA(J)
C      IS A RADIAL SCALING FACTOR FOR SHELL J.
C     IF I IS SET NEGATIVE, ZEFF IS SET TO THE JUCYS FIT TO THE
C      THOMAS-FERMI EFFECTIVE CHARGE.
C     IF I IS SET POSITIVE OR ZERO, SLATER-TYPE ORBITALS ARE USED, AND
C      ZEFF IS THE EFFECTIVE CHARGE, AS SEEN BY A (SPECTATOR) ELECTRON
C      IN SHELL I, DUE TO THE NUCLEUS PLUS ALL THE OTHER ELECTRONS IN
C      SHELLS J=1,NSHELL.
C      N.B. I MAY BE SET GREATER THAN NSHELL (OR ZERO), IN WHICH
C      CASE THERE IS SCREENING BY ALL THE ELECTRONS IN SHELLS J=1,NSHELL
C      SINCE THE SPECTATOR ELECTRON IS NOT ONE OF THEM.
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
       DIMENSION N(10),NUMEL(10),ALFA(10)
       T=-Z0
       IF(NSHELL)5,5,6
    6  IF(I)7,8,8
    7  T2=(-Z0)**(1.0D0/3.0D0)
       DO 9 J=1,NSHELL
       T1=NUMEL(J)
       RHO=ALFA(J)*R*T2
       T3=0.2075D0*RHO
       T4=0.0D0
       IF(T3-140.0D0)12,13,13
   12  T4=DEXP(-T3)
   13  T=T-T1+T1*T4/(1.0D0+1.19D0*RHO)
    9  CONTINUE
       GO TO 5
    8  Z=-Z0
       DO 4 J=1,NSHELL
       X=ALFA(J)
       X=X*R
       T1=NUMEL(J)-1
       Z=Z-0.5D0*T1
       EN=N(J)
       RHO=2.0D0*Z*X/EN
       IMAX=N(J)+N(J)-1
       T2=EN+EN
       T3=1.0D0
       T4=1.0D0/T2
       DO 1 I1=1,IMAX
       T5=I1
       T2=T2-1.0D0
       T4=T4*RHO/T5
       T3=T3+T4*T2
   1   CONTINUE
       IF(I-J)2,3,2
   2   T1=NUMEL(J)
    3  T4=0.0D0
       IF(RHO-140.0D0)10,11,11
   10  T4=DEXP(-RHO)
   11  T=T-T1*(1.0D0-T3*T4)
       T1=NUMEL(J)+1
       Z=Z-0.5D0*T1
    4  CONTINUE
    5  ZEFF=-T
       RETURN
      END
