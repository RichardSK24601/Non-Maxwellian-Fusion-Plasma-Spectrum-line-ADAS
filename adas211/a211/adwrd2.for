CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/adwrd2.for,v 1.3 2017/08/21 19:52:14 mog Exp $ Date $Date: 2017/08/21 19:52:14 $
CX
       FUNCTION ADWRD2(N,L,L1)
       IMPLICIT REAL*8(A-H,O-Z)
C ---------------------------------------------------------------------
C
C  VERSION OF DWRD2 FOR USE BY ADASRRC WHICH MAKES USE OF ADWLPOL
C  *********  H.P. SUMMERS, JET    30 JUNE 1992 ***********************
C
C  PURPOSE: CALCULATES SQUARE OF BOUND FREE DIPOLE INTEGRAL IN
C  DISTORTED WAVE APPROXIMATION.
C
C  THIS FUNCTION ACTS AS AN INTERFACE BETWEEN GIIDW AND DWLPOL. 
C  ATOMIC STRUCTURE AND POTENTIAL DATA IS BROUGHT IN VIA LABELLED
C  COMMON BLOCK  /DWPARS/
C  *********  H.P. SUMMERS, JET 19 AUGUST 1985  ***********************
C  INPUT
C      N=PRINCIPAL QUANTUM NUMBER OF BOUND ELECTRON
C      L=ORBITAL QUANTUM NUMBER OF BOUND ELECTRON
C      L1=ORBITAL QUANTUM NUMBER OF FREE ELECTRON
C  OUTPUT
C      ADWRD2=SQUARED RADIAL DIPOLE INTEGRAL
C ---------------------------------------------------------------------
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
C VERSION: 1.2                          DATE: 16-05-07
C MODIFIED: Allan Whiteford
C               - Modified comments as part of subroutine documentation
C                 procedure.
C
C VERSION : 1.3                          
C DATE    : 21-08-2017
C MODIFIED: Martin O'Mullane
C               - Add (internal only) open17 logical to adwlpol to
C                 make the number of arguments match. 
C-----------------------------------------------------------------------
       DIMENSION NLQS(10),NA(2),LA(2),EA(2),QDA(2),ALFAA(2,10)
       COMMON /DWPARS/Z0,EA,QDA,ALFAA,ACC,XMAX,H,NLQS,NSHELL,NA,LA,JSN,
     &JEALFA,IONCE
       LOGICAL OPEN17
       open17=.false.
       IREPT=0
       IEXT=0
       LAM=1
       IF(IONCE.LE.0)GO TO 20
       IF(NA(1).NE.N.OR.LA(1).NE.L.OR.LA(2).NE.L1)GO TO 20
       IREPT=2
       GO TO 21
   20  IONCE=0
   21  CALL ADWLPOL(Z0,NLQS,NSHELL,NA,LA,EA,QDA,ALFAA,JSN,JEALFA,ACC,
     &XMAX,H,LAM,IREPT,IEXT,ANS,open17)
       IONCE=1
       ADWRD2=ANS*ANS
       RETURN
      END
