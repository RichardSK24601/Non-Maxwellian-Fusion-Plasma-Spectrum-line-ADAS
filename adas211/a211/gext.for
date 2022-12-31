CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/gext.for,v 1.2 2007/07/20 09:46:35 allan Exp $ Date $Date: 2007/07/20 09:46:35 $
CX
       FUNCTION GEXT(X,N,L)
       IMPLICIT REAL*8(A-H,O-Z)
C
C  PURPOSE: PRODUCES ONE ELECTRON ORBITALS FROM SPECIFIED
C           FUNCTIONAL FORMS
C
C  FOR USE IN DWBES, RDWBES,DWDIP WITH EXTERNAL OPTION IEXT=1
C  ____________________________________________________________________
C  HIBBERT (CIV3 PROGRAM) ORBITALS FOR OII 24/4/85
C  ____________________________________________________________________
C  INDEXING OF WAVE FUNCTIONS BY I=(N*(N-1))/2+L+1
C UNIX-IDL PORT:
C
C-----------------------------------------------------------------------
C AUTHOR:  WILLIAM OSBORN (TESSELLA SUPPORT SERVICES PLC)
C
C DATE:    4TH JULY 1996
C
C VERSION: 1.1                          DATE: 04-07-96
C MODIFIED: WILLIAM OSBORN
C               - FIRST VERSION.
C
C VERSION: 1.2                          DATE: 20-07-07
C MODIFIED: Allan Whiteford
C               - Small modification to comments to allow for
C                 automatic documentation preparation.
C-----------------------------------------------------------------------
C
       I=(N*(N-1))/2+L+1
       GO TO (1,2,3,4,5,6),I
    1  GEXT=X*(38.1978304D0*DEXP(-7.4780300*X)+4.9817906D0*DEXP(-12.6307
     &000D0*X))+X*X*(0.0928714D0*DEXP(-3.1009000D0*X)+2.1368144D0*DEXP(-
     &6.3727700D0*X)-0.0087901D0*DEXP(-2.0732300D0*X))
       RETURN
    2  GEXT=X*(-9.6934267D0*DEXP(-7.4780300D0*X)-0.5036558D0*DEXP(-12.63
     &07000D0*X))+X*X*(9.2494101D0*DEXP(-3.1009000D0*X)-11.2771775D0*DEX
     &P(-6.3727700D0*X)+4.5358980D0*DEXP(-2.0732300D0*X))
       RETURN
    3  GEXT=X*X*(4.5603425D0*DEXP(-2.2378000D0*X)+7.9197229D0*DEXP(-3.82
     &44700D0*X)+1.1697093D0*DEXP(-1.6770200D0*X)+2.6575680D0*DEXP(-8.58
     &10500D0*X))
       RETURN
    4  GEXT=X*(3.7467103D0*DEXP(-6.4449722D0*X))+X*X*(-5.0254465D0*DEXP(
     &-2.4960885D0*X))+X*X*X*(0.5162317D0*DEXP(-1.0396983D0*X))
       RETURN
    5  GEXT=X*X*(4.0687494D0*DEXP(-2.6750348D0*X))+X*X*X*(-0.2285374D0*D
     &EXP(-0.8361112D0*X))
       RETURN
    6  GEXT=X*X*X*(0.1289291D0*DEXP(-0.7128119D0*X))
       RETURN
      END
