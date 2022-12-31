       SUBROUTINE H9TRNI( NDLEV  , NDTRN  , NDTEM , ndmet ,
     &                    IL     , ISTRN  , NV    ,
     &                    IA     , WA     , XJA   ,
     &                    I1A    , I2A    , AVAL  ,
     &                    SCOM   , zpla   , bwnoa , ipla ,
     &                    IUPPER , ILOWER ,
     &                    LUPPER , LLOWER ,
     &                    WUPPER , WLOWER ,
     &                    EUPPER , ELOWER ,
     &                    AA     , GAMMA  ,
     &                    zeta   , ip
     &                  )
       IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: H9TRNI *********************
C
C  PURPOSE:  TO SET UP SELECTED IONISATION TRANSITION PARAMETERS
C
C  CALLING PROGRAM: ADAS809
C
C  SUBROUTINE:
C
C  INPUT : (I*4)  NDLEV     = MAXIMUM NUMBER OF INDEX LEVELS
C  INPUT : (I*4)  NDTRN     = MAXIMUM NUMBER OF TRANSITIONS
C  INPUT : (I*4)  NDTEM     = MAXIMUM NUMBER OF INPUT FILE TEMPERATURES
C
C  INPUT : (I*4)  IL        = NUMBER OF INDEX LEVELS
C  INPUT : (I*4)  ISTRN     = SELECTED TRANSITION INDEX.
C  INPUT : (I*4)  NV        = INPUT DATA FILE: NUMBER OF GAMMA/TEMPERATURE
C                             PAIRS FOR THE SELECTED TRANSITION.
C
C  INPUT : (I*4)  IA()      = LEVEL INDEX NUMBER ARRAY
C  INPUT : (R*8)  WA()      = LEVEL ENERGIES RELATIVE TO LEVEL 1 (CM-1)
C  INPUT : (R*8)  XJA()     = QUANTUM NUMBER (J-VALUE) FOR LEVEL
C                             NOTE: (2*XJA)+1 = STATISTICAL WEIGHT
C
C  INPUT : (I*4)  I1A()     = LOWER LEVEL INDEX FOR ELECTRON IMPACT
C                             TRANSITION
C  INPUT : (I*4)  I2A()     = UPPER LEVEL INDEX FOR ELECTRON IMPACT
C                             TRANSITION
C  INPUT : (I*4)  AVAL()    = A-VALUE FOR ELECTRON IMPACT TRANSITION
C  INPUT : (I*4)  SCOM(,)   = GAMMA VALUES FOR ELECTRON IMPACT 
C                             (DE-)EXCITATION
C                             1st DIMENSION: TEMPERATURE INDEX
C                             2nd DIMENSION: TRANSITION INDEX
C  OUTPUT: (I*4)  IUPPER    = SELECTED TRANSITION: UPPER LEVEL ARRAY INDEX
C  OUTPUT: (I*4)  ILOWER    = SELECTED TRANSITION: LOWER LEVEL ARRAY INDEX
C
C
C  OUTPUT: (I*4)  LUPPER    = SELECTED TRANSITION: UPPER INDEX LEVEL
C  OUTPUT: (I*4)  LLOWER    = SELECTED TRANSITION: LOWER INDEX LEVEL
C
C  OUTPUT: (R*8)  WUPPER    = SELECTED TRANSITION: UPPER LEVEL STAT. WT.
C  OUTPUT: (R*8)  WLOWER    = SELECTED TRANSITION: LOWER LEVEL STAR. WT.
C                             (NOTE: STAT. WT. = STATISTICAL WEIGHT)
C
C  OUTPUT: (R*8)  EUPPER    = SELECTED TRANSITION: UPPER ENERGY LEVEL
C                             RELATIVE TO INDEX LEVEL 1. (CM-1)
C  OUTPUT: (R*8)  ELOWER    = SELECTED TRANSITION: LOWER ENERGY LEVEL
C                             RELATIVE TO INDEX LEVEL 1. (CM-1)
C  OUTPUT: (R*8)  AA        = SELECTED TRANSITION A-VALUE (SEC-1)
C  OUTPUT: (R*8)  GAMMAUP() = INPUT DATA FILE: SELECTED EXCITATION -
C                                              GAMMAUP VALUE AT 'TEMP()'
C  OUTPUT: (R*8)  GAMMADN() = INPUT DATA FILE: SELECTED DE-EXCITATION -
C                                              GAMMADN VALUE AT 'TEMP()'
C
C          (I*4)  I         = GENERAL USE.
C
C ROUTINES: NONE
C
C AUTHOR:  HUGH SUMMERS (UNIVERSITY OF STRATHCLYDE)
C          JA7.08
C          EXT. 4196
C
C DATE:    30/11/01
C
C UPDATE:  Paul Bryans
C	   24/11/04
C	   Added extra parameters needed for ionisation transition
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER  NDLEV  , NDTRN  , NDTEM , ndmet ,
     &         IL     , ISTRN  , NV    ,
     &         IUPPER , ILOWER ,
     &         LUPPER , LLOWER
      INTEGER  I      , ipla(ndmet,ndlev)
C-----------------------------------------------------------------------
      REAL*8   WUPPER , WLOWER ,
     &         EUPPER , ELOWER ,
     &         AA
C-----------------------------------------------------------------------
      INTEGER  IA(NDLEV)   ,
     &         I1A(NDTRN)  , I2A(NDTRN)
C-----------------------------------------------------------------------
      REAL*8   WA(NDLEV)   , XJA(NDLEV)   ,
     &         AVAL(NDTRN) , GAMMA(NDTEM) ,
     &         SCOM(NDTEM,NDTRN)
      real*8   bwnoa(ndmet), zpla(ndmet,ndlev)
      real*8   zeta, ip
C-----------------------------------------------------------------------
C
C **********************************************************************
C
C-----------------------------------------------------------------------
C SET UPPER AND LOWER ENERGY LEVEL INDEX SPECIFICATIONS
C-----------------------------------------------------------------------
C
      LLOWER=I1A(ISTRN)
      LUPPER=I2A(ISTRN)
         DO 1 I=1,IL
           IF ( IA(I).EQ.LLOWER ) ILOWER=I
           IF ( IA(I).EQ.LUPPER ) IUPPER=I
    1    CONTINUE
C
C-----------------------------------------------------------------------
C SET UP REQUIRED A-VALUE AND GAMMA COEFFICIENTS.
C-----------------------------------------------------------------------
C
      AA   = AVAL(ISTRN)
         DO 2 I=1,NV
           GAMMA(I) = SCOM(I,ISTRN)
    2    CONTINUE
C
C-----------------------------------------------------------------------
C ASSIGN INDEX LEVEL ENERGIES AND STATISTICAL WEIGHTS.
C-----------------------------------------------------------------------
C
      ELOWER=WA(ILOWER)
      EUPPER=WA(IUPPER)
      WLOWER=(2.d0*XJA(ILOWER))+1.d0
      WUPPER=(2.d0*XJA(IUPPER))+1.d0

C-----------------------------------------------------------------------
      zeta = zpla(ilower,iupper)
      ip = bwnoa(ipla(ilower,iupper)) - wa(iupper)
      


C-----------------------------------------------------------------------
      RETURN
      END
