C UNIX-IDL PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas215/bfttyp.for,v 1.2 2011/12/02 00:17:26 mog Exp $ Date $Date: 2011/12/02 00:17:26 $
C
      SUBROUTINE BFTTYP( NDLEV  , NDTRN  ,
     &                   IZ1    , IL     ,
     &                   IA     , CSTRGA , ISA    , ILA   , XJA   , WA ,
     &                   ITRAN  , TCODE  , I1A    , I2A   , AVAL  ,
     &                   ICNTE  , ICNTP  , ICNTR  , ICNTH , ICNTI ,
     &                   IETRN  , PECODE , TECODE , IE1A  , IE2A  , AA ,
     &                   CEA
     &                 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: BFTTYP *********************
C
C  PURPOSE:  TO SORT TRANSITION ARRAYS INTO FOUR TRANSITION/RECOMB TYPES
C            AND ASSIGN INITIAL TYPES TO ELECTRON IMPACT TRANSITIONS
C
C  CALLING PROGRAM: ADAS215
C
C  SUBROUTINE:
C
C  INPUT : (I*4)  NDLEV   = MAXIMUM NUMBER OF LEVELS THAT CAN BE READ
C  INPUT : (I*4)  NDTRN   = MAX. NUMBER OF TRANSITIONS THAT CAN BE READ
C
C  INPUT : (I*4)  IZ1     = RECOMBINING ION CHARGE READ
C  INPUT : (I*4)  IL      = INPUT DATA FILE: NUMBER OF ENERGY LEVELS
C
C  INPUT : (I*4)  IA()    = ENERGY LEVEL INDEX NUMBER
C  INPUT : (C*(*))CSTRGA()= NOMENCLATURE/CONFIGURATION FOR LEVEL 'IA()'
C  INPUT : (I*4)  ISA()   = MULTIPLICITY FOR LEVEL 'IA()'
C                           NOTE: (ISA-1)/2 = QUANTUM NUMBER (S)
C  INPUT : (I*4)  ILA()   = QUANTUM NUMBER (L) FOR LEVEL 'IA()'
C  INPUT : (R*8)  XJA()   = QUANTUM NUMBER (J-VALUE) FOR LEVEL 'IA()'
C                           NOTE: (2*XJA)+1 = STATISTICAL WEIGHT
C  INPUT : (R*8)  WA()    = ENERGY RELATIVE TO LEVEL 1 (CM-1) FOR LEVEL
C                           'IA()'
C
C  INPUT : (I*4)  ITRAN   = INPUT DATA FILE: NUMBER OF TRANSITIONS
C  INPUT : (C*1)  TCODE() = TRANSITION: DATA TYPE POINTER:
C                           ' ' => Electron Impact   Transition
C                           'P' => Proton   Impact   Transition
C                           'H' => Charge   Exchange Recombination
C                           'R' => Free     Electron Recombination
C                           'I' => Electron Impact Ionisation
C  INPUT : (I*4)  I1A()   = TRANSITION:
C                            LOWER ENERGY LEVEL INDEX (CASE ' ' & 'P')
C                            NOT USED                 (CASE 'H' & 'R')
C  INPUT : (I*4)  I2A()   = TRANSITION:
C                            UPPER ENERGY LEVEL INDEX (CASE ' ' & 'P')
C                            CAPTURING    LEVEL INDEX (CASE 'H' & 'R')
C  INPUT : (R*8)  AVAL()  = TRANSITION:
C                            A-VALUE (SEC-1)          (CASE ' ')
C                            NEUTRAL BEAM ENERGY      (CASE 'H')
C                            NOT USED                 (CASE 'P' & 'R')
C
C  OUTPUT: (I*4)  ICNTE   = NUMBER OF ELECTRON IMPACT TRANSITIONS INPUT
C  OUTPUT: (I*4)  ICNTP   = NUMBER OF PROTON IMPACT TRANSITIONS INPUT
C  OUTPUT: (I*4)  ICNTR   = NUMBER OF FREE ELECTRON RECOMBINATIONS INPUT
C  OUTPUT: (I*4)  ICNTH   = NO. OF CHARGE EXCHANGE RECOMBINATIONS INPUT
C  OUTPUT: (I*4)  ICNTI   = NO. OF INNNER SHELL IONISATION INPUT
C
C  OUTPUT: (I*4)  IETRN() = INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           1ST. DIM.: EL-TRANS. INDEX
C                           REPRESENT ELECTRON IMPACT TRANSITIONS.
C  OUTPUT: (C*1)  PECODE()= ELECTRONIC TRANSITION PLOT SELECTOR:
C                           ' ' => do not plot
C                           'P' or 'p' => plot
C                           1ST. DIM.: EL-TRANS. INDEX
C  OUTPUT: (C*1)  TECODE()= ELECTRONIC TRANSITION: DATA TYPE POINTER:
C                           ' ' => unassigned
C                           '1' => dipole
C                           '2' => non-dipole, non-spin change
C                           '3' => spin change
C                           '4' => small oscillator strength
C                           1ST. DIM.: EL-TRANS. INDEX
C
C  OUTPUT: (I*4)  IE1A()  = EL-TRANS. LOWER ENERGY LEVEL INDEX
C                           1ST. DIM.: EL-TRANS. INDEX
C  OUTPUT: (I*4)  IE2A()  = EL-TRANS. UPPER ENERGY LEVEL INDEX
C                           1ST. DIM.: EL-TRANS. INDEX
C  OUTPUT: (R*8)  AA()    = EL-TRANS.  A-VALUE (SEC-1)
C                           1ST. DIM.: EL-TRANS. INDEX
C  OUTPUT: (R*8)  CEA()   = EL-TRANS. BURGESS & TULLY C-VALUE
C                           1ST. DIM.: EL-TRANS. INDEX
C
C
C          (R*8)  CEREF   = PARAMETER = REFERENCE VALUE FOR B&T C-VAL.
C          (R*8)  FZERO   = PARAMETER = EFF. ZERO FOR F-VALUES IN 
C                                       BURGESS & TULLY TYPE SELECTION.
C          (R*8)  FBIG    = PARAMETER = F-VALUE FOR TYPE SWITCH 1-4 IN 
C                                       BURGESS & TULLY TYPE SELECTION.
C
C
C ROUTINES: NONE
C
C AUTHOR:  HP SUMMERS, UNIVERSITY OF STRATHCLYDE
C          JA8.08
C          TEL.  0141-553-4196
C
C DATE  :  04/06/98
C
C UPDATE:
C
C VERSION:  1.1
C DATE   : 09/08/98
C MODIFIED: RICHARD MARTIN
C		- PUT UNDER SCCS CONTROL.							 
C
C VERSION : 1.2
C DATE    : 04-11-2011
C MODIFIED: Hugh Summers
C           - tidied cstrga string length descriptor to variable field 
C             length and verified correct in declarations.
c
C 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      REAL*8     CEREF       , FZERO     , FBIG
C-----------------------------------------------------------------------
      PARAMETER (CEREF = 1.5 , FZERO = 1.0D-4  , FBIG = 0.01 )
      INTEGER    NDLEV       , NDTRN     , 
     &           ICNTE       , ICNTP     , ICNTR       , ICNTH     ,
     &           ICNTI
      INTEGER    IZ1         , IL        , ITRAN          
      INTEGER    I           , IE
C-----------------------------------------------------------------------
      REAL*8     WVNOU       , WVNOL     , ELU         ,
     &           AIN         , S         , FIN         ,
     &           WTL         , WTU
      REAL*8     Z1
C-----------------------------------------------------------------------
      INTEGER    IA(NDLEV)   , ISA(NDLEV)     , ILA(NDLEV) ,
     &           I1A(NDTRN)  , I2A(NDTRN)     
      INTEGER    IETRN(NDTRN),
     &           IE1A(NDTRN) , IE2A(NDTRN)
C-----------------------------------------------------------------------
      REAL*8     XJA(NDLEV)  , WA(NDLEV)      
      REAL*8     AA(NDTRN)   , AVAL(NDTRN)  , CEA(NDTRN)
C-----------------------------------------------------------------------
      CHARACTER  TCODE(NDTRN)*1  , CSTRGA(NDLEV)*(*)
      CHARACTER  PECODE(NDTRN)*1 , TECODE(NDTRN)*1
C-----------------------------------------------------------------------
C
C **********************************************************************
C
C-----------------------------------------------------------------------
C ZERO TRANSITION TYPE COUNTS AND ACTIVE METASTABLE COUNTS.
C-----------------------------------------------------------------------
C
      ICNTE = 0
      ICNTP = 0
      ICNTR = 0
      ICNTH = 0
      ICNTI = 0
C
C-----------------------------------------------------------------------
C ADD TRANSITION PARAMETERS TO RELEVANT ARRAYS.
C-----------------------------------------------------------------------
C
         DO 1 I=1,ITRAN
               IF ((TCODE(I).EQ.' ').OR.
     &             (TCODE(I).EQ.'1').OR.
     &             (TCODE(I).EQ.'2').OR.
     &             (TCODE(I).EQ.'3').OR.
     &             (TCODE(I).EQ.'4')) THEN
                   ICNTE        = ICNTE + 1
                   IETRN(ICNTE) = I
                   PECODE(ICNTE)= ' '
                   TECODE(ICNTE)= TCODE(I)
                   IE1A(ICNTE)  = I1A(I)
                   IE2A(ICNTE)  = I2A(I)
                   AA(ICNTE)    = AVAL(I)
                   CEA(ICNTE)   = CEREF
               ELSEIF (TCODE(I).EQ.'P') THEN
                  ICNTP        = ICNTP + 1
               ELSEIF (TCODE(I).EQ.'R') THEN
                  ICNTR        = ICNTR + 1
               ELSEIF (TCODE(I).EQ.'H') THEN
                  ICNTH        = ICNTH + 1
               ELSEIF (TCODE(I).EQ.'I') THEN
                  ICNTI        = ICNTI + 1
               ENDIF
    1    CONTINUE
C
C-----------------------------------------------------------------------
C SET UP INITIAL CLASSIFICATION OF ELECTRON IMPACT TRANSITIONS.
C-----------------------------------------------------------------------
C
       Z1 = IZ1
C
       DO 2 IE = 1, ICNTE
C
         IF (TCODE(IETRN(IE)).EQ.' ') THEN
C
             WVNOU = WA(IE2A(IE))
             WVNOL = WA(IE1A(IE))
             ELU   = DABS(WVNOU-WVNOL)/109737.26
             WTU   = 2.0*XJA(IE2A(IE))+1.0
             WTL   = 2.0*XJA(IE1A(IE))+1.0
             AIN   = AVAL(IE)
             S     =3.73491E-10*WTU*AIN/(ELU*ELU*ELU)                        
             FIN   =3.333333D-1*ELU*S/WTL                                     
C
             IF(ISA(IE2A(IE)).EQ.ISA(IE1A(IE))) THEN
                  IF((IABS(ILA(IE2A(IE))-ILA(IE1A(IE))).LE.1).AND.
     &                    (FIN.GT.FBIG))THEN
                     TECODE(IE) = '1'
                  ELSE
                     TECODE(IE) = '2'
                  ENDIF
             ELSE
                 IF((FIN.GT.FZERO).AND.(FIN.LT.FBIG))THEN
                     TECODE(IE) = '4'
                 ELSE
                     TECODE(IE) = '3'
                 ENDIF
             ENDIF
          ENDIF
    2  CONTINUE
C
C-----------------------------------------------------------------------
      RETURN
      END
