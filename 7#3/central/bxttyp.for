CX UNIX PORT - SCCS Info : Module @(#)$Header: /home/adascvs/fortran/adas2xx/adaslib/bxttyp.for,v 1.3 2005/03/17 11:37:40 allan Exp $ Date $Date: 2005/03/17 11:37:40 $
CX
      SUBROUTINE BXTTYP( NDLEV  , NDMET  , NDTRN  , NPLR  , NPLI  ,
     &                   ITRAN  , TCODE  , I1A    , I2A   , AVAL  ,
     &                   ICNTE  , ICNTP  , ICNTR  , ICNTH , ICNTI , 
     &                   ICNTL  , ICNTS  ,
     &                   IETRN  , IPTRN  , IRTRN  , IHTRN , IITRN ,
     &                   ILTRN  , ISTRN  ,
     &                                     IE1A   , IE2A  , AA    ,
     &                                     IP1A   , IP2A  ,
     &                                     IA1A   , IA2A  , AUGA  ,
     &                                     IL1A   , IL2A  , WVLA  ,
     &                                     IS1A   , IS2A  , LSS04A  
     &                 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: BXTTYP *********************
C
C  PURPOSE:  TO SORT TRANSITION ARRAYS INTO SEVEN TRANSITION/RECOMB TYPES
C
C  CALLING PROGRAM: General
C
C  SUBROUTINE:
C
C  INPUT : (I*4)  NDLEV   = MAXIMUM NUMBER OF LEVELS THAT CAN BE READ
C  INPUT : (I*4)  NDMET   = MAXIMUM NUMBER OF METASTABLES
C  INPUT : (I*4)  NDTRN   = MAXIMUM NUMBER OF TRANS. THAT CAN BE READ
C
C  OUTPUT: (I*4)  NPLR    = NO. OF ACTIVE METASTABLES OF (Z+1) ION
C  OUTPUT: (I*4)  NPLI    = NO. OF ACTIVE METASTABLES OF (Z-1) ION
C
C  INPUT : (I*4)  ITRAN   = INPUT DATA FILE: NUMBER OF TRANSITIONS
C  INPUT : (C*1)  TCODE() = TRANSITION: DATA TYPE POINTER:
C                           ' ' => Electron Impact   Transition
C                           'P' => Proton   Impact   Transition
C                           'H' => Charge   Exchange Recombination
C                           'R' => Free     Electron Recombination
C                           'I' => Electron Impact Ionisation to z
C                           'L' => Satellites from DR Recombination
C                           'S' => Electron Impact Ionisation to z+1 
C  INPUT : (I*4)  I1A()   = TRANSITION:
C                            LOWER ENERGY LEVEL INDEX (CASE ' ' & 'P')
C                            PARENT ENERGY LEVEL INDEX(CASE 'H' & 'R')
C                                                     (         & 'L')
C                            FINAL PARENT LEVEL INDEX (CASE 'S')
C  INPUT : (I*4)  I2A()   = TRANSITION:
C                            UPPER ENERGY LEVEL INDEX (CASE ' ' & 'P')
C                            CAPTURING    LEVEL INDEX (CASE 'H' & 'R')
C                                                     (         & 'L')
C                            IONISING     LEVEL INDEX (CASE 'S')
C  INPUT : (R*8)  AVAL()  = TRANSITION:
C                            A-VALUE (SEC-1)          (CASE ' ')
C                            NEUTRAL BEAM ENERGY      (CASE 'H')
C                            AUGER VALUE(SEC-1)       (CASE 'R')
C                            PARENT WAVLENGTH (A)     (CASE 'L')
C                            NOT USED                 (CASE 'P' & 'S')
C
C  OUTPUT: (I*4)  ICNTE   = NUMBER OF ELECTRON IMPACT TRANSITIONS INPUT
C  OUTPUT: (I*4)  ICNTP   = NUMBER OF PROTON IMPACT TRANSITIONS INPUT
C  OUTPUT: (I*4)  ICNTR   = NUMBER OF FREE ELECTRON RECOMBINATIONS INPUT
C  OUTPUT: (I*4)  ICNTH   = NO. OF CHARGE EXCHANGE RECOMBINATIONS INPUT
C  OUTPUT: (I*4)  ICNTI   = NO. OF IONISATIONS TO Z INPUT
C  OUTPUT: (I*4)  ICNTL   = NO. OF SATELLITE DR RECOMBINATIONS INPUT
C  OUTPUT: (I*4)  ICNTS   = NO. OF IONISATIONS TO Z+1 INPUT
C
C  OUTPUT: (I*4)  IETRN() = ELECTRON IMPACT TRANSITION:
C                           INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           REPRESENT ELECTRON IMPACT TRANSITIONS.
C  OUTPUT: (I*4)  IPTRN() = PROTON IMPACT TRANSITION:
C                           INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           REPRESENT PROTON IMPACT TRANSITIONS.
C  OUTPUT: (I*4)  IRTRN() = FREE ELECTRON RECOMBINATION:
C                           INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           REPRESENT FREE ELECTRON RECOMBINATIONS.
C  OUTPUT: (I*4)  IHTRN() = CHARGE EXCHANGE RECOMBINATION:
C                           INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           REPRESENT CHARGE EXCHANGE RECOMBINATIONS.
C  OUTPUT: (I*4)  IITRN() = ELECTRON IMPACT IONISATION:
C                           INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           REPRESENT IONISATIONS FROM LOWER STAGE ION.
C  OUTPUT: (I*4)  ILTRN() = SATELLITE DR RECOMBINATION:
C                           INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           REPRESENT SATELLITE DR RECOMBINATIONS.
C  OUTPUT: (I*4)  ISTRN() = ELECTRON IMPACT IONISATION:
C                           INDEX VALUES IN MAIN TRANSITION ARRAYS WHICH
C                           REPRESENT IONISATIONS TO UPPER STAGE ION.
C
C  OUTPUT: (I*4)  IE1A()  = ELECTRON IMPACT TRANSITION:
C                            LOWER ENERGY LEVEL INDEX
C  OUTPUT: (I*4)  IE2A()  = ELECTRON IMPACT TRANSITION:
C                            UPPER ENERGY LEVEL INDEX
C  OUTPUT: (R*8)  AA()    = ELECTRON IMPACT TRANSITION: A-VALUE (SEC-1)
C
C
C  OUTPUT: (I*4)  IP1A()  = PROTON IMPACT TRANSITION:
C                            LOWER ENERGY LEVEL INDEX
C  OUTPUT: (I*4)  IP2A()  = PROTON IMPACT TRANSITION:
C                            UPPER ENERGY LEVEL INDEX
C
C  OUTPUT: (I*4)  IA1A()  = AUGER TRANSITION:
C                            PARENT ENERGY LEVEL INDEX
C  OUTPUT: (I*4)  IA2A()  = AUGER TRANSITION:
C                            RECOMBINED ION ENERGY LEVEL INDEX
C  OUTPUT: (R*8)  AUGA()  = AUGER TRANSITION: AUG-VALUE (SEC-1)
C                            RECOMBINED ION ENERGY LEVEL INDEX
C  OUTPUT: (I*4)  IL1A()  = SATELLITE DR TRANSITION:
C                            RECOMBINING ION  INDEX
C  OUTPUT: (I*4)  IL2A()  = SATELLITE DR TRANSITION:
C                            RECOMBINED ION INDEX
C  OUTPUT: (R*8)  WVLA()  = SATELLITE DR TRANSITION: PARENT WVLGTH.(A)
C                            DR SATELLITE LINE INDEX
C  OUTPUT: (I*4)  IS1A()  = IONISING TRANSITION:
C                            IONISED ION  INDEX
C  OUTPUT: (I*4)  IS2A()  = IONISING TRANSITION:
C                            IONISING ION INDEX
C  OUTPUT: (L*4)  LSS04A(,)= .TRUE. => IONIS. RATE SET IN ADF04 FILE:
C                            .FALSE.=> NOT SET IN ADF04 FILE
C                            1ST DIM: LEVEL INDEX
C                            2ND DIM: PARENT METASTABLE INDEX
C
C          (I*4)  I       = GENERAL USE.
C
C
C ROUTINES: NONE
C
C AUTHOR:  HP SUMMERS (REVISION OF BXTTYP BY PE BRIDEN)
C          K1/1/57
C          JET EXT. 4941
C
C DATE  :  11/06/92
C
C-----------------------------------------------------------------------
C PUT UNDER SCCS CONTROL: 
C
C VERSION: 1.1				DATE: 10/05/96
C MODIFIED: WILLIAM OSBORN (TESSELLA SUPPORT SERVICES PLC)
C	    - FIRST PUT UNDER SCCS
C 
C VERSION: 1.2				DATE: 13/09/99
C MODIFIED: HUGH SUMMERS, UNIVERSITY OF STRATHCLYDE
C	    - ADDED DETECTION OF L-LINES AND S-LINES
C 
C-----------------------------------------------------------------------
C
C VERSION: 1.2				DATE: 01/05/2003
C MODIFIED: Martin O'Mullane
C	    - Replaced original bxttyp with b8ttyp version 1.2. 
C		Hence the 1.2 version no.
C 
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
C VERSION: 1.3				DATE: 17/03/2005
C MODIFIED: Allan Whiteford
C	    - Made the routine accept that transition codes of '1',
C             '2' and '3' as well as ' ' correspond to electron
C             impact excitation.
C 
C-----------------------------------------------------------------------
      INTEGER    NDLEV       , NDMET          , NDTRN      ,       
     &           NPLR        , NPLI           , ITRAN      ,
     &           ICNTE       , ICNTP          , ICNTR      , ICNTH     ,
     &           ICNTI       , ICNTL          , ICNTS
      INTEGER    I           , J
C-----------------------------------------------------------------------
      INTEGER    I1A(NDTRN)  , I2A(NDTRN)
      INTEGER    IETRN(NDTRN), IPTRN(NDTRN)   ,
     &           IRTRN(NDTRN), IHTRN(NDTRN)   ,
     &           IITRN(NDTRN), ILTRN(NDTRN)   ,
     &           ISTRN(NDTRN),
     &           IE1A(NDTRN) , IE2A(NDTRN)    ,
     &           IP1A(NDTRN) , IP2A(NDTRN)    ,
     &           IA1A(NDTRN) , IA2A(NDTRN)    ,
     &           IL1A(NDLEV) , IL2A(NDLEV)    ,
     &           IS1A(NDLEV) , IS2A(NDLEV)
C-----------------------------------------------------------------------
      REAL*8     AA(NDTRN)   , AVAL(NDTRN)    , AUGA(NDTRN) ,
     &           WVLA(NDLEV)
C-----------------------------------------------------------------------
      CHARACTER  TCODE(NDTRN)*1
C-----------------------------------------------------------------------
      LOGICAL    LSS04A(NDLEV,NDMET)
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C ZERO TRANSITION TYPE COUNTS AND ACTIVE METASTABLE COUNTS.
C-----------------------------------------------------------------------

      ICNTE = 0
      ICNTP = 0
      ICNTR = 0
      ICNTH = 0
      ICNTI = 0
      ICNTL = 0
      ICNTS = 0

      NPLR  = 0
      NPLI  = 0
      
        DO I=1,NDLEV
          DO J=1,NDMET
            LSS04A(I,J)=.FALSE.
          ENDDO
        ENDDO  

C-----------------------------------------------------------------------
C ADD TRANSITION PARAMETERS TO RELEVANT ARRAYS.
C-----------------------------------------------------------------------

         DO 20 I=1,ITRAN
               IF     (TCODE(I).EQ.' ' .OR. TCODE(I).EQ.'1'
     &         .OR.    TCODE(I).EQ.'2' .OR. TCODE(I).EQ.'3') THEN
                  ICNTE        = ICNTE + 1
                  IETRN(ICNTE) = I
                  IE1A(ICNTE)  = I1A(I)
                  IE2A(ICNTE)  = I2A(I)
                  AA(ICNTE)    = AVAL(I)
               ELSEIF (TCODE(I).EQ.'P') THEN
                  ICNTP        = ICNTP + 1
                  IPTRN(ICNTP) = I
                  IP1A(ICNTP)  = I1A(I)
                  IP2A(ICNTP)  = I2A(I)
               ELSEIF (TCODE(I).EQ.'R') THEN
                  ICNTR        = ICNTR + 1
                  IRTRN(ICNTR) = I
                  NPLR=MAX0(NPLR,I1A(I))
                  IA1A(ICNTR)  = I1A(I)
                  IA2A(ICNTR)  = I2A(I)
                  AUGA(ICNTR)  = AVAL(I)
               ELSEIF (TCODE(I).EQ.'H') THEN
                  ICNTH        = ICNTH + 1
                  IHTRN(ICNTH) = I
                  NPLR=MAX0(NPLR,I1A(I))
               ELSEIF (TCODE(I).EQ.'I') THEN
                  ICNTI        = ICNTI + 1
                  IITRN(ICNTI) = I
                  NPLI=MAX0(NPLI,IABS(I1A(I)))
               ELSEIF (TCODE(I).EQ.'L') THEN
                  ICNTL        = ICNTL + 1
                  ILTRN(ICNTL) = I
                  IL1A(ICNTL)  = I1A(I)
                  IL2A(ICNTL)  = I2A(I)
                  WVLA(I2A(I))  = AVAL(I)
                  NPLI=MAX0(NPLI,IABS(I1A(I)))
               ELSEIF (TCODE(I).EQ.'S') THEN
                  ICNTS        = ICNTS + 1
                  ISTRN(ICNTS) = I
                  IS1A(ICNTS)  = I1A(I)
                  IS2A(ICNTS)  = I2A(I)
                  LSS04A(I2A(I),I1A(I))=.TRUE.
               ENDIF
   20    CONTINUE
C-----------------------------------------------------------------------
      RETURN
      END
