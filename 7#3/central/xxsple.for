CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/maths/xxsple.for,v 1.4 2007/04/11 13:01:53 allan Exp $ Date $Date: 2007/04/11 13:01:53 $
CX
      SUBROUTINE XXSPLE( LSETX , IOPT   , FINTX ,
     &                   NIN   , XIN    , YIN   ,
     &                   NOUT  , XOUT   , YOUT  ,
     &                   DY    , LINTRP
     &                 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXSPLE *********************
C
C  PURPOSE:           TO INTERPOLATE/EXTRAPOLATE USING CUBIC SPLINES
C
C                     (IF IOPT < 0 NO EXTRAPOLATION TAKES PLACE = VALUES
C                      SET TO ZERO).- LOGICAL ARRAY 'LINTRP()' SPECIFIES
C                      WHETHER OUTPUT SPLINE IS INTERPOLATED '.TRUE.' OR
C                      EXTRAPOLATED '.FALSE.'.
C
C                      (AS FOR 'XXSPLN' EXCEPT 'LINTRP' ARGUMENT ADDED).
C
C  CALLING PROGRAMS:  GENERAL USE
C
C  SUBROUTINE:
C
C  I/O   : (L*4)  LSETX   = .TRUE.  => SET UP SPLINE PARAMETERS RELATING
C                                      TO 'XIN' AXIS.
C                           .FALSE. => DO NOT SET UP SPLINE PARAMETERS
C                                      RELATING TO 'XIN' AXIS.
C                                      (I.E. THEY WERE SET IN A PREVIOUS
C                                            CALL )
C                           ( 'LSETX' IS ALWAYS RETURN AS '.FALSE.'  ON
C                             RETURN FROM THE SUBROUTINE ).
C                           ** IMPORTANT: SEE NOTES BELOW ON 'LSETX' **
C  INPUT : (I*4)  IOPT    = SPLINE END CONDITIONS/EXTRAPOLATION CONTROL
C                           SWITCH - SEE NOTES BELOW
C                           I.E. DEFINES THE BOUNDARY DERIVATIVES.
C                                (VALID VALUES = 0, 1, 2, 3, 4)
C                           IF IOPT < 0 THEN NO EXTRAPOLATION TAKES
C                           - ANY VALUES REQUIRING EXTRAPOLATION WILL BE
C                             SET TO ZERO (END CONDITIONS AS FOR IOPT=0)
C  INPUT : (R*8)  FINTX   = INTERPOLATING X-COORDINATE TRANSFORMATION.
C                           EXTERNAL FUNCTION (SEE ROUTINES BELOW)
C
C  INPUT : (I*4)  NIN     = NUMBER OF KNOTS
C  INPUT : (R*8)  XIN()   = X-VALUES OF KNOTS
C  INPUT : (R*8)  YIN()   = Y-VALUES OF KNOTS
C
C  INPUT : (I*4)  NOUT    = NUMBER OF OUTPUT VALUES TO BE INTERPOLATED
C                           EXTRAPOLATED.
C  INPUT : (R*8)  XOUT()  = X-VALUES AT WHICH INTERPOLATION/EXTRAPOLA-
C                           TION REQUIRED
C  OUTPUT: (R*8)  YOUT()  = INTERPOLATED/EXTRAPOLATED Y-VALUES FOR
C                           REQUESTED 'XOUT()' VALUES.
C
C  OUTPUT: (R*8)  DY()    = DERIVATIVES AT INPUT KNOTS (ARRAY SIZE: NIN)
C  OUTPUT: (L*4)  LINTRP()= .TRUE.  => 'YOUT()' VALUE INTERPOLATED.
C                           .FALSE. => 'YOUT()' VALUE EXTRAPOLATED.
C                           (ARRAY SIZE: NOUT)
C
C          (I*4)  NKNOTS  = PARAMETER = MAXIMUM  NUMBER OF KNOTS ALLOWED
C          (I*4)  NIOPT   = PARAMETER = MAXIMUM  VALUE OF IOPT ALLOWED
C
C          (I*4)  I       = GENERAL ARRAY USE
C          (I*4)  K       = INDEX OF 'XOUT()' VALUE FOR INTERPOLATION/
C                           EXTRAPOLATION.
C          (I*4)  NIN0    = 'NIN' - 1
C          (I*4)  INTER   = INDEX OF CLOSEST/NEXT HIGHEST VALUE OF
C                           'XIN()' TO THE VALUE OF 'XOUT()' BEING
C                           INTERPOLATED/EXTRAPOLATED.
C          (I*4)  NOPT    = VALUE OF  'IOPT' USED IN  CALCULATING  END-
C                           CONDITIONS   FOR  STORED  'X-VALUE'  SPLINE
C                           PARAMETERS.   (NOTE:  IF  'IOPT < 0',  THEN
C                           'NOPT = 0'.) - I.E. 'NOPT = MAX( 0, IOPT )'.
C
C          (R*8)  XK      = VALUE OF 'XOUT(K)' BEING INTERPOLATED/
C                           EXTRAPOLATED
C          (R*8)  XKK     = TRANSFORMED VALUE OF 'XOUT(K)' BEING
C                           INTERPOLATED/EXTRAPOLATED.
C          (R*8)  T1      = INVERSE OF SEPARATION OF KNOTS EITHER
C                           SIDE OF CURRENT KNOT.
C          (R*8)  T2      = (CURRENT KNOT POSITION TO NEXT HIGHEST KNOT
C                            POSITION) DIVIDED BY 'T1'
C          (R*8)  T3      = (CURRENT KNOT POSITION TO NEXT LOWEST  KNOT
C                            POSITION) DIVIDED BY 'T1'
C          (R*8)  T4      = INTERPOLATION FACTOR FOR CURRENT KNOT
C          (R*8)  DL1     = (REQUESTED 'XOUT()' VALUE TO NEXT HIGHEST
C                            KNOT POSITION) DIVIDED BY SEPERATION OF
C                            KNOTS EITHER SIDE OF 'XOUT(K)'.
C          (R*8)  DL2     = (REQUESTED 'XOUT()' VALUE TO NEXT LOWEST
C                            KNOT POSITION) DIVIDED BY SEPERATION OF
C                            KNOTS EITHER SIDE OF 'XOUT(K)'.
C          (R*8)  DL2     = (REQUESTED 'XOUT()' VALUE TO NEXT LOWEST
C          (R*8)  DL3     =  SEPERATION OF KNOTS EITHER SIDE OF
C                            'XOUT(K)' * 'DL1' * 'DL2'.
C
C          (L*4)  LEXTRP  = .TRUE.  => 'EXTRAPOLATION SWITCHED ON'.
C                           .FALSE. => 'EXTRAPOLATION SWITCHED OFF'.
C
C          (R*8)  QVAL()  = VALUE OF 'Q(1)'   : FUNCTION OF 'NOPT'
C          (R*8)  D2VAL() = VALUE OF 'D2(1)'  : FUNCTION OF 'NOPT'
C          (R*8)  D3VAL() = VALUE OF 'D3(1)'  : FUNCTION OF 'NOPT'
C          (R*8)  UVAL() =  VALUE OF 'U(NIN)' : FUNCTION OF 'NOPT'
C          (R*8)  AGRL() =  POLYNOMIAL CONSTANTS FOR CUBIC SPLINE FOR
C                           GIVEN 'XOUT(K)' VALUE.
C          (R*8)  X()    =  TRANSFORMED VALUES OF 'XIN()'
C          (R*8)  H()    =  SEPERATION, ALONG X-AXIS, OF KNOT FROM NEXT
C                           HIGHEST KNOT.
C          (R*8)  Q()    =  SECOND DERIVATIVE FOR KNOT
C          (R*8)  U()    =  TEMPORARY STORAGE OF DECOMPOSED FACTORS
C          (R*8)  DELY() =  SEPERATION, ALONG Y-AXIS, OF KNOT FROM NEXT
C                           HIGHEST KNOT.
C          (R*8)  D1()   =  MULTIPLICATION FACTOR USED IN CALCULATING
C                           'U()'.
C          (R*8)  D2()   =  MULTIPLICATION FACTOR USED IN CALCULATING
C                           'U()'.
C          (R*8)  D3()   =  MULTIPLICATION FACTOR USED IN CALCULATING
C                           'U()'.
C
C          (L*4)  LUVAL()=  .TRUE. => VALUE OF 'UVAL()' REFERS TO RATE
C                                     OF CHANGE OF SLOPE AT FINAL POINT.
C                           .FALSE.=> VALUE OF 'UVAL()' REFERS TO FINAL
C                                     SLOPE
C                            FUNCTION OF 'NOPT'
C
C NOTES: 'LSETX': SET TO .TRUE. ON ENTRY IF A NEW 'XIN' ARRAY IS BEING
C                 USED.  IF THE 'XIN' AXIS IS THE SAME FOR A NUMBER OF
C                 CALLS THEN DO NOT RESET 'LSETX'  -  THIS  SUBROUTINE
C                 SETS IT TO .FALSE. FOR YOU.   IF THE VALUE OF 'NOPT'
C                 IS CHANGED BETWEEN CALLS THEN THE VALUE  OF  'LSETX'
C                 ON  ENTRY IS TAKEN AS BEING EQUAL TO .TRUE. .
C
C                 THEREFORE 'LSETX' NEED ONLY BE SET TO .TRUE. ON ENTRY
C                 IF EITHER IT IS ITS FIRST CALL OR IF ANY ONE  OF  THE
C                 FOLLOWING VALUES HAS CHANGED:
C
C                 'NIN' , 'FINTX' , 'XIN(I), I=1,NIN'
C
C                 CARE: A VARIABLE MUST BE USED FOR 'LSETX', A CONSTANT,
C                       I.E.  .TRUE. ,  CANNOT BE DIRECTLY TYPED AS  AN
C                       ARGUMENT BECAUSE IT WILL BE CHANGED TO  .FALSE.
C                       ON RETURN.
C
C         SPLINE  END CONDITIONS AND EXTRAPOLATION DEPEND ON 'IOPT' AS
C         FOLLOWS:
C
C         --------------------------------------------------------------
C         | IOPT  | NOPT |  DY(1)  DDY(1)  |  DY(N)   DDY(N)  |EXTRAP'N|
C         |-------|------|-----------------|------------------|--------|
C         | < 0   |   0  |    -     0.0    |    -      0.0    |  NO    |
C         |   0   |   0  |    -     0.0    |    -      0.0    |  YES   |
C         |   1   |   1  |    -     0.0    |  -1.5      -     |  YES   |
C         |   2   |   2  |   0.0     -     |   1.0      -     |  YES   |
C         |   3   |   3  |  -0.5     -     |  -1.5      -     |  YES   |
C         |   4   |   4  |   0.0     -     |    -      0.0    |  YES   |
C         |   5   |   5  |  -4.5     -     |  -1.5      -     |  YES   |
C         |   6   |   6  |  +0.5     -     |    -      0.0    |  YES   |
C         |   7   |   7  |  -3.5     -     |    -      0.0    |  YES   |
C         --------------------------------------------------------------
C
C            NB. OPTIONS TO BE EXTENDED FOR POWER AND CX APPLICATION
C
C         -------------------------------------------------------------
C          IF ( IOPT.LT.0 ) - NO EXTRAPOLATION TAKES PLACE VALUES SET
C                             TO ZERO (CARE IF LOG OF OUTPUT IS NEEDED).
C          IF ( IOPT.GT.7 ) PROGRAM STOPS
C         -------------------------------------------------------------
C
C          THIS SUBROUTINE IS AN AMENDED  AND STRUCTURED VERSION OF  THE
C          SUBROUTINE  'ESPLINE'  WRITTEN BY  H.P. SUMMERS,   JET   26TH
C          OCTOBER 1989.   IT REMOVES THE COMMON BLOCK  /IONSPL/ ,   THE
C          SWITCHES 'ISW & ISW2' AND ALSO THE CASE FOR THE INTERPOLATION
C          OF CHARGE STATE VALUES.   IT INTRODUCES THE FEATURE  THAT  AN
C          ARRAY OF INPUT  'X-VALUES'  CAN BE  INTERPOLATED/EXTRAPOLATED
C          IN ONE CALL.
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          FINTX      ------    EXTERNAL  REAL*8  FUNCTION,   USED  TO
C                               TRANSFORM X-COORDINATES.
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     14/01/91 - ADAS91: AS FOR 'XXSPLN' BUT WITH 'LINTRP()' ADDED
C
C VERSION: 	1.2
C 
C MODIFIED: LORNE HORTON (JET) 				DATE: 25/10/97
C           - ADDED IOPT CHOICES 5, 6 AND 7
C
C VERSION: 	1.3
C 
C MODIFIED: Martin O'Mullane (JET) 			DATE: 2/6/99
C           - SAVE nin0 and inter variables also. All compilers, ie
C             especially g77, do not automatically save (or initialise 
C             variables to zero).
C
C VERSION  : 1.4                         
C DATE     : 10-04-2007
C MODIFIED : Allan Whiteford
C               - Modified documentation as part of automated
C		  subroutine documentation preparation.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER    NKNOTS       , NIOPT
C-----------------------------------------------------------------------
      PARAMETER( NKNOTS = 101 , NIOPT = 7 )
C-----------------------------------------------------------------------
      INTEGER    IOPT         , NIN       , NOUT      , NOPT
      INTEGER    I            , NIN0      , K         , INTER
C-----------------------------------------------------------------------
      REAL*8     FINTX
      REAL*8     XK           , XKK
      REAL*8     T1           , T2         , T3       , T4   ,
     &           DL1          , DL2        , DL3
C-----------------------------------------------------------------------
      LOGICAL    LSETX        , LEXTRP
C-----------------------------------------------------------------------
      REAL*8     XIN(NIN)       , YIN(NIN)       ,
     &           XOUT(NOUT)     , YOUT(NOUT)     ,
     &           DY(NIN)
      REAL*8     QVAL(0:NIOPT)  , D2VAL(0:NIOPT) ,
     &           D3VAL(0:NIOPT) , UVAL(0:NIOPT)
      REAL*8     AGRL(4)
      REAL*8     X(NKNOTS)      , DELY(NKNOTS)   ,
     &           H(NKNOTS)      , Q(NKNOTS)      , U(NKNOTS)  ,
     &           D1(NKNOTS)     , D2(NKNOTS)     , D3(NKNOTS)
C-----------------------------------------------------------------------
      LOGICAL    LINTRP(NOUT)   , LUVAL(0:NIOPT)
C-----------------------------------------------------------------------
      DATA  QVAL / -0.5    , -0.5    ,  0.0    ,  0.0    ,  0.0    ,
     &              0.0    ,  0.0    ,  0.0    / ,
     &      D2VAL/  1.5    ,  1.5    ,  0.0    ,  0.0    ,  0.0    ,
     &              0.0    ,  0.0    ,  0.0    / ,
     &      D3VAL/  0.0    ,  0.0    ,  0.0    , -0.5    ,  0.0    ,
     &             -4.5    ,  0.5    , -3.5    / ,
     &      UVAL /  0.0    , -1.5    ,  1.0    , -1.5    ,  0.0    ,
     &             -1.5    ,  0.0    ,  0.0    /
      DATA  LUVAL/ .TRUE.  , .FALSE. , .FALSE. , .FALSE. , .TRUE.  ,
     &             .FALSE. , .TRUE.  , .TRUE.  /
      DATA  NOPT / -1     /
C-----------------------------------------------------------------------
      SAVE QVAL , D2VAL , D3VAL , UVAL , LUVAL
      SAVE NOPT
      SAVE X    , H     , Q     , D1   , D2   , D3
      SAVE nin0 , inter
C-----------------------------------------------------------------------
C***********************************************************************
      IF (NKNOTS.LT.NIN)
     &                 STOP ' XXSPLE ERROR: TOO MANY KNOTS REQUESTED'
      IF (IOPT.GT.NIOPT)
     &                 STOP ' XXSPLE ERROR: INVALID IOPT VALUE ENTERED'
C-----------------------------------------------------------------------
      IF (IOPT.LT.0) THEN
         LEXTRP = .FALSE.
         IF (NOPT.NE.0)    LSETX=.TRUE.
      ELSE
         LEXTRP = .TRUE.
         IF (IOPT.NE.NOPT) LSETX=.TRUE.
      ENDIF
C-----------------------------------------------------------------------
         IF (LSETX) THEN
C
C***********************************************************************
C SET UP PARAMETERS RELATING TO 'XIN'-AXIS
C***********************************************************************
C
            NOPT = MAX0( 0 , IOPT )
C-----------------------------------------------------------------------
               DO 1 I=1,NIN
                  X(I)=FINTX(XIN(I))
    1          CONTINUE
C-----------------------------------------------------------------------
            H(2)  = X(2)-X(1)
            Q(1)  = QVAL(NOPT)
            D2(1) = D2VAL(NOPT)/H(2)
            D3(1) = D3VAL(NOPT)
            NIN0=NIN-1
C-----------------------------------------------------------------------
               DO 2 I=2,NIN0
                  H(I+1) = X(I+1) - X(I)
                  T1     = 1.0 / ( H(I+1) + H(I) )
                  T2     = H(I+1) * T1
                  T3     = 1.0 - T2
                  T4     = 1.0 / ( ( T2 * Q(I-1) ) + 2.0 )
                  Q(I)   = -T3 * T4
                  D1(I)  = ( 3.0 * T4 * T2 ) / H(I)
                  D2(I)  = ( 3.0 * T4 * T3 ) / H(I+1)
                  D3(I)  =  T2 * T4
    2          CONTINUE
C-----------------------------------------------------------------------
            T4      = 1.0 / ( Q(NIN0) + 2.0 )
            D1(NIN) = ( 3.0 * T4 ) / H(NIN)
            D3(NIN) = T4
C-----------------------------------------------------------------------
            LSETX=.FALSE.
C***********************************************************************
         ENDIF
C***********************************************************************
C SET UP CUBIC SPLINE DERIVATIVES
C***********************************************************************
      DELY(2) = YIN(2) - YIN(1)
      U(1)    = ( D2(1) * DELY(2) ) + D3(1)
C-----------------------------------------------------------------------
         DO 3 I=2,NIN0
            DELY(I+1) = YIN(I+1) - YIN(I)
            U(I) = (D1(I)*DELY(I)) + (D2(I)*DELY(I+1)) - (D3(I)*U(I-1))
    3    CONTINUE
C-----------------------------------------------------------------------
         IF (LUVAL(NOPT)) THEN
            U(NIN) = ( D1(NIN)*DELY(NIN) ) - ( D3(NIN)*U(NIN0) )
         ELSE
            U(NIN) = UVAL(NOPT)
         ENDIF
C-----------------------------------------------------------------------
      DY(NIN) = U(NIN)
C-----------------------------------------------------------------------
         DO 4 I=NIN0,1,-1
            DY(I)  = ( Q(I) * DY(I+1) ) + U(I)
    4    CONTINUE
C***********************************************************************
C SET UP PARAMETERS RELATING TO THE REQUESTED 'XOUT' ARRAY VALUES
C***********************************************************************
         DO 5 K=1,NOUT
            XK  = XOUT(K)
            XKK = FINTX(XK)
C-----------------------------------------------------------------------
C EXTRAPOLATE: HIGH 'XOUT' VALUE - IF EXTRAPOLATION SWITCHED ON
C-----------------------------------------------------------------------
               IF     (XK.GE.XIN(NIN)) THEN
                  INTER     = NIN
                  LINTRP(K) = .FALSE.
                  AGRL(1)   = 0.0
                  AGRL(3)   = 0.0
                     IF (XK.EQ.XIN(NIN)) THEN
                        LINTRP(K) = .TRUE.
                        AGRL(2)   = 0.0
                        AGRL(4)   = 1.0
                     ELSEIF (LEXTRP) THEN
                        AGRL(2)   = XKK-X(NIN)
                        AGRL(4)   = 1.0
                     ELSE
                        AGRL(2)   = 0.0
                        AGRL(4)   = 0.0
                     ENDIF
C-----------------------------------------------------------------------
C EXTRAPOLATE: LOW 'XOUT' VALUE  - IF EXTRAPOLATION SWITCHED ON
C-----------------------------------------------------------------------
               ELSEIF (XK.LE.XIN(1))   THEN
                  INTER     = 2
                  LINTRP(K) = .FALSE.
                  AGRL(2)   = 0.0
                  AGRL(4)   = 0.0
                     IF (XK.EQ.XIN(1)) THEN
                        LINTRP(K) = .TRUE.
                        AGRL(1)   = 0.0
                        AGRL(3)   = 1.0
                     ELSEIF (LEXTRP) THEN
                        AGRL(1)   = XKK-X(1)
                        AGRL(3)   = 1.0
                     ELSE
                        AGRL(1)   = 0.0
                        AGRL(3)   = 0.0
                     ENDIF
C-----------------------------------------------------------------------
C INTERPOLATE:
C-----------------------------------------------------------------------
               ELSE
                     DO 6 I=NIN,1,-1
                        IF (XK.LT.XIN(I)) INTER=I
    6                CONTINUE
                  LINTRP(K) = .TRUE.
                  DL1       =  ( X(INTER) - XKK    ) / H(INTER)
                  DL2       =  1.0 - DL1
                  DL3       =  H(INTER) * DL1 * DL2
                  AGRL(1)   =  DL1 * DL3
                  AGRL(2)   = -DL2 * DL3
                  AGRL(3)   =  DL1 * DL1 * ( 1.0 + DL2 + DL2 )
                  AGRL(4)   =  1.0 - AGRL(3)
               ENDIF
C-----------------------------------------------------------------------
C EVALUATE 'YOUT'-VALUE FOR REQUESTED 'XOUT'-VALUE
C-----------------------------------------------------------------------
            YOUT(K) = ( AGRL(1)* DY(INTER-1) ) + ( AGRL(2)* DY(INTER) )
     &              + ( AGRL(3)*YIN(INTER-1) ) + ( AGRL(4)*YIN(INTER) )
C-----------------------------------------------------------------------
    5    CONTINUE
C***********************************************************************
      RETURN
      END
