CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/system/i4unit.for,v 1.1 2004/07/06 14:08:05 whitefor Exp $ Date $Date: 2004/07/06 14:08:05 $
CX
      FUNCTION I4UNIT( IUNIT )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ************** FORTRAN77 INTEGER*4 FUNCTION: I4UNIT *****************
C
C  PURPOSE: TO RESET OR RETURN A STORED INTEGER*4 VALUE GREATER THAN OR
C           EQUAL TO ZERO.
C           THIS IS USED WITHIN ADAS TO STORE THE STREAM/UNIT NUMBER
C           FOR THE OUTPUT OF ERROR MESSAGES (TO THE SCREEN).
C
C           BY DEFAULT THE STORED VALUE WILL BE 6, AND WILL BE RETURNED
C           BY THE FUNCTION IF IUNIT ON INPUT < 0.
C
C           TO RESET THE STORED VALUE THEN SET IUNIT TO THE REQUIRED
C           POSITIVE INTEGER (INC. ZERO). THIS VALUE WILL ALSO BE
C           RETURNED BY THE FUNCTION.
C
C                 IUNIT VALUE               RETURNED FUNCTION VALUE
C                 -----------               -----------------------
C                 IUNIT <  0            = CURRENT STORED INTEGER VALUE
C                                         (6 BY DEFAULT).
C                 IUNIT >= 0            = IUNIT , AND RESETS THE STORED
C                                                 VALUE TO IUNIT.
C
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  O     : (I*4)  I4UNIT   = FUNCTION NAME - (SEE ABOVE)
C
C  I     : (I*4)  IUNIT    = FUNCTION ARGUMENT - (SEE ABOVE)
C
C          (I*4)  IDEFLT   = PARAMETER = DEFAULT STORED INTEGER VALUE
C
C          (I*4)  ICURNT   = CURRENT STORED INTEGER VALUE
C
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 5023
C
C DATE:    23/04/93
C
C UPDATE:  24/05/93 - PE BRIDEN - ALLOWED 0 TO BE A VALID STORED NUMBER
C
C-----------------------------------------------------------------------
      INTEGER    I4UNIT     , IDEFLT
C-----------------------------------------------------------------------
      PARAMETER( IDEFLT = 6 )
C-----------------------------------------------------------------------
      INTEGER    IUNIT      , ICURNT
C-----------------------------------------------------------------------
      SAVE       ICURNT
C-----------------------------------------------------------------------
      DATA       ICURNT / IDEFLT /
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C RETRIEVE OR RESET STORED INTEGER VALUE ACCORDINGLY
C-----------------------------------------------------------------------
C
      IF (IUNIT.LT.0) THEN
        I4UNIT = ICURNT
      ELSE
        ICURNT = IUNIT
        I4UNIT = IUNIT
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
