CX UNIX PORT - SCCS Info : Module @(#)$Header: /home/adascvs/fortran/adaslib/utility/xxslen.for,v 1.1 2004/07/06 15:39:12 whitefor Exp $ Date $Date: 2004/07/06 15:39:12 $
CX      
      SUBROUTINE XXSLEN( CSTRNG , IFIRST , ILAST )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXSLEN *********************
C
C  PURPOSE: TO IDENTIFY THE FIRST AND LAST NON-BLANK CHARACTER IN A
C           STRING. (IF INPUT STRING IS BLANK IFIRST=ILAST=0)
C
C  CALLING PROGRAM: GENERAL USE
C
C  SUBROUTINE:
C
C  INPUT : (C*(*)) CSTRNG   = INPUT STRING FOR INTERROGATION
C
C  OUTPUT: (I*4)   IFIRST   = BYTE POSITION OF FIRST NON-BLANK CHARACTER
C                             IN INPUT STRING.
C  OUTPUT: (I*4)   ILAST    = BYTE POSITION OF LAST  NON-BLANK CHARACTER
C                             IN INPUT STRING.
C
C          (I*4)   I        = GENERAL USE
C          (I*4)   ILEN     = LENGTH OF 'CSTRNG' STRING IN BYTES
C
C ROUTINES: NONE
C
C NOTE:
C
C
C AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C          K1/0/37
C          JET EXT. 6023
C
C DATE  :  06/07/93
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      INTEGER     IFIRST     , ILAST    , ILEN    , I
C-----------------------------------------------------------------------
      CHARACTER   CSTRNG*(*)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
      ILEN   = LEN(CSTRNG)
C-----------------------------------------------------------------------
      IFIRST = 0
      ILAST  = 0
C-----------------------------------------------------------------------
C
         DO 1 I=1,ILEN
C
            IF (CSTRNG(I:I).NE.' ') THEN
               IF (IFIRST.EQ.0) IFIRST = I
               ILAST = I
            ENDIF
C
    1    CONTINUE
C
C-----------------------------------------------------------------------
C
       RETURN
       END
