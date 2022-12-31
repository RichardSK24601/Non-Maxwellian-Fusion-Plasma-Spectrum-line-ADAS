CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/bbspf1.for,v 1.1 2004/07/06 11:36:33 whitefor Exp $ Date $Date: 2004/07/06 11:36:33 $
CX
      SUBROUTINE BBSPF1( OPEN10 , OPEN11 ,
     &                   LCONT  , LPASS  , LPAPER , LBTSEL ,
     &                   DSNOUT , DSNPAS , DSNPAP ,
     &                   LPEND  , TITLE  , *
     &                 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: BBSPF1 *********************
C
C  PURPOSE: TO RETRIEVE OUTPUT DATA SET SPECIFICATIONS FROM IDL OR
C	    THE INFO FILE (WHICH HAS BEEN REDIRECTED TO STANDARD INPUT).
C
C  CALLING PROGRAM: ADAS211
C
C  SUBROUTINE:
C
C  I/O   : (L*4)   OPEN10   = .TRUE.  => ADF08 PASSING FILE OPENED
C                             .FALSE. => ADF08 PASSING FILE NOT OPENED
C  I/O   : (L*4)   OPEN11   = .TRUE.  => ADF04   PASSING FILE OPENED
C                             .FALSE. => ADF04   PASSING FILE NOT OPENED
C
C  OUTPUT: (L*4)   LCONT    = .TRUE.  => OUTPUT DATA TO ADF08 PASSING
C                                        FILE.
C                             .FALSE. => NO OUTPUT OF CURRENT DATA TO
C                                        ADF08 PASSING FILE.
C  OUTPUT: (L*4)   LPASS    = .TRUE.  => OUTPUT DATA TO ADF04 PASSING
C                                        FILE.
C                             .FALSE. => NO OUTPUT OF CURRENT DATA TO
C                                        ADF04 PASSING FILE.
C  OUTPUT: (L*4)   LPAPER   = .TRUE.  => OUTPUT TEXT TO 'PAPER.TXT'
C                                        FILE.
C                             .FALSE. => NO TEXT OUTPUT
C  OUTPUT: (L*4)   LBTSEL   = .TRUE.  => USER HAS SELECTED BATCH EXECUTION
C                             .FALSE. => 'RUN NOW' OR 'CANCEL' OPTION
C
C  OUTPUT: (C*80)  DSNOUT   = OUTPUT ADF08 DATA SET NAME
C  OUTPUT: (C*80)  DSNPAS   = OUTPUT PASSING FILE NAME
C  OUTPUT: (C*80)  DSNPAP   = OUTPUT TEXT FILE NAME
C
C  OUTPUT: (L*4)   LPEND    = .TRUE.  => USER SELECTED 'CANCEL'
C                             .FALSE. => USER DID NOT
C
C  OUTPUT: (C*40)  TITLE    = USER SUPPLIED TITLE FOR THIS RUN
C
C          (I*4)   PIPEIN   = STANDARD INPUT
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C
C AUTHOR:  H. P. SUMMERS, UNIVERSITY OF STRATHCLYDE
C          JA8.08
C          TEL. 0141-553-4196
C
C DATE:    21/06/96
C
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
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
      INTEGER      PIPEIN, LOGIC
      PARAMETER   (PIPEIN = 5)
C-----------------------------------------------------------------------
      CHARACTER    DSNOUT*80   , DSNPAS*80 , DSNPAP*80 ,  TITLE*40
C-----------------------------------------------------------------------
      LOGICAL      OPEN10      , OPEN11    , LBTSEL    ,  LPAPER  ,
     &             LCONT       , LPASS     , LPEND
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      READ(PIPEIN,*) LOGIC
      IF(LOGIC.NE.1)THEN
         LPEND = .FALSE.
         READ(PIPEIN,'(A)')DSNOUT
         READ(PIPEIN,'(A)')DSNPAS
         READ(PIPEIN,'(A)')DSNPAP
         READ(PIPEIN,'(A)')TITLE
         READ(PIPEIN,*)LOGIC
         IF(LOGIC.EQ.0)THEN
            LCONT=.FALSE.
         ELSE
            LCONT=.TRUE.
         ENDIF
         READ(PIPEIN,*)LOGIC
         IF(LOGIC.EQ.0)THEN
            LPASS=.FALSE.
         ELSE
            LPASS=.TRUE.
         ENDIF
         READ(PIPEIN,*)LOGIC
         IF(LOGIC.EQ.0)THEN
            LPAPER=.FALSE.
         ELSE
            LPAPER=.TRUE.
         ENDIF
         READ(PIPEIN,*)LOGIC
         IF(LOGIC.EQ.0)THEN
            LBTSEL=.FALSE.
         ELSE
            LBTSEL=.TRUE.
         ENDIF
      ELSE
         LPEND = .TRUE.
      ENDIF

      RETURN
      END
