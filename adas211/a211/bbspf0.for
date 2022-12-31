      SUBROUTINE BBSPF0( REP , DSFULL , DIST , DPARAM , DS37)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: BBSPF0 *********************
C
C  PURPOSE: TO DISPLAY AND FETCH VALUES FROM IDL OR INFO FILE (STANDARD
C	    (INPUT) IF BATCH - INPUT DATA SET SPECIFICATIONS.
C
C  CALLING PROGRAM: ADAS211
C
C  SUBROUTINE:
C
C  OUTPUT: (C*3)   REP     = 'YES' => TERMINATE PROGRAM EXECUTION.
C                          = 'NO ' => CONTINUE PROGRAM EXECUTION.
C
C  OUTPUT: (C*80)  DSFULL  = INPUT DATA SET NAME (FULL MVS DSN)
C                            (IN FORM SUITABLE FOR DYNAMIC ALLOCATION)
C
C  OUTPUT: (I*4)   DIST	   = The distribution type from the widget
C	                        0 - Maxwellian
C	                        1 - Kappa
C	                        2 - Numerical
C	                        3 - Druvesteyn
C
C  OUTPUT: (R*8)   DPARAM  = Non-Maxwellian parameter:
C	                         Maxwellian: N/A
C	                         Kappa:      Kappa
C	                         Numerical:  N/A
C	                         Druvesteyn: x
C
C  OUTPUT: (C*80)  DS37    = Numerical distribution file (only
C                            applicable for DIST=2)
C
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C
C AUTHOR:  H. P. SUMMERS, UNIVERSITY OF STRATHCLYDE
C          JA 8.08
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
C VERSION: 1.2                          DATE: 02-02-05
C MODIFIED: ALLAN WHITEFORD
C               - MODIFIED TO READ DIST, DPARAM AND DS37
C                 FOR NON-MAXWELLIAN PROCESSING
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER     PIPEIN , DIST
      REAL*8      DPARAM
      PARAMETER ( PIPEIN=5          )
C
      CHARACTER   REP*3         , DSFULL*80 , DS37*80
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C READ IN DATA FROM STANDARD INPUT (IDL OR INFO FILE)
C-----------------------------------------------------------------------
C

      READ(PIPEIN,'(A)')REP
      READ(PIPEIN,'(A)')DSFULL
      READ(PIPEIN,*) DIST
      READ(PIPEIN,*) DPARAM
      READ(PIPEIN,'(A)')DS37

      RETURN
      END

