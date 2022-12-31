CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/maths/eei.for,v 1.3 2007/04/23 13:40:45 allan Exp $ Date $Date: 2007/04/23 13:40:45 $
CX
      FUNCTION EEI(X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 FUNCTION: EEI **************************
C
C PURPOSE: EVALUATES EXP(X)E1(X) WHERE E1 IS THE 1ST EXPONENTIAL
C          INTEGRAL
C
C CALLING PROGRAMS: GENERAL
C
C INPUT:  (R*8)  X       = INDEPENDENT VARIABLE
C
C OUTPUT: (R*8)  EEI     = EXP(X)E1(X) 
C
C ROUTINES: NONE
C
C UNIX-IDL PORT:
C
C VERSION: 1.1                          DATE: 11-7-95
C MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
C               - FIRST VERSION
C
C VERSION: 1.2                          DATE: 16-1-96
C MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
C               - TIDIED UP COMMENTS
C
C VERSION: 1.3                          DATE: 17-4-07
C MODIFIED: HUGH SUMMERS
C           - COMPLETED COMMENT BLOCK DESCRIPTION
C
C-----------------------------------------------------------------------
C
      IF (X.LE.1.0D0) THEN
C
         A = -LOG(X) - 0.57721566D0 + X * (0.99999193D0 - X *
     +       (0.24991055D0 - X * (0.05519968D0 - X *
     +       (0.00976004D0 - X * 0.00107857D0))))
         Y = 0.5D0 * X
         Z = 1.0D0 - Y * (0.9998684D0 - Y * (0.4982926D0 - Y *
     +       (0.1595332D0 - Y * 0.0293641D0)))
         EEI = A / (Z * Z)
C
      ELSE
C
         EEI = (0.2677737343D0 + X * (8.6347608925D0 + X *
     +         (18.059016973D0 + X * (8.5733287401D0 + X)))) /
     +         (X * (3.9584969228D0 + X * (21.0996530827D0 + X *
     +         (25.6329561486D0 + X * (9.5733223454D0 + X)))))
C
      ENDIF
C
      RETURN
      END
