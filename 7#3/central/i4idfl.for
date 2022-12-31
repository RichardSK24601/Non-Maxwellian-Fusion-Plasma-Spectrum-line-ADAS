CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/atomic/i4idfl.for,v 1.3 2007/04/11 13:01:47 allan Exp $ Data $Date: 2007/04/11 13:01:47 $
CX
      FUNCTION I4IDFL( N , L )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  *************** FORTRAN77 INTEGER*4 FUNCTION: I4INDL ****************
C
C  PURPOSE:  RETURNS A UNIQUE INDEX NUMBER BASED ON THE VALUE OF THE
C            N AND L QUANTUM NUMBERS PASSED TO IT. THE INDEX IS USED TO
C            REFERENCE ARRAYS CONTAINING DATA DEPENDENT ON THE N AND L
C            QUANTUM NUMBERS.
C
C  CALLING PROGRAM: GENERAL USE
C
C  FUNCTION:
C
C  FUNC:   (I*4)   I4IDFL  = INDEX
C
C  INPUT:  (I*4)   N       = N QUANTUM NUMBER.
C  INPUT:  (I*4)   L       = L QUANTUM NUMBER.
C
C AUTHOR:   JONATHAN NASH (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 5183
C
C DATE:     10/09/93
C
C VERSION  : 1.2                          
C DATE     : 20-12-2001
C MODIFIED : Martin O'Mullane
C               - Removed mainframe listing information beyond column 72.
C
C VERSION  : 1.3                          
C DATE     : 10-04-2007
C MODIFIED : Allan Whiteford
C               - Modified documentation as part of automated
C		  subroutine documentation preparation.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER    I4IDFL
C-----------------------------------------------------------------------
      INTEGER    N  , L
C-----------------------------------------------------------------------
C
C***********************************************************************
C CALCULATE INDEX.
C***********************************************************************
C
      I4IDFL = ((N * (N - 1)) / 2) + L + 1
C
C-----------------------------------------------------------------------
C
      RETURN
      END
