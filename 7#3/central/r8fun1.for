CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/maths/r8fun1.for,v 1.3 2007/04/11 13:01:53 allan Exp $ Data $Date: 2007/04/11 13:01:53 $
CX
      FUNCTION R8FUN1 ( Z )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  **************** FORTRAN77 REAL*8 FUNCTION: R8FUN1 ******************
C
C  PURPOSE: Returns argument
C
C  CALLING PROGRAM: GENERAL USE
C
C  FUNCTION:
C          (R*8)  R8FUN1  = FUNCTION NAME
C          (R*8)  Z       = INPUT VALUE
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     13/08/90
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
      REAL*8 R8FUN1 , Z
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      R8FUN1=Z
C-----------------------------------------------------------------------
      RETURN
      END
