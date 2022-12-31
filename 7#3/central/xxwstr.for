C UNIX-IDL PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/utility/xxwstr.for,v 1.1 2004/07/06 15:40:46 whitefor Exp $ Date $Date: 2004/07/06 15:40:46 $
C
      subroutine xxwstr(iunit,string)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXWSTR *********************
C
C  PURPOSE: Writes a string with format (a) to iunit with no 
C           trailing blanks.
C
C  CALLING PROGRAM: GENERAL USE.
C
C  INPUT    : (C*(*)) STRING = STRING
C
C  INPUT    : (I*4)   IUNIT  = OUTPUT UNIT
C
C  ROUTINES :
C            ROUTINE    SOURCE    BRIEF DESCRIPTION
C            ------------------------------------------------------------
C            XXSLEN     ADAS      GETS NON-BLANK LENGTH OF A STRING
C
C  AUTHOR   : Martin O'Mullane,
C             K1/1/43,
C             JET
C
C  VERSION  : 1.1                          DATE: 17/03/1999
C  MODIFIED : Martin O'Mullane
C             First version.
C
C-----------------------------------------------------------------------
      INTEGER    IUNIT  , L1, L2
C----------------------------------------------------------------------
      CHARACTER  STRING*(*) 
C-----------------------------------------------------------------------
                  
      call xxslen(string,L1,L2)
      write(iunit,'(A)')string(1:L2)
 
      end
