CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/utility/lenstr.for,v 1.2 2007/04/11 13:02:01 allan Exp $ Date $Date: 2007/04/11 13:02:01 $
CX
      FUNCTION LENSTR(ASTR)
C
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 FUNCTION: LENSTR ***********************
C
C     PURPOSE : RETURNS THE EFFECTIVE LENGTH OF A GIVEN STRING
C               (IGNORING TRAILING BLANKS)
C
C NOTES: THIS ROUTINE IS NOT YET PROPERLY ANNOTATED
C
C UNIX-IDL PORT:
C
C VERSION: 1.1                          DATE: 18-1-96
C MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
C               - PUT UNDER SCCS CONTROL
C
C VERSION  : 1.2                          
C DATE     : 10-04-2007
C MODIFIED : Allan Whiteford
C               - Modified documentation as part of automated
C		  subroutine documentation preparation.
C-----------------------------------------------------------------------
C
      CHARACTER*(*) ASTR
C
      DO 10 I = LEN(ASTR),1,-1
        IF (ASTR(I:I) .NE. ' ') GO TO 20
   10 CONTINUE
C
   20 CONTINUE
      LENSTR = I
C
      RETURN
      END
