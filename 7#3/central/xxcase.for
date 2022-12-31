      subroutine xxcase(input,output,type)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXCASE *********************
C
C  PURPOSE: Change a string of arbitrary size into all upper case
C           or all lower case
C
C  CALLING PROGRAM: GENERAL USE.
C
C  INPUT    : (C*(*)) INPUT = Input String
C  INPUT    : (C*2)   TYPE = Type of case to convert to:
C                       'UC' -> Convert to Upper Case
C                       'LC' -> Convert to Lower Case
C                       'CA' -> Convert case of first character to 
C                               upper case
C                       Anything else -> No conversion
C
C  OUTPUT   : (C*(*)) OUTPUT = Output string in selected case
C
C  ROUTINES : NONE
C
C  AUTHOR   : Allan Whiteford,
C             University of Strathclyde
C
C  VERSION  : 1.1                          
C  DATE     : 05/09/2001
C  MODIFIED : Allan Whiteford
C             First version.
C
C  VERSION  : 1.2                          
C  DATE     : 05/05/2005
C  MODIFIED : Martin O'Mullane
C             The routine converted length-1 rather than the whole
C             input string.
C
C  VERSION  : 1.3                          
C  DATE     : 08-06-2021
C  MODIFIED : Martin O'Mullane
C              - Add capitalization of first character option (CA).
C              - Ignore the case of type setting.
C
C-----------------------------------------------------------------------
       integer i
       integer size
C----------------------------------------------------------------------
       character*(*) input
       character*(*) output
       character*2 type
C----------------------------------------------------------------------

       size=len(input)

       write(output(1:size),'(A)') input(1:size)
       i=1

       if (type.eq.'UC' .or. type.eq.'uc'.or.
     &     type.eq.'uC'.or.type.eq.'Uc') then
10        if    ( ichar(input(i:i)) .ge. ichar('a')
     &    .and.   ichar(input(i:i)) .le. ichar('z')
     &          ) output(i:i)=char(ichar(input(i:i))+ichar('A')
     &                                              -ichar('a'))

          i=i+1
          if (i .le. size) goto 10
       endif

       if (type.eq.'LC' .or. type.eq.'lc'.or.
     &     type.eq.'lC'.or.type.eq.'Lc') then
20        if    ( ichar(input(i:i)) .ge. ichar('A')
     &    .and.   ichar(input(i:i)) .le. ichar('Z')
     &          ) output(i:i)=char(ichar(input(i:i))+ichar('a')
     &                                              -ichar('A'))

          i=i+1
          if (i .le. size) goto 20
       endif

       if (type.eq.'CA'.or.type.eq.'ca'.or.
     &     type.eq.'cA'.or.type.eq.'Ca') then
          if    ( ichar(input(1:1)) .ge. ichar('a')
     &    .and.   ichar(input(1:1)) .le. ichar('z')
     &          ) output(1:1)=char(ichar(input(i:i))+ichar('A')
     &                                              -ichar('a'))
       endif

       end
