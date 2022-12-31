      subroutine xxcomm(iunit, ndcomm, comments, ncomm)

      implicit none
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXCOMM *********************
C
C  PURPOSE: Read the comments from an ADAS adf dataset.
C
C  CALLING PROGRAM: General use.
C
C  INPUT:  (C*(*)) dsname    = dataset name
C  INPUT:  (I*4)   ndcomm    = maximum number of comment lines
C
C  OUTPUT: (I*4)   ncomm     = nunber of comment lines
C  OUTPUT: (C*(*)) comments  = array of comment strings
C
C  ROUTINES : NONE
C
C  AUTHOR   : Martin O'Mullane,
C
C  VERSION  : 1.1
C  DATE     : 14-04-2005
C  MODIFIED : Martin O'Mullane
C              - First version.
C
C  VERSION  : 1.2
C  DATE     : 28-10-2010
C  MODIFIED : Martin O'Mullane
C              - Permit longer strings - but only fill up to length of
C                the input comments array.
C
C  VERSION  : 1.3
C  DATE     : 03-12-2012
C  MODIFIED : Martin O'Mullane
C              - Clarify warning on number of comment lines.
C
C-----------------------------------------------------------------------
      integer   ndcomm     , ncomm , iunit , i4unit  , 
     &          ilen       , L1    , L2
C----------------------------------------------------------------------
      character str*180
C----------------------------------------------------------------------
      character comments(ndcomm)*(*)
C----------------------------------------------------------------------

C Only fill length of comment array

      ilen = len(comments(1))

C Skip the first line to avoid carbon pecs etc.

      read(iunit,'(A)')str

      ncomm = 1
  10  read(iunit,'(A)', end=20)str
      if (str(1:1).eq.'C'.or.str(1:1).eq.'c') then
         call xxslen(str, L1, L2)
         if (L2.GT.ilen) L2 = ilen
         comments(ncomm) = str(1:L2)
         ncomm = ncomm + 1
         if (ncomm.GT.ndcomm) then
            write(i4unit(-1),*)
     &          'Too many comment lines - stop reading at',ncomm-1
            ncomm = ndcomm
            return
         endif
       endif
       goto 10
  20   continue
       ncomm = ncomm - 1

       end
