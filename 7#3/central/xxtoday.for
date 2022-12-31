      subroutine xxtoday(str_today)

      implicit none
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: xxtoday ********************
C
C PURPOSE: To return today's date.
C
C CALLING PROGRAM: General use
C
C OUTPUT: (C*10)  str_today    = Today's data in DD-MM-YYYY format.
C
C ROUTINES:
C          ROUTINE        SOURCE
C          ------------------------------------------------------------
C          date_and_time  intrinsic fortran routine
C
C
C AUTHOR   : Martin O'Mullane
C DATE     : 19-04-2005
C
C
C VERSION  : 1.1
C DATE     : 19-04-2005
C MODIFIED : Martin O'Mullane
C               - First version.
C
C VERSION  : 1.2
C DATE     : 23-05-2019
C MODIFIED : Martin O'Mullane
C               - Use the intrinsic date_and_time routine rather
C                 than linking against the C routine in today.c
C
C-----------------------------------------------------------------------
      integer       values(8)
C-----------------------------------------------------------------------
      character     date*8, time*10, zone*5, str_today*10
C-----------------------------------------------------------------------

      call date_and_time(date, time, zone, values)
      write(str_today, 100)values(3), values(2), values(1)

  100 format(i2.2, '-', i2.2, '-', i4.4)

      end
