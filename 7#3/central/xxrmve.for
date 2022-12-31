      subroutine xxrmve( cstrg1 , cstrg2 , crmve )
      implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 subroutine: xxrmve *********************
C
C  purpose: to remove all occurrences of a selected character from a
C           string and concatenate.  Output string tail is blank filled
C
C  calling program: general use
C
C  subroutine:
C
C  input : (c*(*)) cstrg1   = input string for conversion
C  input : (c*1)   crmve    = character to be removed
C
C  output: (c*(*)) cstrg2   = output string after conversion
C
C          (i*4)   i        = general use
C          (i*4)   ilen     = length of 'cstrng' string in bytes
C
C routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          i4unit     adas      fetch unit number for output of messages
C
C
C
C author:  H. P. Summers, university of strathclyde
C          ja7.08
C          tel. 0141-548-4196
C
C date:    06/09/01
C
C version : 1.1
C date    : 06/09/2001
C modified: Hugh Summers
C               - first edition.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      integer     ilen1       , ilen2      , i       , ic
      integer     i4unit
C-----------------------------------------------------------------------
      character   cstrg1*(*)  , cstrg2*(*)  , crmve*1
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      ilen1   = len(cstrg1)
      ilen2   = len(cstrg2)

      ic = 1
      cstrg2 = ' '

         do i=1,ilen1
            if (cstrg1(i:i).ne.crmve) then
                if(ic.le.ilen2) then
                   cstrg2(ic:ic)=cstrg1(i:i)
                   ic=ic+1
                else
                   write(i4unit(-1),1001)
                   write(i4unit(-1),1002)
                endif
            endif

          enddo

       return

C-----------------------------------------------------------------------
 1001 format(1x,32('*'),' xxrmve warning ',32('*')//
     &       1x,'output string too short and returned incomplete: ')
 1002 format(/1x,29('*'),' subroutine terminated ',29('*'))
C-----------------------------------------------------------------------
       end
