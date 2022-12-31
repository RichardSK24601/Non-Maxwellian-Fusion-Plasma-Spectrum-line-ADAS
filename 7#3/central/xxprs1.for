       subroutine xxprs1(ndmet,string,wno,cpl,npt,ipla,zpla,ifail)

       implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 subroutine: xxprs1 *********************
C
C  purpose:  to analyse the tail character string of an level data line
C            of an adf04 specific ion file into wave-number and sets of
C            (parent identifier, effective zeta for the parent) pairs.
C
C            unified version of baprs1, b9prs1, bbprs1, g5prs1 which is
C            a replacement for these subroutines
C
C  calling program: various
C
C  notes: detect  -  level wave number which preceeds first '{'
C                 -  sets of   parent index contained in '{.}'
C                              followed by effective zeta
C         nb. 'x' as first parent assignment means exclude ionisation
C             from this level.
C             no parent assignment means take lowest parent with
C             zeta =1.
C             lowest parent but no zeta means take zeta =1.
C             if there is more than one parent then zeta's must be in.
C
C
C  subroutine:
C
C  input : (i*4)  ndmet    =  maximum number of parents
C  input : (c*(*))string   =  string to be parsed
C
C  output: (r*8)  wno      =  excitation wave number of level relative
C                             to lowest parent
C  output: (c*1)  cpl      =  lead parent for ionisation  or 'x'
C  output: (i*4)  npt      =  number of parents detected
C  output: (i*4)  ipla()   =  parent indices.
C  output: (r*8)  zpla()   =  effective zeta for parent ipla()
C  output: (i*4)  ifail    =  0 - subroutine concludes correctly
C                             1 - fault detected in subroutine
C                             2 - single ionisation potential detected
C
C          (i*4)  maxwrd   =  maximum number of words sought initially
C                             initially, finally number actually found
C          (i*4)  nfirst   =  first word to be extracted from string
C          (i*4)  ifirst() =  index of first char. of word () in string
C          (i*4)  ilast()  =  index of last  char. of word () in string
C          (i*4)  iwords   =  number of words found in string
C
C          (i*4)  ic       =  general use
C          (i*4)  iabt     =  failure number from r8fctn
C          (c*15) sstrng   =  isolated substring
C
C routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          i4unit     adas      fetch unit number for output of messages
C          r8fctn     adas      converts from character to real variable
C          i4fctn     adas      converts from char. to integer  variable
C          xxword     adas      parses a string into separate words
C                               for ' ()<>{}' delimiters
C
C AUTHOR:  HP Summers
C          JA7.08, University of Strathclyde
C          Tel: 0141-548-4196
C
C DATE:    04/12/02
C
C
C VERSION : 1.2
C DATE    : 09-04-2010
C MODIFIED: Martin O'Mullane
C           - Change integer*4 to integer.
C
C VERSION : 1.3
C DATE    : 27-10-2016
C MODIFIED: Martin O'Mullane
C           - Set idwds as internal parameter and check that maxwrd
C             is smaller than it, otherwise stop.
C           - Increase idwds from 12 to 30.
C
C VERSION : 1.4
C DATE    : 19-01-2019
C MODIFIED: Martin O'Mullane
C           - The specialized adf04 files holding NIST data have
C             extra information on the level given as a comment
C             delineated by a '#' after the ionization parentage.
C             Exclude this before parsing.
C
C-----------------------------------------------------------------------
       integer   idwds
C-----------------------------------------------------------------------
       parameter(idwds = 30 )
C-----------------------------------------------------------------------
       integer   ndmet
       integer   npt          , iabt       , ic        , i
       integer   ifail
       integer   nfirst       , maxwrd     , iwords
       integer   i4fctn       , i4unit
       integer   ihash
C-----------------------------------------------------------------------
       real*8    wno
       real*8    r8fctn
C-----------------------------------------------------------------------
       character string*(*)   , cpl*1
       character cdelim*7
C-----------------------------------------------------------------------
       integer   ipla(ndmet)
       integer   ifirst(idwds)   , ilast(idwds)
C-----------------------------------------------------------------------
       real*8    zpla(ndmet)
C-----------------------------------------------------------------------
       data cdelim/' ()<>{}'/
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C Set and check dimensions are sufficient
C-----------------------------------------------------------------------

       maxwrd=2*ndmet+1

       if (maxwrd.gt.idwds)then
           write(i4unit(-1),1000)'increase number allowed words ',
     &           maxwrd, idwds
           write(i4unit(-1),1002)
           stop
       endif

C-----------------------------------------------------------------------

       npt    = 0
       ifail  = 0
       wno    = 0.0D0

       do ic=1,ndmet
         ipla(ic)=0
         zpla(ic)=0.0d0
       end do

C Remove any text from the first # onwards

       ihash = index(string, '#')
       if (ihash.GT.0) then
          do i = ihash, len(string)
             string(i:i) = ' '
          end do
       endif

       nfirst=1
       call xxword(string,cdelim,nfirst,maxwrd,ifirst,ilast,iwords)
       if(iwords.eq.0)then
           write(i4unit(-1),1001)'no excitation energy'
           write(i4unit(-1),1002)
           ifail = 1
           return
       elseif (iwords.eq.1) then
           wno=r8fctn(string(ifirst(1):ilast(1)),iabt)
           if(iabt.gt.0)then
               write(i4unit(-1),1001)'fault in excit. energy'
               write(i4unit(-1),1002)
               ifail = 1
               return
           else
               npt=0
               cpl='x'
               ipla(1)=0
               zpla(1)=0.0D0
               ifail = 2
               return
           endif
       elseif (iwords.eq.2) then
           if((string(ifirst(2):ilast(2)).eq.'x').or.
     &        (string(ifirst(2):ilast(2)).eq.'X'))  then
               wno=r8fctn(string(ifirst(1):ilast(1)),iabt)
               if(iabt.gt.0)then
                   write(i4unit(-1),1001)'fault in excit. energy'
                   write(i4unit(-1),1002)
                   ifail = 1
                   return
               endif
               npt=1
               cpl='x'
               ipla(1)=1
               zpla(1)=0.0D0
               ifail = 2
               return
           else
               wno=r8fctn(string(ifirst(1):ilast(1)),iabt)
               if(iabt.gt.0)then
                   write(i4unit(-1),1001)'fault in excit. energy'
                   write(i4unit(-1),1002)
                   ifail = 1
                   return
               endif
               npt=1
               ipla(1)=i4fctn(string(ifirst(2):ilast(2)),iabt)
               if(iabt.gt.0.or.ipla(1).gt.ndmet)then
                   write(i4unit(-1),1001)'fault in parent index '
                   write(i4unit(-1),1002)
                   ifail = 1
                   return
               endif
               cpl=string(ifirst(2):ilast(2))
               zpla(1)=1.0D0
               ifail = 2
               return
           endif
       endif
       wno=r8fctn(string(ifirst(1):ilast(1)),iabt)
       if(iabt.gt.0)then
           write(i4unit(-1),1001)'fault in excit. energy'
           write(i4unit(-1),1002)
           ifail = 1
           return
       endif

       if (mod(iwords-1,2).ne.0)then
           write(i4unit(-1),1001)'mismatch of parents'
           write(i4unit(-1),1002)
           ifail = 1
           return
       endif

       npt=(iwords-1)/2

       do 100 i=1,npt
        ipla(i)=i4fctn(string(ifirst(2*i):ilast(2*i)),iabt)
            if(iabt.gt.0)then
                write(i4unit(-1),1001)'fault in parent index '
                write(i4unit(-1),1002)
                ifail = 1
                return
            endif
        zpla(i)=r8fctn(string(ifirst(2*i+1):ilast(2*i+1)),iabt)
            if(iabt.gt.0)then
                write(i4unit(-1),1001)'fault in zeta value   '
                write(i4unit(-1),1002)
                ifail = 1
                return
            endif
  100  continue
       cpl=string(ifirst(2):ilast(2))
       return
C-----------------------------------------------------------------------
 1000 format(1x,32('*'),' xxprs1 error ',32('*')//
     &       1x,'dimension fault: ',a,2i3)
 1001 format(1x,32('*'),' xxprs1 error ',32('*')//
     &       1x,'fault in input data file: ',a,i3,a)
 1002 format(/1x,27('*'),' subroutine terminated ',28('*'))
C-----------------------------------------------------------------------
      end
