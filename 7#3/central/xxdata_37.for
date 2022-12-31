       subroutine xxdata_37( iunit  ,
     &                       nemax  , ntmax  ,
     &                       title  , icateg , nenerg , nblock ,
     &                       nform1 , param1 , nform2 , param2 ,
     &                       ea     , fa     , teff   , mode   ,
     &                       median , filnam , filout , calgeb ,
     &                       ealgeb
     &                      )

       implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 subroutine: xxdata_37 ******************
C
C  purpose:  To fetch data from an adf37 data set and detect its main
C            characteristics.
C
C  calling program: various
C
C  input : (i*4)  iunit    = unit to which input file is allocated
C  input : (i*4)  nemax    = max no of energy points that can be read in
C  input : (i*4)  ntmax    = max no of effective temps that can be read in
C
C  output: (c*80) title    = header for file
C  output: (i*4)  icateg   = category of file
C                              1 => superposition
C                              2 => numerical
C  output: (i*4)  nenerg   = type 1 => number of distribution families
C                            type 2 => number of energy points
C  output: (i*4)  nblock   = type 1 => number of members in output family
C                            type 2 => number of effective temperatures
C  output: (i*4)  nform1   = type of threshold behaviour
C                              1 => cutoff
C                              2 => energy^param1
C  output: (r*8)  param1   = parameter of threshold form
C  output: (i*4)  nform2   = type of high-energy behaviour
C                              1 => cutoff
C                              2 => energy^-param2(1)
C                              3 => exp(-param2(1)*energy)
C                              4 => exp(-param2(1)*energy^param2(2))
C  output: (r*8)  param2() = parameter of high-energy form
C  output: (r*8)  ea(,)    = energy points of tabulation
C  output: (r*8)  fa(,)    = distribution function tabulation
C  output: (r*8)  teff()   = effective temperature (eV)
C  output: (r*8)  mode()   = most probable energy (eV)
C  output: (r*8)  median() = median energy (eV)
C  output: (c*120)filnam() = file names of input families
C  output: (c*120)filout   = file name of output family
C  output: (c*25) calgeb(,)= distribution function algebra
C  output: (c*25) ealgeb() = energy parameter algebra
C
C  local : (i*4)  ieunit   = energy units of distribution function
C                              1 => kelvin
C                              2 => eV
C  local : (i*4)  i        = general use
C  local : (i*4)  j        = general use
C  local : (i*4)  med_index= energy index of median
C  local : (i*4)  mode_index()  = energy index of mode
C  local : (i*4)  dummy    = general use
C  local : (i*4)  ie       = general use
C  local : (i*4)  iblock   = general use
C  local : (r*8)  sum      = average energy contribution from i -> i+1
C  local : (r*8)  contrib()= average energy contribution from i -> i+1
C  local : (r*8)  de       = energy step from i -> i+1
C  local : (i*4)  ifirst   = position of first non-blank character in string
C  local : (i*4)  ilast    = position of last non-blank character in string
C  local : (i*4)  indx()   = index of algebra
C  local : (c*80) blank    = dummy string
C
C  routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          i4unit     ADAS      fetch unit number for output of messages
C          xxslen     ADAS      finds string length excluding leading and 
C                               trailing blanks
C
C  author: Paul Bryans, University of Strathclyde
C
C  date:   04/02/04
C
C  update:
C
C  version  : 1.1
C  date     : 04-02-2004
C  modified : Paul Bryans
C               - first version
C
C  version  : 1.2
C  date     : 06-09-2010
C  modified : Martin O'Mullane
C               - adjust reading of header to match adf37 definition.
C
C  version  : 1.3
C  date     : 27-07-2018
C  modified : Martin O'Mullane
C               - check that the size of data in the file does not 
C                 exceed the dimensions.
C
C-----------------------------------------------------------------------     
      integer   icateg           , ieunit          , nenerg  , nblock  ,
     &          nform1           , nform2          , iunit   , i4unit  ,
     &          dummy            , ntmax           , nemax   , i       ,
     &          med_index        , indx(ntmax)     , ie      , iblock  ,
     &          ifirst           , ilast           , j       ,
     &          mode_index(ntmax)
C-----------------------------------------------------------------------
      real*8    param1           , param2(2)       , sum     , de      ,
     &          ea(ntmax,nemax)  , fa(ntmax,nemax) , mode(ntmax)       ,
     &          median(ntmax)    , teff(ntmax)     , contrib(nemax)
C-----------------------------------------------------------------------     
      character title*80         , blank*80        , filout*120        ,
     &          ealgeb(ntmax)*25 , calgeb(ntmax,4)*25                  ,
     &          filnam(nemax)*120
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  read title and primary file specifications
C-----------------------------------------------------------------------
      read(iunit,1000) title
      read(iunit,*) icateg    
      read(iunit,*) ieunit
      read(iunit,*) nenerg
      read(iunit,*) nblock      

C-----------------------------------------------------------------------
C  check data in file does not exceed dimensions
C-----------------------------------------------------------------------

      if (nenerg.gt.nemax) then
        write(i4unit(-1),1002)'Too many energies in adf37 file'
        write(i4unit(-1),'(a,i5)')' Maximum number accepted : ', nemax
        write(i4unit(-1),'(a,i5)')' Number in file          : ', nenerg
        write(i4unit(-1),1003)
        stop
      endif

      if (nblock.gt.ntmax) then
        write(i4unit(-1),1002)'Too many blocks in adf37 file'
        write(i4unit(-1),'(a,i5)')' Maximum number accepted : ', ntmax
        write(i4unit(-1),'(a,i5)')' Number in file          : ', nblock
        write(i4unit(-1),1003)
        stop
      endif


C-----------------------------------------------------------------------      
C  type 1 file: superposition
C-----------------------------------------------------------------------
      if (icateg.eq.1) then
      
C-----------------------------------------------------------------------
C  read file names of input families
C-----------------------------------------------------------------------
        read(iunit,1000) blank, blank, blank
        do 80 i=1,nenerg
          read(iunit,1008) filnam(i)
          call xxslen(filnam(i),ifirst,ilast)
          filnam(i)=filnam(i)(ifirst:ilast)
   80   enddo
        
C-----------------------------------------------------------------------
C  read name of output file, remove trailing and leading blanks
C-----------------------------------------------------------------------
        read(iunit,1000) blank, blank, blank
        read(iunit,1008) filout
        call xxslen(filout,ifirst,ilast)
        filout=filout(ifirst:ilast)

C-----------------------------------------------------------------------
C  read algebra, remove trailing and leading blanks
C-----------------------------------------------------------------------
        read(iunit,1000) blank, blank, blank, blank, blank
        do 60 i=1,nblock
          read(iunit,1009) indx(i), (calgeb(i,j), j=1,nenerg), ealgeb(i)
          do 70 j=1,nenerg
            call xxslen(calgeb(i,j),ifirst,ilast)
            calgeb(i,j)=calgeb(i,j)(ifirst:ilast)
   70     enddo
          call xxslen(ealgeb(i),ifirst,ilast)
          ealgeb(i)=ealgeb(i)(ifirst:ilast)
   60   enddo

C-----------------------------------------------------------------------
C  type 2 file: numerical
C-----------------------------------------------------------------------
      elseif (icateg.eq.2) then
        
C-----------------------------------------------------------------------
C  read threshold behaviour parameters
C-----------------------------------------------------------------------        
        read(iunit,1001) nform1
        if (nform1.eq.2) then
          backspace(iunit)
          read(iunit,1005) dummy, param1          
        elseif (nform1.ne.1) then
          write(i4unit(-1),1002)'(nform1 .ne. 1 or 2) -'
          write(i4unit(-1),1004)
     &    'file contains invalid threshold behaviour ',
     &    '(must equal 1 or 2)'
          write(i4unit(-1),1003)
          stop
        endif
        
C-----------------------------------------------------------------------
C  read high energy behaviour parameters
C-----------------------------------------------------------------------
        read(iunit,1001) nform2
        if ((nform2.eq.2).or.(nform2.eq.3)) then
          backspace(iunit)
          read(iunit,1005) dummy, param2(1)
        elseif (nform2.eq.4) then
          backspace(iunit)
          read(iunit,1006) dummy, param2(1), param2(2)
        elseif (nform2.ne.1) then
          write(i4unit(-1),1002)'(nform2 .ne. 1, 2, 3 or 4) -'
          write(i4unit(-1),1004)
     &    'file contains invalid high energy behaviour ',
     &    '(must equal 1, 2, 3 or 4)'
          write(i4unit(-1),1003)
          stop
        endif
        
C-----------------------------------------------------------------------
C  read energy points and distribution
C-----------------------------------------------------------------------
        do 10 iblock=1,nblock
                  
          read(iunit,1007) (ea(iblock,i), i=1,nenerg)
          read(iunit,1007) (fa(iblock,i), i=1,nenerg)
                  
C-----------------------------------------------------------------------
C  convert to eV if units used are kelvin
C----------------------------------------------------------------------
          if (ieunit.eq.1) then
            do 20 ie=1,nenerg
              ea(iblock,ie)=ea(iblock,ie)/11605.4d0
              fa(iblock,ie)=fa(iblock,ie)*11605.4d0
   20       enddo
          elseif (ieunit.ne.2) then
            write(i4unit(-1),1002)'(ieunit .ne. 1 or 2) -'
            write(i4unit(-1),1004)
     &      'file contains invalid energy unit specifier ',
     &      '(must equal 1 or 2)'
            write(i4unit(-1),1003)
            stop
          endif
          
C-----------------------------------------------------------------------
C  calculate effective temperature and most probable energy
C----------------------------------------------------------------------
          mode(iblock) = fa(iblock,1)
          mode_index(iblock) = 1
          sum = 0.0d0
          if (nform1.eq.2)
     &      sum = fa(iblock,1)*ea(iblock,1)**2/(param1+2.0d0)
          contrib(1) = sum
          do 30 ie=1,nenerg-1
            de = ea(iblock,ie+1)-ea(iblock,ie)
            sum = sum + 0.5d0*de*( fa(iblock,ie)*ea(iblock,ie) +
     &                             fa(iblock,ie+1)*ea(iblock,ie+1) )
            contrib(ie+1) = sum
            if (fa(iblock,ie+1).gt.mode(iblock)) then
              mode(iblock) = fa(iblock,ie+1)
              mode_index(iblock) = ie+1
            endif
   30     enddo  
          if (nform2.eq.2) sum = sum + fa(iblock,nenerg)*
     &      ea(iblock,nenerg)**2/(param2(1)-2.0d0)
          if (nform2.eq.3) sum = sum + fa(iblock,nenerg)*
     &      (ea(iblock,nenerg)+1.0d0/param2(1))/param2(1)
c          if (nform2.eq.4) sum = sum + fa(iblock,nenerg)*
c     &     dexp(param2(1)*ea(iblock,nenerg)**param2(2)+
c     &      lngama(2.d0/param2(2)))/
c     &      (param2(2)*param2(1)**(2.d0/param2(2)))*
c     &      (1.d0-ingama(2.d0/param2(2),
c     &      param2(1)*ea(iblock,nenerg)**param2(2)))
          teff(iblock) = 2.0d0*sum/3.0d0

C-----------------------------------------------------------------------
C  calculate median effective temperature
C-----------------------------------------------------------------------

          do 40 ie=1,nenerg
            if (contrib(ie).ge.(sum/2.0d0)) then
              med_index = ie
              goto 50
            endif
   40     enddo
   50     median(iblock) = (sum/2.0d0-contrib(med_index-1)) /
     &                     (contrib(med_index)-contrib(med_index-1)) *
     &                     (ea(iblock,med_index)-ea(iblock,med_index-1))
     &                     + ea(iblock,med_index-1)
C-----------------------------------------------------------------------
   10   enddo

C-----------------------------------------------------------------------
C  icateg not equal to 1 or 2     
C-----------------------------------------------------------------------
      else
        write(i4unit(-1),1004)
     &  'file contains invalid category of file (must equal 1 or 2)'
        write(i4unit(-1),1003)
        stop
      endif


C-----------------------------------------------------------------------

 1000 format(1a80)
 1001 format(i5)
 1002 format(1x,31('*'),' xxdata_37 error ',30('*')//
     &       1x,'Fault in input data file: ',a)
 1003 format(/1x,29('*'),' program terminated ',29('*'))
 1004 format(1x,a)
 1005 format(i5,1d10.2)
 1006 format(i5,2d10.2)
 1007 format(7d10.2)
 1008 format(1a120)
 1009 format(i5,5a25)
 
C-----------------------------------------------------------------------
      return
      end
