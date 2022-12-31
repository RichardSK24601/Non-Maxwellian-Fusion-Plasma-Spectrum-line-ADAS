       subroutine xxwrto_04( iunit  ,
     &                       ndlev  , ndtrn  , nvmax , ndmet  , ndqdn ,
     &                       ndtext ,
     &                       iz     , iz0    , iz1   ,
     &                       npl    , bwnoa  , cprta ,
     &                       il     , qdorb  , lqdorb, qdn    , iorb  ,
     &                       cstrga , isa    , ila   , xja    , wa    ,
     &                       cpla   , npla   , ipla  , zpla   ,
     &                       nv     , scef   ,
     &                       itran  ,
     &                       tcode  , i1a    , i2a   , aval  , scom   ,
     &                       iadftyp,
     &                       lbeth  , beth   ,
     &                       ntext  , ctext  ,
     &                       lx_untied
     &                  )
       implicit none
c----------------------------------------------------------------------
c
c  ***************** fortran77 subroutine: xxwrto_04 ******************
c
c  purpose:  write out a valid adf04 dataset.
c
c  subroutine:
c
c  input : (i*4)  iunit    = unit to which output file is allocated
c  input : (i*4)  ndlev    = maximum number of levels
c  input : (i*4)  ndtrn    = max. number of transitions
c  input : (i*4)  nvmax    = max. number of temp/energy/partial wave-l/
c                            non-mxwl ener parameter values (according
c                            to iadftyp)
c  input : (i*4)  ndmet    = max. number of metastables
c  input : (i*4)  ndqdn    = max. number of nshells for quantum defects
c  input : (i*4)  ndtext   = max. number of comment text lines
c
c  input : (i*4)  iz       =  recombined ion charge read
c  input : (i*4)  iz0      =         nuclear charge read
c  input : (i*4)  iz1      = recombining ion charge read
c                            (note: iz1 should equal iz+1)
c
c  input : (i*4)  npl      = no. of parents on first line of adf04 file
c  input : (r*8)  bwnoa()  = ionisation potential (cm-1) of parents
c                            1st.dim.: parent index
c  input : (c*9)  cprta()  = parent name in brackets
c                            1st dim.: parent index
c
c  input : (i*4)  il       = input data file: number of energy levels
c  input : (r*8)  qdorb()  = quantum defects for orbitals
c                            1st dim: index for nl orbital (cf i4idfl.for)
c  input : (l*4)  lqdorb() = .true.  => source data available for qd.
c                          = .false. => source data not availabe qd.=0.0
c  input : (r*8)  qdn()    = quantum defect for n-shells.  non-zero only
c                            for adf04 files with orbital energy data
c                            1st. dim: n-shell (1<=n<=ndqdn)
c  input : (i*4)  iorb     = input data file: number of orbital energies
c
c  input : (c*(8))cstrga() = nomenclature/configuration for level
c                            1st dim.: level index
c  input : (i*4)  isa()    = multiplicity for level 'ia()'
c                            note: (isa-1)/2 = quantum number (s)
c                            1st dim.: level index
c  input : (i*4)  ila()    = quantum number (l) for level 'ia()'
c                            1st dim.: level index
c  input : (r*8)  xja()    = quantum number (j-value) for level 'ia()'
c                            note: (2*xja)+1 = statistical weight
c                            1st dim.: level index
c  input : (r*8)  wa()     = energy relative to level 1 (cm-1) for level
c                            1st dim.: level index
c
c  input : (c*1)  cpla()   = char. specifying 1st parent for level 'ia()'
c                            integer - parent in bwnoa() list
c                            'blank' - parent bwnoa(1)
c                              'X'   - do not assign a parent
c                            1st dim.: level index
c  input : (i*4)  npla()   = no. of parent/zeta contributions to ionis.
c                            of level
c                            1st dim.: parent index
c  input : (i*4)  ipla(,)  = parent index for contributions to ionis.
c                            of level
c                            1st dim.: parent index
c                            2nd dim.: level index
c  input : (r*8)  zpla(,   = eff. zeta param. for contributions to ionis.
c                            of level
c                            1st dim.: parent index
c
c  input : (i*4)  nv       = number of temp/energy/partial wave-l/
c                            non-mxwl ener parameter values (according
c                            to iadftyp
c  input : (r*8)  scef()   = temp/energy/partial wave-l/
c                            non-mxwl ener parameter values (according
c                            to iadftyp
c                            1st dim.: temp/energy parameter index
c  input : (i*4)  itran   = input data file: number of transitions
c
c  input : (c*1)  tcode()  = transition: data type pointer:
c                            ' ' => electron impact   transition
c                            'p' => proton   impact   transition
c                            'h' => charge   exchange recombination
c                            'r' => free     electron recombination
c                            1st dim. transition index
c  input : (i*4)  i1a()    = transition:
c                             lower energy level index (case ' ' & 'p')
c                             signed parent index (case 'h','r','s' & 'i')
c                            1st dim. transition index
c  input : (i*4)  i2a()    = transition:
c                             upper energy level index (case ' ' & 'p')
c                             capturing (ionising) level index
c                            (case 'h','r','s' & 'i')
c                            1st dim. transition index
c  input : (r*8)  aval()   = transition:
c                             a-value (sec-1)          (case ' ')
c                             neutral beam energy      (case 'h')
c                             not used                 (case 'p' & 'r')
c                            1st dim. transition index
c  input : (r*8)  scom(,)  = transition rate values
c                            gamma values             (case ' ' & 'p')
c                            rate coefft. (cm3 sec-1) (case 'h' & 'r')
c                            1st dim.: temp/energy/l index (cf. 'scef()')
c                                      according to iadftyp
c                            2nd dim.: - transition index
c  input : (i*4)  iadftyp  = adf04 type: 1=omega vs X, 3=upsilon vs te,
c                            4=non-maxwl. vs energy, 5=omega vs ef,
c                            6=omega_l vs l
c  input : (l*4)  lbeth    = .true.  => Bethe limit points available
c                            .false. => Bethe limit points not available
c  input : (r*8)  beth()   = Bethe infinite energy limit point
c                            1st dim.: transition index
c
c  input : (i*4)  ntext    = number of comment lines to be output
c  input : (c*(*)) ctext() = comment lines
c                            1st dim.: comment text line index
c
c  input : (l*4)  lx_untied= .true.  => exclude untied levels from output dataset
c                          = .false. => do not exclude untied levels
c
c
c  routines:
c          routine    source    brief description
c          -------------------------------------------------------------
c          xfesym     adas      converts nuclear charge to element symbol
c          xxwstr     adas      writes string to a unit with trailing
c                               blanks removed
c
c  author : Martin O'Mullane
c  date   : 28-01-2011
c
c
c  modification history:
c
c  version : 1.1
c  date    : 28-01-2011
c  modified: Martin O'Mullane
c            - first version.
c            - restricted to type 3 and 4 adf04 types.
c            - restricted to excitation effective collision strengths.
c
c  version : 1.2
c  date    : 22-08-2011
c  modified: Hugh Summers
c            - implemented orbital energy line output.
c
c  version : 1.3
c  date    : 04-11-2011
c  modified: Hugh Summers
c            - allowed for extended configuration string field length
c              obtained from cstrga input variable.
c            - allowed for extended xja field length and protected
c              isa, ila and xja values of zero.
c
c  version : 1.4
c  date    : 14-11-2011
c  modified: Hugh Summers
c            - restructured to a form compatible with variants in use
c              adas705 and elsewhere and so a substitute for them.
c            - added ndtext, ntext,ctext, iadftyp paramters.
c            - enabled correct signed output of 's', 't', 'r'
c              transition lines.
c
c  version : 1.5
c  date    : 18-12-2012
c  modified: Martin O'Mullane
c            - set a field width of 5 characters if the ionisation
c              potential is 0.0.
c
c  version : 1.6
c  date    : 27-03-2013
c  modified: Martin O'Mullane
c            - do not write blanks at end of comment lines.
c
c  version : 1.7
c  date    : 19-04-2016
c  modified: Hugh Summers
c            - allowed deletion of untied levels, controlled by logical
c              parameter lx_untied.
c            - set orb_string length to match ndqdn and il
c            - test internal idlev and idqdn match ndlev and ndqdn
c
c  version : 1.8
c  date    : 25-10-2016
c  modified: Martin O'Mullane
c            - Do not write A-value for R and S lines.
c            - The position of the '+' for S-lines depends on the
c              type of adf04 file. Assume the input transitions for
c              type 4 come in up/down pairs.
c
c  version : 1.9
c  date    : 27-07-2018
c  modified: Martin O'Mullane
c            - change ctext() to c*(*) rather than a fixed c*80.
c
c  version : 1.10
c  date    : 03-10-2018
c  modified: Martin O'Mullane
c            - Fix the calculation of the length of the configuration
c              field to be either 19 or the length of the longest
c              string in cstrga.
c
c  version : 1.11
c  date    : 03-10-2018
c  modified: Martin O'Mullane
c            - Define length of orb_string using idqdn in way it
c              is used to dimension qdorb. 
c
c  version : 1.12
c  date    : 18-02-2020
c  modified: Martin O'Mullane
c            - Allow 2 extra spaces between the end of cstrga and isa. 
c
c----------------------------------------------------------------------
       integer   idlev         , idqdn
       integer   len_cstrga_min
c----------------------------------------------------------------------
       real*8    dzero
c-----------------------------------------------------------------------
       parameter( idlev = 2000 , idqdn = 8   , len_cstrga_min = 19 )
       parameter( dzero = 1.0d-30 )
c-----------------------------------------------------------------------
       integer   ndlev       , ndtrn          , nvmax     , ndmet     ,
     &           ndqdn       , ndtext
       integer   i4unit      , lenstr         , iunit     , npl
       integer   iz          , iz0            , iz1       ,
     &           il          , nv             , itran     , iorb
       integer   it          , i              , j         , ipos
       integer   ndecb       , ndeca          , dtype     , prtmax
       integer   len_l       , len_s          , nfmt      , l_gap     ,
     &           len_j       , len_cstrga
       integer   ilen_orb_string              , iorb_s    , iorb_e
       integer   n           , l              , k
       integer   len_cstrga_max               , len_i1    , len_i2
       integer   iadftyp     , ntext
       integer   il_new      , im
c-----------------------------------------------------------------------
       logical   lbeth
       logical   lx_untied   , lsfirst
c----------------------------------------------------------------------
       real*8    wnomax      , z1             , eorb
c----------------------------------------------------------------------
       character esym*2      , c1*1           , xfesym*2  ,
     &           c9*9        , string*426     , blanks*426,
     &           strg2*44
       character fmt_02*7    , fmt_01*29      , fmt_03*20 ,
     &           fmt_04*48   , fmt_05*7       , fmt_06*10
       character orbfmt*9    , wavfmt*9       , 
     &           orb_string*(8*((idqdn*(idqdn+1))/2)+19)
c-----------------------------------------------------------------------
       integer   ietrn(ndtrn), istrn(ndtrn)
       integer   isa(ndlev)  , ila(ndlev)     ,
     &           i1a(ndtrn)  , i2a(ndtrn)
       integer   ipla(ndmet,ndlev)            , npla(ndlev)
       integer   il_mapa(idlev)               , il_unmapa(idlev)
c----------------------------------------------------------------------
       real*8    scef(nvmax)
       real*8    xja(ndlev)          , wa(ndlev)            ,
     &           aval(ndtrn)         , scom(nvmax,ndtrn)    ,
     &           zpla(ndmet,ndlev)   ,
     &           bwnoa(ndmet)        , beth(ndtrn)
       real*8    qdorb((ndqdn*(ndqdn+1))/2)       , qdn(ndqdn)
c-----------------------------------------------------------------------
       character tcode(ndtrn)*1 , cstrga(ndlev)*(*)
       character cpla(ndlev)*1  , cprta(ndmet)*9
       character ctext(ndtext)*(*)
c-----------------------------------------------------------------------
       logical   lqdorb((ndqdn*(ndqdn+1))/2)
c-----------------------------------------------------------------------
       data      blanks/' '/
c-----------------------------------------------------------------------
       external  xfesym
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c  initial internal dimension check
c-----------------------------------------------------------------------
         if (idlev.lt.ndlev) then
            write(i4unit(-1),1001) '(ndlev.gt.idlev) -'
            write(i4unit(-1),1002)
            stop
         endif
         if (idqdn.ne.ndqdn) then
            write(i4unit(-1),1001) '(ndqdn.ne.idqdn) -'
            write(i4unit(-1),1002)
            stop
         endif
c-----------------------------------------------------------------------
c  first line
c-----------------------------------------------------------------------

       fmt_01 = '(a2, 1h+,i2,2i10,?(f15.?,a?))'

       wnomax = bwnoa(1)
       prtmax = len(cprta(1))
       do i = 2, npl
         if (bwnoa(i).gt.wnomax) wnomax=bwnoa(i)
         if (len(cprta(i)).gt.prtmax) prtmax=len(cprta(i))
       enddo
       if (wnomax.GT.0.0D0) then
          ndecb = int(log10(float(int(wnomax))))
       else
          ndecb = 5
       endif
       ndeca = min(15-ndecb-2,4)

       write(fmt_01(18:18),'(i1)')npl
       write(fmt_01(24:24),'(i1)')ndeca
       write(fmt_01(27:27),'(i1)')prtmax

       esym = xfesym(iz0)
       call xxcase(esym(1:1), c1, 'uc')
       esym(1:1) = c1

       string = blanks
       write(string,fmt_01)esym,iz,iz0,iz1,
     &                     (bwnoa(i),cprta(i),i=1,npl)
       call xxwstr(iunit,string)

c-----------------------------------------------------------------------
c  identify untied levels for exclusion
c-----------------------------------------------------------------------
       if(lx_untied) then

           il_new = 0

           do i=1,il
             il_unmapa(il)=-1
           enddo

           do i=1,il
             do j=1,itran
                if(((tcode(j).eq.' ').or.(tcode(j).eq.'1').or.
     &              (tcode(j).eq.'3').or.(tcode(j).eq.'4').or.
     &              (tcode(j).eq.'p').or.(tcode(j).eq.'P')).and.
     &             ((i1a(j).eq.i).or.(i2a(j).eq.i))) then
                   il_new = il_new+1
                   il_mapa(i)=il_new
                   il_unmapa(il_new)=i
                   go to 5
                endif
                if(((tcode(j).eq.'h').or.(tcode(j).eq.'H').or.
     &              (tcode(j).eq.'r').or.(tcode(j).eq.'R').or.
     &              (tcode(j).eq.'i').or.(tcode(j).eq.'I').or.
     &              (tcode(j).eq.'s').or.(tcode(j).eq.'S')).and.
     &             (i2a(j).eq.i)) then
                   il_new = il_new+1
                   il_mapa(i)=il_new
                   il_unmapa(il_new)=i
                   go to 5
                endif
             enddo
             il_mapa(i) = -1
    5        continue

           enddo
       else
           do i=1,il
             il_mapa(i)=i
           enddo
       endif

c-----------------------------------------------------------------------
c  level information
c-----------------------------------------------------------------------

       len_cstrga = lenstr(cstrga(1))

       fmt_02 = '(f15.?)'
       fmt_03 = '(7x,"{",i1,"}",f4.2)'
       fmt_04 = "(1x,i?,1x,a??,3x,'(',i?,')',i?,'(',f?.1,')',a)"
       fmt_05 = '(1x,i?)'
       fmt_06 = '(a1,i?,i?)'


       nfmt = max(4, 3+int(dlog10(dfloat(il))))
       nfmt = min(9, nfmt)

       write(fmt_04(6:6),'(i1)')nfmt
       write(fmt_05(6:6),'(i1)')nfmt
       write(fmt_06(6:6), '(i1)')nfmt-1
       write(fmt_06(9:9), '(i1)')nfmt

       wnomax  = wa(il)
       len_l   = 1
       len_s   = 1
       len_j   = 4
       do i = 1, il
         if (wa(i).gt.wnomax) wnomax=wa(i)
         len_l = max(int(dlog10(dfloat(ila(i))+dzero))+1,len_l)
         len_s = max(int(dlog10(dfloat(isa(i))+dzero))+1,len_s)
         len_j = max(int(dlog10(xja(i)+dzero))+4,len_j)
         len_cstrga = max(lenstr(cstrga(i)),len_cstrga)
       enddo
       len_cstrga_max = max(len_cstrga, len_cstrga_min)

       write(fmt_04(12:13),'(i2.2)')len_cstrga_max
       write(fmt_04(6:6),'(i1)')nfmt
       write(fmt_04(23:23),'(i1)')len_s
       write(fmt_04(30:30),'(i1)')len_l
       write(fmt_04(37:37),'(i1)')len_j

       ndecb = int(log10(float(int(wnomax))))
       ndeca = min(15-ndecb-2,4)
       write(fmt_02(6:6),'(i1)')ndeca


       do i = 1,il
         im = il_mapa(i)
         if(im.gt.0) then

              string = blanks
              strg2  = blanks(1:44)
              write(strg2(1:16),fmt_02) wa(i)
              if(npla(i).le.0) then
                  strg2(20:22)='{X}'
              elseif((npla(i).eq.1).and.(cpla(i).eq.' ')) then
                  strg2(20:22)='{X}'
              elseif((npla(i).eq.1).and.((cpla(i).eq.'X').or.
     &                                   (cpla(i).eq.'x'))) then
                  strg2(20:22)='{X}'
              elseif((npla(i).ge.1).and.
     &            ((cpla(i).ne.' ').or.(cpla(i).ne.'X').or.
     &                                 (cpla(i).ne.'X'))) then
                  ipos=19
                  do j=1,npla(i)
                    write(strg2(ipos:ipos+8),'('' {'',i1,''}'',1f5.3)')
     &                     ipla(j,i),zpla(j,i)
                    ipos=ipos+9
                  enddo
              endif
              write(string,fmt_04)im,cstrga(i),isa(i),ila(i),xja(i),
     &                            strg2(1:lenstr(strg2))
              call xxwstr(iunit,string)

         endif

       end do

c-----------------------------------------------------------------------
c  level terminator and orbital line
c-----------------------------------------------------------------------

       z1  = iz1
       orb_string = blanks
       orbfmt = '(1x,f7.?)'
       write(orb_string, fmt_05)-1
       iorb_s = 2*nfmt+1
       do k = 1,iorb
         iorb_e = iorb_s + 7

         n=1
   10    if((n*(n+1))/2.lt.k) then
             n=n+1
             go to 10
         endif
         l=k-(n*(n-1))/2-1

         eorb = (z1/(float(n)-qdorb(k)))**2
         write(orbfmt(8:8),'(i1)') max(2,5-max(0,int(log10(
     &                            max(eorb,dzero)))))

         write(orb_string(iorb_s:iorb_e),orbfmt)eorb

         iorb_s = iorb_e+1
       enddo

       call xxwstr(iunit,orb_string)

c----------------------------------------------------------------------
c  zeff, adf04 type and temperature line
c----------------------------------------------------------------------

       if(iadftyp.eq.4) then
           dtype = 1
       else
           dtype = 0
       endif


       string = blanks
       write(string,'(f5.2,i5)')z1,iadftyp

       l_gap = 2 *(nfmt-4)

       ipos = 0
       do it = 1,nv
         write(c9,'(1pd9.2)') scef(it)
         string(17+l_gap+ipos:21+l_gap+ipos) = c9(1:5)
         string(22+l_gap+ipos:24+l_gap+ipos) = c9(7:9)
         ipos =ipos+8
       end do
       call xxwstr(iunit,string)

c----------------------------------------------------------------------
c  transition lines
c----------------------------------------------------------------------

       wavfmt  = '(1x,f7.?)'
       lsfirst = .TRUE.

       do i = 1,itran

         string = blanks
         call xxcase(tcode(i), c1, 'uc')
         write(string,fmt_06)c1,il_mapa(i2a(i)),il_mapa(i1a(i))
         len_i1 = int(dlog10(dfloat(il_mapa(i1a(i)))))
         len_i2 = int(dlog10(dfloat(il_mapa(i2a(i)))))

         if ((c1.eq.'R').or.(c1.eq.'L')) then
             write(string(2*nfmt-len_i1-1:2*nfmt-len_i1-1),'(a1)')'+'
         elseif ((c1.eq.'S').and.(il_mapa(i1a(i)).gt.0)) then
             if (iadftyp.eq.3) then
                write(string(2*nfmt-len_i1-1:2*nfmt-len_i1-1),'(a1)')'+'
             else
                if (lsfirst) then
                   write(string(2*nfmt-len_i1-1:
     &                          2*nfmt-len_i1-1),'(a1)')'+'
                   lsfirst = .FALSE.
                else
                   write(string(nfmt-len_i2-1:nfmt-len_i2-1),'(a1)')'+'
                   lsfirst = .TRUE.
                endif
             endif
         elseif ((c1.eq.'T').and.(il_mapa(i2a(i)).gt.0)) then
             write(string(2*nfmt-len_i1-1:2*nfmt-len_i1-1),'(a1)')'+'
         endif

         if(c1.eq.'L') then
             write(wavfmt(8:8),'(i1)') max(0,5-max(0,int(log10(
     &                            max(aval(i),dzero)))))
             write(string(9+l_gap:16+l_gap),wavfmt)aval(i)
         elseif (c1.eq.'R'.or.c1.eq.'S') then
         else
             write(c9,'(1pd9.2)') aval(i)
             string(9+l_gap:13+l_gap) = c9(1:5)
             string(14+l_gap:16+l_gap) = c9(7:9)
         endif

         ipos = 0
         do it = 1,nv
           write(c9,'(1pd9.2)') scom(it,i)
           string(17+l_gap+ipos:21+l_gap+ipos) = c9(1:5)
           string(22+l_gap+ipos:24+l_gap+ipos) = c9(7:9)
           ipos =ipos+8
         end do

         if (lbeth.and.((iadftyp.ne.4).or.(iadftyp.ne.6))) then
           write(c9,'(1pd9.2)') beth(i)
           string(17+l_gap+ipos:21+l_gap+ipos) = c9(1:5)
           string(22+l_gap+ipos:24+l_gap+ipos) = c9(7:9)
         endif

         call xxwstr(iunit,string)

       end do

c----------------------------------------------------------------------
c  terminators
c----------------------------------------------------------------------

       write(iunit,fmt_06)' ',-1
       write(iunit,fmt_06)' ',-1,-1

c----------------------------------------------------------------------
c write out comments
c----------------------------------------------------------------------

       do i=1,ntext
         call xxwstr(iunit, ctext(i))
       enddo

c----------------------------------------------------------------------
 1001 format(1x,31('*'),' xxwrto_04 error ',30('*')//
     &       1x,'Fault in input data file: ',a,i5,a)
 1002 format(/1x,29('*'),' program terminated ',29('*'))
c----------------------------------------------------------------------

       return

       end
