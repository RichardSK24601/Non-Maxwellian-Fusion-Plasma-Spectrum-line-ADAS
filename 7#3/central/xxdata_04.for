       subroutine xxdata_04( iunit  ,
     &                       ndlev  , ndtrn , ndmet   , ndqdn , nvmax ,
     &                       titled , iz    , iz0     , iz1   , bwno  ,
     &                       npl    , bwnoa , lbseta  , prtwta, cprta ,
     &                       il     , qdorb , lqdorb  , qdn   , iorb  ,
     &                       ia     , cstrga, isa     , ila   , xja   ,
     &                       wa     ,
     &                       cpla   , npla  , ipla    , zpla  ,
     &                       nv     , scef  ,
     &                       itran  , maxlev,
     &                       tcode  , i1a   , i2a     , aval  , scom  ,
     &                       beth   ,
     &                       iadftyp, lprn  , lcpl    , lorb  , lbeth ,
     &                       letyp  , lptyp , lrtyp   , lhtyp , lityp ,
     &                       lstyp  , lltyp , itieactn, ltied ,
     &                       irestyp
     &                  )

       implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 subroutine: xxdata_04 ******************
C
C  PURPOSE:  To fetch data from an adf04 data set and detect its main
C            characteristics. This is a fully inclusive version, based
C            on badata.for, detecting the following:
C
C                1. Multiple parent data on the first line including
C                   the j-resolved case
C                2. Supplementary parent assignment data on level
C                   lines for improved automatic ionisation calculation
C                3. Orbital energy data on the level terminator line
C                4. First bethe coefft. at end of e-transition lines for
C                   improved asymptotics
C                5. All transition line qualifiers , 'h',r',s',i','p'
C                   in upper or lower case;  ' ','1','2','3' electron
C                   impact transition types; multiple parents in 'r',
C                   'i',s' transition lines.
C                6. Doubly excited 'r' lines with Auger rate and resonance
C                   capture.
C                7. 'l' lines for dielectronic power correction to singly
C                   excited levels, including effective mean wavelength.
C
C  calling program: various
C
C  data:
C           The 'real' data in the file is represented in an abbreviated
C           form which omits the "d" or "e" exponent specifier.
C           e.g. 1.23d-06 or 1.23e-06 is represented as 1.23-06
C                6.75d+07 or 6.75e+07 is represented as 6.75+07
C
C           Therefore the form of each 'real' number in the data set is:
C                          n.nn+nn or n.nn-nn
C
C           The units used in the data file are taken as follows:
C
C           ionisation potential: wave number (cm-1)
C           index level energies: wave number (cm-1)
C           temperatures        : kelvin
C           a-values            : sec-1
C           gamma-values        :
C           rate coefft.        : cm3 sec-1
C
C
C  subroutine:
C
C  input : (i*4)  iunit   = unit to which input file is allocated
C  input : (i*4)  ndlev   = maximum number of levels that can be read
C  input : (i*4)  ndtrn   = max. number of transitions that can be read
C  input : (i*4)  nvmax   = max. number of temperatures that can be read in.
C
C  input : (i*4)  itieactn= 1 return data even if some levels are untied.
C                           0 default behaviour - terminate if untied
C                             levels are present.
C                           On output 1 if untied levels present
C                                     0 for no untied levels.
C
C  output: (c*3)  titled  = element symbol.
C  output: (i*4)  iz      =  recombined ion charge read
C  output: (i*4)  iz0     =         nuclear charge read
C  output: (i*4)  iz1     = recombining ion charge read
C                           (note: iz1 should equal iz+1)
C  output: (r*8)  bwno    = ionisation potential (cm-1) of lowest parent
C  output: (i*4)  npl     = number of parents on first line and used
C                           in level assignments
C  output: (r*8)  bwnoa() = ionisation potential (cm-1) of parents
C  output: (l*4)  lbseta()= .true.  - parent weight set for bwnoa()
C                           .false. - parent weight not set for bwnoa()
C  output: (r*8)  prtwta()= parent weight for bwnoa()
C  output: (c*9)  cprta() = parent name in brackets
C
C  output: (i*4)  il      = input data file: number of energy levels
C  output: (r*8)  qdorb() = quantum defects for orbitals. Negative quantum
C                           defects are forced to zero.
C                           1st dim: index for nl orbital (cf i4idfl.for)
C                           for missing l orbitals below the maximum l for
C                           an n-shell a value is calculated but lqdorb
C                           is set to .false.
C  output: (l*4)  lqdorb()= .true.  => source data available for qd.
C                         = .false. => source data not availabe qd.=0.0
C  output: (r*8)  qdn()   = quantum defect for n-shells.  non-zero only
C                           for adf04 files with orbital energy data
C                           1st. dim: n-shell (1<=n<=ndqdn)
C  output: (i*4)  iorb    = input data file: number of orbital energies
C
C  output: (i*4)  ia()    = energy level index number
C  output: (c*18) cstrga()= nomenclature/configuration for level 'ia()'
C  output: (i*4)  isa()   = multiplicity for level 'ia()'
C                           note: (isa-1)/2 = quantum number (s)
C  output: (i*4)  ila()   = quantum number (l) for level 'ia()'
C  output: (r*8)  xja()   = quantum number (j-value) for level 'ia()'
C                           note: (2*xja)+1 = statistical weight
C  output: (r*8)  wa()    = energy relative to level 1 (cm-1) for level
C                           'ia()'
C  output: (c*1)  cpla()  = char. specifying 1st parent for level 'ia()'
C                                integer - parent in bwnoa() list
C                                'blank' - parent bwnoa(1)
C                                  'x'   - do not assign a parent
C  output: (i*4)  npla()  = no. of parent/zeta contributions to ionis.
C                           of level
C  output: (i*4)  ipla(,) = parent index for contributions to ionis.
C                           of level
C                           1st dimension: parent index
C                           2nd dimension: level index
C  output: (i*4)  zpla(,) = eff. zeta param. for contributions to ionis.
C                           of level
C                           1st dimension: parent index
C                           2nd dimension: level index
C
C  output: (i*4)  nv      = input data file: number of gamma/temperature
C                           pairs for a given transition.
C  output: (r*8)  scef()  = input data file: electron temperatures (k)
C                           (initially just the mantissa. see 'itpow()')
C                           (note: te=tp=th is assumed)
C
C  output: (i*4)  itran   = input data file: number of transitions
C  output: (i*4)  maxlev  = highest index level in read transitions
C
C  output: (c*1)  tcode() = transition: data type pointer:
C                           ' ','1','2','3' => elec. impact trans.
C                           'p','P' => proton   impact   transition
C                           'h','H' => charge   exchange recombination
C                           'r','R' => free     electron recombination
C                           'i','I' => coll. ionis. from lower stage ion
C                           's','S' => Ionisation from current ion
C                           'l','L' => L-line for unresolved DR emissivity
C  output: (i*4)  i1a()   = transition:
C                            lower energy level index (case ' ' & 'p')
C                            signed parent index (case 'h','r','s' & 'i')
C  output: (i*4)  i2a()   = transition:
C                            upper energy level index (case ' ' & 'p')
C                            capturing level index (case 'h','r','s' & 'i')
C  output: (r*8)  aval()  = transition:
C                            a-value (sec-1)          (case ' ')
C                            neutral beam energy      (case 'h')
C                            not used             (case 'p','r' & 'i')
C  output: (r*8)  scom(,) = transition:
C                            gamma values             (case ' ' & 'p')
C                            rate coefft.(cm3 sec-1)(case 'h','r' & 'i')
C                            scaled rate coefft.(cm3 sec-1)(case 's')
C                           1st dimension - temperature 'scef()'
C                           2nd dimension - transition number
C  output: (i*4)  beth()  = transition
C                           1st Bethe coefficient     (case ' ','1','2')
C  output: (i*4)  iadftyp = adf04 type: 0=no transition data (NIST energies),
C                                       1=omega vs X, 3=upsilon vs Te,
C                                       4=non-maxwl. vs energy param.,
C                                       5=omega vs Ef, 6 omega_l vs l
C  output: (l*4)  lprn    = .true.  => multiple parent data on 1st line
C                         = .false. => multiple parent data not present
C  output: (l*4)  lcpl    = .true.  => parent assignment on level lines
C                         = .false. => parent assignment not present
C  output: (l*4)  lorb    = .true.  => orbital data on level terminator
C                         = .false. => orbital data not present
C  output: (l*4)  lbeth   = .true.  => bethe data on e-transition lines
C                         = .false. => bethe data not present
C  output: (l*4)  letyp   = .true.  => e- excitation transitions present
C  output: (l*4)  lptyp   = .true.  => p- excitation transitions present
C  output: (l*4)  lrtyp   = .true.  => recombination transitions present
C  output: (l*4)  lhtyp   = .true.  => cx transitions present
C  output: (l*4)  lityp   = .true.  => ionis. trans. from z-1 ion present
C  output: (l*4)  lstyp   = .true.  => ionis. or 3-body recom trans.
c                                      from current ion present
C  output: (l*4)  lltyp   = .true.  => 'l'-line for unresolved DR emissivity
C  output: (l*4)  ltied() = .true.  => specified level tied
C                         = .false. => specified level is untied
C                                      1st dimension - level index
C  output: (i*4)  irestyp = adf04 resolution:
C                                       0=cannot be determined
C                                       1=ls, 2=ic, 3=ca, 4=n
C
C          (i*4)  ndmet   = parameter = max. number of metastables allowed
C          (i*4)  ndqdn   = parameter = max. number of n-shells for quantum
C                                       defects
C          (i*4)  ntdim   = parameter = max. number of internal temperatures
C                                       (must equal nvmax)
C          (r*8)  dzero   = parameter = minimum value for 'aval()' and
C                                       'scom()' arrays = 1.0d-30
C
C          (i*4)  i4unit  = function (see routine selection below)
C          (i*4)  iqs     = x-sect data format selector
C          (i*4)  ifail   = failure number from xxpars and xxprs1
C          (i*4)  i       = general use.
C          (i*4)  iabt    = return code from 'r(fctn' (0 => no error)
C                           or from interrogation of 'c10'
C          (i*4)  j       = general use.
C          (i*4)  j1      = input data file - selected transition:
C                            lower energy level index (case ' ' & 'p')
C          (i*4)  j2      = input data file - selected transition:
C                            upper energy level index (case ' ' & 'p')
C                            capturing    level index (case 'h' & 'r')
C          (i*4)  lencst  = byte length of string cstrga()
C          (i*4)  iline   = energy level index for current line
C          (i*4)  irecl   = record length of input dataset (<=128)
C          (i*4)  itype   = resolution of parent metastables
C                             1 - ls resolved
C                             2 - lsj resolved
C                             3 - arbitrary resolution
C          (i*4)  iapow   = exponent of 'avalm'
C          (i*4)  igpow() = exponent of 'gamma()'
C          (i*4)  itpow() = temperatures - exponent
C                           note: mantissa initially kept in 'scef()'
C
C          (r*4)  zf      = should be equivalent to 'iz1'
C
C          (r*8)  avalm   = input data file - selected transition:
C                           mantissa of:   ('iapow' => exponent)
C                            a-value (sec-1)          (case ' ')
C                            neutral beam energy      (case 'h')
C                            not used              (case 'p','s' & 'i')
C          (r*8)  gamma() = input data file - selected transition:
C                           mantissa of: ('igpow()' => exponent)
C                            gamma values      (case ' ','1','2','3' & 'p')
C                            rate coefft.(cm3 sec-1)(case 'h','r','s' & 'i')
C                           dimension => temperature 'scef()'
C
C          (c*10) c10     = used to parse value for xja()
C          (c*7)  cdelim  = delimiters for input of data from headers
C          (c*90) c25     = used to parse value to cstrga()
C                            - nb: internal limit is now 90
C          (c*80) cline   = current energy level index parameter line
C          (c*75) string  = tail string of 1st data line for parsing
C          (c*44) strg1   = tail string of level spec lines for parsing
C          (c*500)buffer  = general string buffer storage
C          (c*3)  citpow()= used to parse values to itpow()
C          (c*5)  cscef() = used to parse values to scef()
C
C          (l*4)  ldata   = identifies  whether  the end of an  input
C                           section in the data set has been located.
C                           (.true. => end of section reached)
C          (l*4)  ltchr   = .true.  => current 'tcode()' = 'h' or 'r'
C                                                          's' or 'i'
C                         = .false. => current 'tcode()'.ne.'h' or 'r'
C                                                           's' or 'i'
C          (l*4)  ltcpr   = .true.  => current 'tcode()' = 'p' or
C                                                          's' or 'i'
C                         = .false. => current 'tcode()'.ne. 'p' or 'r'
C                                                            's' or 'i'
C          (l*4) lerror   = .true.  => untied level found
C                         = .false. => all levels tied
C          (l*4) lneg     = .true.  => negative gamma found
C                         = .false. => all are positive
C          (l*4) ltran    = .true.  => transitions are present
C                           .false. => no transition in file (probably NIST)
C
C
C note:        ltchr        ltcpr            tcode()
C              ----------------------------------------
C              .true.       .true.      =>   'r','i','s'
C              .true.       .false.     =>   'h', 'l'
C              .false.      .true.      =>   'p'
C              .false.      .false.     =>   ' ','1','2','3'
C
C        for a-values & gamma-values entries less than 'dzero' are taken
C        as being equal to dzero. this affects the 'aval()' and 'scom()'
C        arrays.
C
C routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          xxpars     ADAS      analyses the adf04 1st string for parents
C          xxprs1     ADAS      analyses the adf04 level string for ionis.
C          i4unit     ADAS      fetch unit number for output of messages
C          r8fctn     ADAS      converts from character to real variable
C          i4fctn     ADAS      converts from char. to integer  variable
C          xxslen     ADAS      finds string length excluding leading and
C                               trailing blanks
C          xxword     ADAS      parses a string into separate words
C                               for ' ()<>{}' delimiters
C
C AUTHOR:  Hugh Summers, University of Strathclyde
C          JA7.08
C          tel. 0141-548-4196
C
C DATE:    27/02/03
C
C UPDATE:
C
C VERSION:  1.2
C DATE:     10/09/2004
C
C MODIFIED: Allan Whiteford
C           - Extended code to handle J values greater than 10,000.
C             Actually, to allow greater spacing between the brackets
C             which delimit the J.
C
C VERSION:  1.3
C DATE:     26/11/2004
C
C MODIFIED: Paul Bryans and Allan Whiteford
C           - Fixed some of the comments.
C           - Do not re-order transition indices for type IV file
C           - Upped dimensions to allow 50 energies/temperatures
C
C VERSION:  1.4
C DATE:     26/11/2004
C
C MODIFIED: Allan Whiteford
C           - Changed dimension checks so that nvmax can be less
C             than ntdim.
C
C VERSION:  1.5
C DATE:     30/11/2004
C
C MODIFIED: Allan Whiteford
C           - Corrected flaw in logic introduced in version 1.4.
C
C VERSION:  1.6
C DATE:     28/07/2008
C
C MODIFIED: Allan Whiteford
C           - Allowed ionisation potential to be zero for bare nuclei
C           - Allowed configurations to be blank (also for bare nuclei)
C
C VERSION:  1.7
C DATE:     03/12/2007
C
C MODIFIED: Adam Foster
C           - Increased length of cline to allow for levels > 7d
C           - Fixed bug in detecting highest level present.
C
C VERSION:  1.8
C DATE:     31/07/2009
C
C MODIFIED: Martin O'Mullane
C           - Use i4 in error format statements to accommodate
C             the larger number of permitted levels.
C           - Return the configuration string (cstrga) as the length
C             of the incoming array. Note that there are internal limits
C             which restrict it to lie between 5 and 99 characters.
C           - Remove legacy warning for 'old argument list'.
C
C VERSION:  1.9
C DATE:     28/01/2010
C
C MODIFIED: Martin O'Mullane
C           - Accept type 5 as valid.
C           - Allow limit point to be negative (ie -4S/3).
C
C VERSION:  2.0
C DATE:     23/12/2010
C
C MODIFIED: Hugh Summers
C           - Modified ltcpr to let through aval for 'r' or 'R' lines.
C             Prevented end of transition list if j2<0 if tcode is not
C             blank.  Modified ltchr to include 'l' or 'L' lines. Change
C             allows through -ve j2 values for unattributed l-lines.
C
C VERSION:  2.1
C DATE:     23/09/2011
C
C MODIFIED: Hugh Summers
C           - Corrected transition index field length deduction using
C             3rd word rather than a + which is unsafe for type 5.
C             Extended to permit adf04 types 5 and 6.
C
C VERSION:  2.2
C DATE:     21/10/2011
C
C MODIFIED: Hugh Summers
C           - following Paul Bryans 22 Apr. 2005 incorporated provision
C             for type 4 files, make ionised index negative to
C             distinguish ionisation from 3-body recombination
C           - rationalised all xxdata_04 variants excluding Jj coupling
C             capability.
C
C VERSION:  2.3
C DATE:     14/11/2011
C
C MODIFIED: Hugh Summers
C           - corrected format specification f_1003 and inconsistency
C             of xja extraction as per specification.
C           - included handling of code 'T' or 't' from adf04 type 4
C             for three-body recombination replacing Bryans procedure
c             of version 2.2.
C           - included check of correct pairing of transition lines in
C             adf04 type 4. Warnings for non-pairing and over-pairing.
C           - modified ltied logic.
C
C VERSION:  2.4
C DATE:     23/12/2011
C
C MODIFIED: Hugh Summers
C           - modify bethe parameter handling, retain true zero, correct
C             false sign on -1.0d-30
C
C VERSION:  1.8
C DATE:     12/01/2012
C
C MODIFIED: Martin O'Mullane
C           - missing brackets in logic when identifying S-lines in
C             type 4 datasets.
C
C VERSION:  1.9
C DATE:     31/08/2012
C
C MODIFIED: Martin O'Mullane
C           - Return -1 for blank S (isa) and L (ila) fields.
C           - Add a check for NIST-style adf04 files which contain only
C             data for the energy levels.
C
C VERSION:  1.10
C DATE:     04/12/2012
C
C MODIFIED: Martin O'Mullane
C           - Incorrect dimensions for lmatcha and imatcha.
C
C VERSION:  1.11
C DATE:     18/12/2012
C
C MODIFIED: Martin O'Mullane
C           - Check if the number of orbital energies exceeds the
C             dimensions set by ndqdn.
C
C VERSION:  1.12
C DATE:     09/03/2015
C
C MODIFIED: Hugh Summers
C           - set iadftyp=0 for adf04 dataset without transition lines.
C
C VERSION:  1.13
C DATE:     02/03/2018
C
C MODIFIED: Matthew Bluteau
C           - Extend internal string to deal with larger top lines.
C           - Extend internal strg1 to deal with longer ionisation
C             pathway lines.
C
C VERSION:  1.14
C DATE:     17/06/2019
C
C MODIFIED: Martin O'Mullane
C           - Do not set lqdorb to .true. when qdorb is 0.0. Continue
C             to calculate a value for qdorb if l is less than the
C             highest l for that n-shell.
C
C VERSION:  1.15
C DATE:     18/06/2019
C
C MODIFIED: Martin O'Mullane
C           - Warn if there are negative gamma and set these to
C             effective zero. Negative values can appear in some
C             autostructure ic neutral files.
C
C VERSION:  1.16
C DATE:     17-02-2020
C
C MODIFIED: Martin O'Mullane
C           - For orbital angular quantum momentum numbers above 9, ie
C             l=n, o, q, r, t, u (limit of Condon and Shortley list) 
C             accept either a 2 digit value or a single letter
C             code (A, B, C etc. is output by autostructure) as valid.
C
C VERSION:  1.17
C DATE:     26-02-2020
C
C MODIFIED: Martin O'Mullane
C           - Add resolution type (ls, ic, ca or n) as an output
C             argument. Set to 0 if it cannot be determined.
C
C VERSION:  1.18
C DATE:     02-03-2020
C
C MODIFIED: Martin O'Mullane
C           - Improve logic when handling type=0 adf04 files. Add a
C             logic switch (ltran) rather than returning from inside
C             an if statement.
C
C VERSION:  1.19
C DATE:     16-09-2020
C
C MODIFIED: Martin O'Mullane
C           - Allow 'N=' when testing for n-resolved datasets. 
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      integer   ntdim
C-----------------------------------------------------------------------
      real*8    dzero
C-----------------------------------------------------------------------
      parameter( ntdim = 50 , dzero = 1.0d-30 )
C-----------------------------------------------------------------------
      integer   i4unit      , ifail          , itype      , irestyp,
     &          itieactn
      integer   iunit       , ndlev          , ndtrn      , nvmax  ,
     &          ndmet       , ndqdn          ,
     &          iz          , iz0            , iz1        ,
     &          il          , nv             , itran      ,
     &          maxlev      , npl            , iadftyp
      integer   i4idfl      , n              , l          , nlast
      integer   iline       , irecl          , iorb
      integer   ifirst(3)   , ilast(3)       , iword      , nwords ,
     &          ilen_index  , ilen_config    , ilen_s     ,
     &          ilen_l      , ilen_j         , ilen_tran
      integer   iqs         , i              , iabt       , j      ,
     &          j1          , j2             , k          ,
     &          lencst      , lcline         , l1         , l2     ,
     &          ibeth       , iplus          ,
     &          iapow       , igpow(ntdim)   , itpow(ntdim)
      integer   ia(ndlev)   , isa(ndlev)     , ila(ndlev) ,
     &          i1a(ndtrn)  , i2a(ndtrn)     ,
     &          ipla(ndmet,ndlev)            , npla(ndlev)
      integer   locpar(4)   , imatcha(ndtrn)
      integer   iwarn_lttyp , iwarn_umatch   , iwarn_omatch
      integer   n_t1        , n_t2           , n_t3       , n_t4   ,
     &          n_t5
C-----------------------------------------------------------------------
      real*8    zf
      real*8    r8fctn
      real*8    bwno        , bwnoa(ndmet)   , prtwta(ndmet)  , bethm  ,
     &          avalm
      real*8    scef(nvmax) , gamma(ntdim)
      real*8    xja(ndlev)  , wa(ndlev)      ,
     &          beth(ndtrn) , aval(ndtrn) , scom(nvmax,ndtrn)  ,
     &          zpla(ndmet,ndlev)
      real*8    qdorb((ndqdn*(ndqdn+1))/2)       , qdn(ndqdn)
C-----------------------------------------------------------------------
      character titled*3    , tcode(ndtrn)*1 , cstrga(ndlev)*(*)
      character c10*20      , cdelim*7       , c25*99       , c1*1     ,
     &          c2*2        , c2a*2
      character cline*296   , string*500     , strg1*150    , buffer*500
      character f_1003*35   , f_1005*28      , f_1015*4
      character f_1008*22   , f_1014*28      , f_1020*10
      character citpow(ntdim)*3              , cscef(ntdim)*5
      character cpla(ndlev)*1,cprta(ndmet)*9
      character caj1*9       , caj2*9        ,
     &          codei*1      , codej*1       , codet*1
C-----------------------------------------------------------------------
      logical   ldata       , ltchr          , ltcpr      , lerror
      logical   lprn        , lcpl           , lorb       , lbeth
      logical   ltele       , lneg           , ltran
      logical   letyp       , lptyp          , lrtyp      , lhtyp     ,
     &          lityp       , lstyp          , lttyp      , lltyp
C-----------------------------------------------------------------------
      logical   ltied(ndlev), lbseta(ndmet)
      logical   lqdorb((ndqdn*(ndqdn+1))/2)
      logical   lmatcha(ndtrn)
C-----------------------------------------------------------------------
      data      cdelim/ ' ()<>{}'/
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C  initial internal dimension check
C-----------------------------------------------------------------------
         if (nvmax.gt.ntdim) then
            write(i4unit(-1),1001) '(nvmax.gt.ntdim) -'
            write(i4unit(-1),1002)
            stop
         endif
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C For configuration strings use the length as passed into routine but
C ignore anything greater than character 99 (internal parsing limit -
C see c25 varaible). Impose a minimum of 5 characters.
C-----------------------------------------------------------------------

      lencst = len(cstrga(1))

      if (lencst.lt.5) then
         write(i4unit(-1),1014)
         write(i4unit(-1),1002)
         write(0,*)'continue regardless'
c         stop
      endif

      lencst = max0( len(cstrga(1)) , 99 )

C----------------------------------------------------------------------
C identify the record length of the input dataset (old-style = 80 ? )
C (maximum currently active usage is 500 bytes)
C----------------------------------------------------------------------
      inquire(iunit,recl=irecl)
      irecl = min0( irecl , 500 )

C-----------------------------------------------------------------------
C  set most logicals to .false. initially but assume that there will
C  be transitions in the file (ltran=.true.).
C  set warnings count
C-----------------------------------------------------------------------
      lprn = .false.
      lcpl = .false.
      lorb = .false.
      lbeth = .false.

      letyp = .false.
      lptyp = .false.
      lrtyp = .false.
      lhtyp = .false.
      lityp = .false.
      lstyp = .false.
      lttyp = .false.
      lltyp = .false.

      iwarn_lttyp = 0

      lneg = .false.
      ltran = .true.

C-----------------------------------------------------------------------
C input ion specifications and parse first data card string
C-----------------------------------------------------------------------
      read(iunit,1000) titled, iz, iz0, iz1, string
      ifail = 0
      call xxpars( ndmet , string , npl  , bwnoa , lbseta  , prtwta  ,
     &             cprta , ifail  , itype )
      if(npl.ge.1) lprn=.true.

      if(ifail.eq.1) then
         write(i4unit(-1),1001) 'ifail=1 from xxpars - stop '
         write(i4unit(-1),1002)
         stop
      endif
      bwno=bwnoa(1)
      if(bwno.le.0.0d0 .and. iz .ne. iz0) then
         write(i4unit(-1),1001) 'ionisation potential .le. 0'
         write(i4unit(-1),1002)
         stop
      endif

C-----------------------------------------------------------------------
C  read in energy level specifications
C  energy level input information is terminated by iline=-1.
C-----------------------------------------------------------------------
      ldata=.true.
      read (iunit,'(a)') cline
      backspace(iunit)

C-----------------------------------------------------------------------
C  pre-analyse formatting of level line
C
C  Search first configuration line from the end to locate
C  the 4 parenthesis bracketing the quantum numbers
C-----------------------------------------------------------------------

      lcline=len(cline)

      k=4
      do i = lcline, 1, -1
         if (k.gt.0.and.(cline(i:i).eq.'('.or.cline(i:i).eq.')')) then
            locpar(k) = i
            k = k - 1
         endif
      end do

      iword = 1
      j     = 1
      call xxword( cline, ' '  , iword   , j       ,
     &             ifirst(1)  , ilast(1) , nwords
     &           )

      ilen_index  = len(cline(1:ilast(1)))
      ilen_config = locpar(1) - ilen_index - 2
      ilen_s      = locpar(2) - locpar(1) - 1
      ilen_l      = locpar(3) - locpar(2) - 1
      ilen_j      = locpar(4) - locpar(3) + 1        ! include the ()

c      write(i4unit(-1),*)
c     &                 'ilen_index,ilen_config,ilen_s,ilen_l,ilen_j=',
c     &                         ilen_index,ilen_config,ilen_s,ilen_l,ilen_j

      write(f_1003,'(a2,i1,a6,i2.2,a5,i2.2,a5,i2.2,a2,i2.2,a3)')
     &             '(i',ilen_index,',1x,1a',ilen_config,',1x,a',ilen_s,
     &             ',1x,a', ilen_l,',a',ilen_j,',a)'

      write(f_1015,'(a2,i1,a1)')'(i',ilen_index,')'

C-------------------------------
C Step through the energy levels
C-------------------------------

      maxlev=0

         do 1 i=1,ndlev

            ltied(i) = .false.

            if (ldata) then
               cline=' '
               read (iunit,'(a)') cline
               read (cline,f_1015)  iline

                  if (iline.le.0) then
                     ldata=.false.
                     il=i-1
C-------------------------------------------------
C  check for orbitals and evaluate quantum defects
C-------------------------------------------------
                     j     = 1
                     iword = 1
                     lcline=len(cline)

                     call xxword( cline(ilen_index+1:lcline),
     &                            cdelim     , iword    ,j      ,
     &                            ifirst(1)  , ilast(1) , iorb   )

                     if (iorb.gt.(ndqdn*(ndqdn+1))/2) then
                         write(i4unit(-1),1001)
     &                          'adf04 data set contains > ',
     &                          (ndqdn*(ndqdn+1))/2,' orbital energies'
                         write(i4unit(-1),1002)
                         stop
                     endif

                     if(iorb.gt.0) then
                         read(cline(ilen_index+1:lcline),*)
     &                        (qdorb(k),k=1,iorb)
                         do 20 l1 = 1,ndqdn
                           l = l1-1
                           nlast = l
                           do 10 n = l1,ndqdn
                             k = i4idfl(n,l)
                             if(dabs(qdorb(k)).gt.0.0d0) then
                                 qdorb(k) = dmax1((-dfloat(iz1)/
     &                             dsqrt(qdorb(k))+dfloat(n)),0.0d0)
                                 lqdorb(k) = .true.
                                 nlast = n
                             else
                                 lqdorb(k) = .false.
                                 if(nlast.gt.l) then
                                     qdorb(k) = qdorb(i4idfl(nlast,l))
                                 else
                                     qdorb(k) = 0.0d0
                                 endif
                             endif
   10                      continue
   20                    continue

                         do 40 n = 1,ndqdn
                           qdn(n) = 0.0d0
                           do 30 l1 = 1,n
                             l=l1-1
                             k = i4idfl(n,l)
                             qdn(n) = qdn(n) + (2.0d0*l+1.0d0)
     &                                /(n-qdorb(k))**2
   30                      continue
                           qdn(n) = -n/dsqrt(qdn(n)) + n
   40                    continue
                         lorb = .true.
                     endif
C-------------------------------------------------
C  energy level index 'iline' should equal 'i'
C-------------------------------------------------
                  elseif (iline.ne.i) then
                     write(i4unit(-1),1001) 'energy level index',iline,
     &                                      ' out of order'
                     write(i4unit(-1),1002)
                     stop
                  else
                     if (ilen_l.EQ.1) then
                        read (cline,f_1003) ia(i)  , c25       ,
     &                                      c2     , c1        ,
     &                                      c10    , strg1
                        if (c1.NE.' ') then
                           k = index('ABCDEF', c1)
                           if (k.eq.0) then
                              read(c1,*)ila(i)
                           else
                              ila(i) = 10 + (k-1)
                           endif
                        else
                           ila(i) = -1
                        endif
                     elseif (ilen_l.EQ.2) then
                        read (cline,f_1003) ia(i)  , c25       ,
     &                                      c2     , c2a        ,
     &                                      c10    , strg1
                        if (c2a.NE.'  ') then
                           k = index('ABCDEF', c2a)
                           if (k.eq.0) then
                              read(c2a,*)ila(i)
                           else
                              ila(i) = 10 + (k-1)
                           endif
                        else
                           ila(i) = -1
                        endif
                     endif

                     call xxslen(c2(1:ilen_s), l1, l2)
                     if (l1+l2.GT.0) then
                        read(c2(1:ilen_s),*)isa(i)
                     else
                        isa(i) = -1
                     endif


                     ifail = 0
                     call xxprs1(ndmet,strg1,wa(i),cpla(i),npla(i),
     &                           ipla(1,i),zpla(1,i),ifail)
                     if(ifail.eq.1) then
                         write(i4unit(-1),1001)
     &                         'ifail=1 from xxprs1 - stop '
                         write(i4unit(-1),1002)
                         stop
                     endif
                     if(ifail.eq.0) lcpl = .true.

C-----------------------------------------------------
C  Reduce cstrga to lencst
C-----------------------------------------------------

                     if (c25(1:2).eq.'  ') c25=c25(3:25)

                     call xxslen(c25,ifirst(1),ilast(1))
                     ifirst(1)=max0(ilast(1)-(lencst-1),ifirst(1))
                     if (ilast(1).gt.0) then
                        cstrga(i)=c25(ifirst(1):ilast(1))
                     else
                        cstrga(i)=''
                     endif
C-------------------------------------------------

                     j = index(c10,')') - 1

                        if (c10(1:j).eq.' ') then
                           iabt   = -1
                        else if (j.gt.1) then
                           xja(i) = r8fctn( c10(2:j) , iabt )
                        else
                           iabt   = j - 1
                        endif

                        if (iabt.ne.0) then
                           write(i4unit(-1),1001)
     &                                   'adf04 data set level ',i,
     &                                   ' has invalid weight'
                              if (iabt.eq.-1) then
                                 write(i4unit(-1),1012)
     &                                    'no number found',c10
                              else if (iabt.eq.-2) then
                                 write(i4unit(-1),1012)
     &                                    'no right bracket found',c10
                              else
                                 write(i4unit(-1),1012)
     &                                    'invalid number found' , c10
                              endif
                           write(i4unit(-1),1002)
                           stop
                        endif
                    maxlev=max0(iline,maxlev)
C-------------------------------------------------
                  endif
            endif
    1    continue

         if (ldata) then
            read (iunit,1004) i
               if (i.gt.0) then
                   write(i4unit(-1),1001) 'adf04 data set contains > ',
     &                                    ndlev,' energy levels'
                   write(i4unit(-1),1002)
                  stop
               else
                  il = ndlev
               endif
         endif

C-----------------------------------------------------------------------
C  read in temperatures (kelvin) and data format selector
C
C  check if this is a NIST dataset, which only contains energy levels
C  if so return setting transition data to zero.
C-----------------------------------------------------------------------

      buffer = ' '
      read(iunit,'(a)') buffer

C End of data check

      iword = 1
      j     = 3
      c1    = ' '
      call xxword( buffer  , c1     , iword     ,
     &             j       ,
     &             ifirst  , ilast  , nwords
     &            )
      if (nwords.eq.2) then

         if (buffer(ifirst(1):ilast(1)).eq.'-1'.and.
     &       buffer(ifirst(2):ilast(2)).eq.'-1') then

            nv    = 0
            itran = 0
            do i = i, ndtrn
               tcode(i) = ' '
               i1a(i)   = 0
               i2a(i)   = 0
               aval(i)  = 0.0
               beth(i)  = 0.0
               do j = 1, nvmax
                  scom(j,i) = 0.0
               end do
            end do

	    iadftyp = 0
            ltran   = .false.

         endif

      endif

C If there are no transitions, eg NIST datasets, skip to end (for the
C coupling resolution determination).

      if (.not.ltran) goto 100
      
C Transition data present so continue

      iword = 3
      j     = 3
      c1    = ' '
      call xxword( buffer  , c1     , iword     ,
     &             j       ,
     &             ifirst  , ilast  , nwords
     &            )

      ilen_tran=ifirst(1)-10

      f_1005='(f5.2,4x,i1,?x,?x,50(a5,a3))'
      f_1008='(a1,i?,i?,51(f5.2,i3))'
      f_1014='(a1,i?,i?,1f8.0,50(f5.2,i3))'
      f_1020='(a1,a?,a?)'
      write(c1,'(i1)')(ilen_tran-2)/2
      f_1005(13:13)=c1
      f_1005(16:16)=c1
      write(c1,'(i1)')ilen_tran/2-1
      f_1008(6:6)=c1
      f_1014(6:6)=c1
      f_1020(6:6)=c1
      write(c1,'(i1)')ilen_tran/2
      f_1008(9:9)=c1
      f_1014(9:9)=c1
      f_1020(9:9)=c1
      read(buffer,f_1005) zf, iqs, (cscef(i),citpow(i),i=1,ntdim)
      iadftyp = iqs

C-------------------------------------------------
C check data format selector is valid
C-------------------------------------------------
         if((iqs.ne.1).and.(iqs.ne.3).and.(iqs.ne.4).and.
     &      (iqs.ne.5).and.(iqs.ne.6)) then
            write(i4unit(-1),1001) '(iqs.ne.1, 3, 4, 5 or 6) -'
            write(i4unit(-1),1006)
     &      'file contains invalid data format selector ',
     &      '(must equal 1, 3, 4, 5 or 6)'
            write(i4unit(-1),1002)
            stop
         endif
C-------------------------------------------------
C check 'zf' equals iz1
C-------------------------------------------------
         if(nint(zf).ne.iz1) then
            write(i4unit(-1),1001) '(zf.ne.iz1) -'
            write(i4unit(-1),1006)
     &      'file contains inconsistent recombining ion charge values'
            write(i4unit(-1),1002)
            stop
         endif
C-------------------------------------------------
C  check that nv isn't greater than nvmax
C-------------------------------------------------
         nv=0
         do 9 j=1,ntdim
            if (cscef(j).ne.' ') then
               nv = j
            endif
    9    continue
         if (nv.gt.nvmax) then
            write(i4unit(-1),1001) '(nv.gt.nvmax) - nv =',nv
            write(i4unit(-1),1002)
            stop
         endif
C-------------------------------------------------
C identify the number of temperatures values input
C-------------------------------------------------
         nv = 0

         do 2 j=1,nvmax
            if (cscef(j).eq.' ') then
               scef(j)  = 0.0d0
               itpow(j) = 0
            else
               read(cscef(j) ,'(f5.2)') scef(j)
               read(citpow(j), '(i3)')  itpow(j)
               if (scef(j).gt.0.0d0) nv = j
            endif
    2    continue
C-------------------------------------------------
C  check that at least one valid temperature exists.
C-------------------------------------------------
         if (nv.le.0) then
            write(i4unit(-1),1001) '(nv.le.0) -'
            write(i4unit(-1),1006)
     &      'no valid temperature values found (none > 0.0)'
            write(i4unit(-1),1002)
            stop
         endif
C-------------------------------------------------
C  combine input mantissa and exponent for temperatures
C-------------------------------------------------
         do 4 i=1,nv
            scef(i)  = scef(i)  * ( 10.0D0**dble(itpow(i)) )
    4    continue
C-------------------------------------------------
C  check temperatures are in a strictly increasing
C  monotonic sequence (required for the later use
C  of the minimax and spline routines)
C-------------------------------------------------
         do 5 i=2,nv
            if (scef(i).le.scef(i-1)) then
               write(i4unit(-1),1007) '(scef(',i,').le.scef(',i-1,')) -'
               write(i4unit(-1),1006)
     &      'temperatures are not in a strict monotonic ascending order'
               write(i4unit(-1),1002)
               stop
             endif
    5    continue

C-----------------------------------------------------------------------
C read in transition specifications
C-----------------------------------------------------------------------
C
      lbeth = .false.
C      maxlev = 0
      ldata = .true.

         do 6 i=1,ndtrn
            if (ldata) then
               buffer = ' '
               read(iunit,'(a)')buffer
               read(buffer(1:1),'(1a1)')tcode(i)
               call xxslen(buffer,ifirst(1),ilast(1))
               if(ilast(1).gt.(8*nv+16+ilen_tran-8)) then
                  read(buffer,f_1008) tcode(i) , j2 , j1 , avalm ,
     &                iapow,(gamma(j),igpow(j),j=1,nv),bethm,ibeth
                  lbeth = .true.
               elseif((tcode(i).ne.'l').and.(tcode(i).ne.'L')) then
                  read(buffer,f_1008) tcode(i) , j2 , j1 , avalm ,
     &                iapow,(gamma(j),igpow(j),j=1,nv)
                  bethm = 0.0d0
                  ibeth = 0
               else
                  read(buffer,f_1014) tcode(i) , j2 , j1 , avalm ,
     &                          (gamma(j),igpow(j),j=1,nv)
                  iapow = 0
               endif
C-------------------------------------------------
C transition input information is terminated by j2=-1
C-------------------------------------------------
                  if ((j2.le.0).and.(tcode(i).eq.' ')) then
                     ldata=.false.
                     itran=i-1
                  else

C-----------------------------------------------------------------------
C set transition type detection logical and identify 'tcode()' type.
C-----------------------------------------------------------------------

                     if((tcode(i).eq.' ').or.(tcode(i).eq.'1')
     &                                   .or.(tcode(i).eq.'2')
     &                                   .or.(tcode(i).eq.'3'))
     &                                   letyp=.true.
                     if((tcode(i).eq.'p').or.(tcode(i).eq.'P'))
     &                                   lptyp=.true.
                     if((tcode(i).eq.'r').or.(tcode(i).eq.'R'))
     &                                   lrtyp=.true.
                     if((tcode(i).eq.'h').or.(tcode(i).eq.'H'))
     &                                   lhtyp=.true.
                     if((tcode(i).eq.'i').or.(tcode(i).eq.'I'))
     &                                   lityp=.true.
                     if((tcode(i).eq.'s').or.(tcode(i).eq.'S'))
     &                                   lstyp=.true.
                     if((tcode(i).eq.'T').or.(tcode(i).eq.'t'))
     &                                   lttyp=.true.
                     if((tcode(i).eq.'l').or.(tcode(i).eq.'L'))
     &                                   lltyp=.true.

                     ltchr = (tcode(i).eq.'h') .or. (tcode(i).eq.'r')
     &                                         .or. (tcode(i).eq.'i')
     &                                         .or. (tcode(i).eq.'s')
     &                                         .or. (tcode(i).eq.'t')
     &                                         .or. (tcode(i).eq.'l')
     &                                         .or. (tcode(i).eq.'H')
     &                                         .or. (tcode(i).eq.'R')
     &                                         .or. (tcode(i).eq.'I')
     &                                         .or. (tcode(i).eq.'S')
     &                                         .or. (tcode(i).eq.'T')
     &                                         .or. (tcode(i).eq.'L')
                     ltcpr = (tcode(i).eq.'p') .or. (tcode(i).eq.'i')
     &                                         .or. (tcode(i).eq.'s')
     &                                         .or. (tcode(i).eq.'t')
     &                                         .or. (tcode(i).eq.'P')
     &                                         .or. (tcode(i).eq.'I')
     &                                         .or. (tcode(i).eq.'S')
     &                                         .or. (tcode(i).eq.'T')
                     ltele = (tcode(i).eq.' ') .or. (tcode(i).eq.'1')
     &                                         .or. (tcode(i).eq.'2')
     &                                         .or. (tcode(i).eq.'3')

C-------------------------------------------------
C for type 4 files, permit three-body recombination lines
C     ensure tcode and indexing are correct for older
C     data sets which use only the 'S' or 's' notation
C-------------------------------------------------
                     if (((tcode(i).eq.'T').or.(tcode(i).eq.'t')).and.
     &                   (iwarn_lttyp.le.0).and.(iqs.ne.4)) then
                         iwarn_lttyp = 1
                         write(i4unit(-1),1003)
     &        'T transition lines only permitted for adf04 type 4'
                         write(i4unit(-1),1005)
                     endif


                     if (((tcode(i).eq.'S').or.(tcode(i).eq.'s')).and.
     &                   (iqs.eq.4)) then
                        read(buffer,f_1020) tcode(i) , caj2, caj1
                        iplus = index(caj1,'+')
                        if (iplus.ne.0) then
                            tcode(i)='t'
                            lttyp   = .true.
                        endif
                        iplus = index(caj2,'-')
                        if (iplus.ne.0) then
                            tcode(i)='t'
                            lttyp   = .true.
                        endif
                     endif

C-------------------------------------------------
C identify that specified excited levels are tied.
C-------------------------------------------------

                     if(.not.ltchr) then
                         ltied(j1) = .true.
                         ltied(j2) = .true.
                     endif

C-------------------------------------------------
C combine input mantissa and exponent for a-value, gammas & bethe coeff.
C-------------------------------------------------

                     if (.not.ltcpr) then
                        aval(i)=avalm*( 10.0D0**dble(iapow) )
                        aval(i)=dmax1(aval(i),dzero)

                        if(lbeth) then
                            beth(i) = bethm*( 10.0D0**dble(ibeth) )
                            if (bethm.gt.0.0d0) then
                                beth(i) = dmax1(beth(i),dzero)
		            elseif(bethm.eq.0.0d0) then
			        beth(i) = 0.0d0
                            else
                                beth(i) = -1.0*dmax1(abs(beth(i)),dzero)
                            endif
                        endif

                     endif

                     do j=1,nv
                        if (gamma(j).lt.0.0) lneg = .true.
                        scom(j,i)=gamma(j)*( 10.0D0**dble(igpow(j)) )
                        scom(j,i)=dmax1(scom(j,i),dzero)
                     end do

C-------------------------------------------------
C  find highest index level for transitions.
C  (tcode() = ' ' or 'p' only)
C-------------------------------------------------

C                     if (.not.ltchr) maxlev=max0(j1,j2,maxlev)

C-------------------------------------------------
C  check if energy levels are in correct order.
C  if not reverse except for type 4 files.
C-------------------------------------------------

                     if (ltchr) then
                        i2a(i)=j2
                        i1a(i)=j1
                     elseif (wa(j2).lt.wa(j1).and.(iqs.ne.4)) then
                        i1a(i)=j2
                        i2a(i)=j1
                        write(i4unit(-1),1009)
     &                    'xxdata_04 message: Upper and ',
     &                    'lower energy levels have been reversed.',
     &                    'Correct levels are - upper = ',j1,
     &                    '  lower = ',j2
                     else
                        i1a(i)=j1
                        i2a(i)=j2
                     endif
C-----------------------------------------------------------------------
                  endif
            endif
    6    continue

C-----------------------------------------------------------------------
C warn if negative gamma are found
C-----------------------------------------------------------------------

         if (lneg) then
            write(i4unit(-1),1003)'Some negative gammas found. '//
     &                            'These are set to effective zero.'
            write(i4unit(-1),1005)
         endif

C-----------------------------------------------------------------------

         if (ldata) then
            read (iunit,1010) i
               if (i.gt.0) then
                  write(i4unit(-1),1001) 'adf04 data set contains > ',
     &                                   ndtrn,' transitions'
                  write(i4unit(-1),1002)
                  stop
               else
                  itran = ndtrn
               endif
         endif

C-----------------------------------------------------------------------
C check that there are no untied levels
C-----------------------------------------------------------------------

      lerror = .false.

      do i=1,maxlev
         if (.not.ltied(i)) then
            lerror = .true.
            write(i4unit(-1),1011) i
         endif
      end do

C return or terminate based on itieact setting

      if (itieactn.ne.0) then
         if (lerror) then
            write(i4unit(-1),1003)'Untied levels present'//
     &                            ' fix before proceeding'
            itieactn = 1
         else
            itieactn = 0
         endif
      else
         if (lerror) then
            write(i4unit(-1),1002)
            stop
         endif
      endif


C-----------------------------------------------------------------------
C check that type 4 dataset comprises matched up/down pairs
C-----------------------------------------------------------------------

       if(iadftyp.eq.4) then

            iwarn_umatch = 0
            iwarn_omatch = 0
            do i=1,itran
              lmatcha(i) = .false.
              imatcha(i) = 0
            enddo
            do i=1,itran-1
              if(.not.lmatcha(i)) then
                  call xxcase(tcode(i),codei,'lc')
                  if((codei.ne.'r').and.(codei.ne.'l').and.
     &               (codei.ne.'i').and.(codei.ne.'h'))then
                        if(codei.eq.'s') then
                            codet='t'
                        elseif (codei.eq.'t') then
                            codet = 's'
                        else
                            codet = codei
                        endif
                        do j = i+1,itran
                           if(.not.lmatcha(j)) then
                               call xxcase(tcode(j),codej,'lc')
                               if((i1a(j).eq.i2a(i)).and.
     &                            (i1a(i).eq.i2a(j)).and.
     &                            (codet.eq.codej))then
                                   lmatcha(i) = .true.
                                   imatcha(i) = imatcha(i)+1
                                   lmatcha(j) = .true.
                                   imatcha(j) = imatcha(i)
                               endif
                           endif
                        enddo
                  else
                      lmatcha(i) = .true.
                      imatcha(i) = 0
                  endif
              endif
            enddo

            do i=1,itran
              if(.not.lmatcha(i)) iwarn_umatch = iwarn_umatch+1
              if(imatcha(i).gt.1) iwarn_omatch = iwarn_omatch+1
            enddo
            if (iwarn_umatch.gt.0)then
                write(i4unit(-1),1001)'There are', iwarn_umatch,
     &          ' unmatched transitions in adf04 type 4 dataset'
                write(i4unit(-1),1005)
            endif
            if (iwarn_omatch.gt.0)then
                write(i4unit(-1),1001)'There are', iwarn_omatch,
     &          ' overmatched transitions in adf04 type 4 dataset'
                write(i4unit(-1),1005)
            endif

       endif

C-----------------------------------------------------------------------
C Try to determine the resolution of the adf04 file - LS, IC, CA or n.
C-----------------------------------------------------------------------

 100   continue

       irestyp = 0

       n_t1 = 0
       n_t2 = 0
       n_t3 = 0
       n_t4 = 0
       n_t5 = 0
       do i = 1, il
          n_t1 = n_t1 + index(cstrga(i), 'n=') + index(cstrga(i), 'N=')
          n_t2 = n_t2 + isa(i)
          n_t3 = n_t3 + ila(i)
          if (xja(i).gt.(ila(i) + (isa(i)-1.0)/2.0)) n_t4 = n_t4 + 1
          if (xja(i).lt.abs(ila(i) - (isa(i)-1.0)/2.0)) n_t5 = n_t5 + 1
       end do

C n-resolved test is based on convention of cstrga having n=?

       if (n_t1.eq.il) irestyp = 4

C configuration average is indicated by all of isa and ila being zero

       if (irestyp.eq.0) then
          if (n_t2.eq.0.and.n_t3.eq.0) irestyp = 3
       endif

C For ic |L-s| <= J <= L+S otherwise it's LS

       if (irestyp.eq.0) then
          if ((n_t4.GT.0).OR.(n_t5.GT.0)) then
             irestyp = 1
          else
             irestyp = 2
          endif
       endif

C-----------------------------------------------------------------------

 1000 format(1a3,i2,2i10,1a300)
 1001 format(1x,31('*'),' xxdata_04 error ',30('*')//
     &       1x,'Fault in input data file: ',a,i7,a)
 1002 format(/1x,29('*'),' program terminated ',29('*'))
 1003 format(1x,30('*'),' xxdata_04 warning ',29('*')//
     &       1x,'Fault in input data file: ',a)
 1004 format(i5)
 1005 format(/1x,29('*'),' program continues ',29('*'))
 1006 format(1x,a)
 1007 format(1x,31('*'),' xxdata_04 error ',30('*')//
     &       1x,'Fault in input data file: ',a,2(i1,a))
 1009 format(a,a/17x,2(a,i4)/)
 1010 format(1x,i3)
 1011 format(1x,31('*'),' xxdata_04 error ',30('*')//
     &       1x,'Error in input data file: level ',i4,' is untied')
 1012 format(1x,a,' - read value = ',a7)
 1014 format(1x,31('*'),' xxdata_04 error ',30('*')//
     &       1x,'Increase configuration string to ',
     &          'more than 5 characters')

C-----------------------------------------------------------------------

      return
      end
