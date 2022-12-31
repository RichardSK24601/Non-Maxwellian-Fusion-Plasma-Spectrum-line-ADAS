       subroutine h9cvrt( adf04       , dsnout  , adf37,
     &                    dtype       , kap_val , dru_val ,
     &                    iadftyp_out ,
     &                    numte       , te
     &                  )

       implicit none
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: H9CVRT *********************
C
C  VERSION:  1.1
C
C  PURPOSE:  Convert type 1 and 5 specific ion files to type 3 or 4.
C
C  CALLING PROGRAM: ADAS809
C
C  SUBROUTINE:
C
C  input : (c*120)adf04      =  type 1 or 5 adf04 dataset for conversion
C  input : (c*120)dsnout     =  type 3 or type 4 adf04 dataset
C  input : (c*120)adf37      =  non-Maxwellian distribution dataset
C  input : (i*4)  dtype      =  0 - Maxwellian, 1 - kappa, 2 - numerical
C                               3 - Druyvesteyn
C  input : (r*8)  kap_val    =  kappa parameter
C  input : (r*8)  dru_val    =  Druyvesteyn x parameter
C  input : (i*4)  iadftyp_out=  type of output adf04 dataset (3 or 4)
C  input : (i*4)  numte      =  number of temperatures for output
C  input : (i*4)  te()       =  temperatures for output (K)
C                                 1st dim: temperature index
C
C
C
C  ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          -------------------------------------------------------------
C          xxdata_04  ADAS      gathers data from adf04 input file
C          xxdata_37  ADAS      gathers data from adf37 input file
C
C  AUTHOR:  Martin O'Mullane
C
C  DATE:    05-09-2010
C
C
C  MODIFICATION HISTORY:
C
C  VERSION : 1.1
C  DATE    : 05-09-2010
C  MODIFIED: Martin O'Mullane
C            - First version.
C
C
C VERSION : 1.2
C DATE    : 02-06-2011
C MODIFIED: Hugh Summers
C           - added parameters iadftyp and lthshft.
C
C VERSION : 1.3
C DATE    : 21-11-2011
C MODIFIED: Hugh Summers
C           - increased ndtrn to 300000.
C           - increase cstrga string length to 31 and propagate
C             variable field length through connected subroutines
C           - switch to revised xxdata_04.for and xxwrto_04.for
C           - removed parameter lthshft since obtained internally
C           - moved comment text lines for reprocesssing into h9cvrt.for
C           - changed iadftyp to iadftyp_out to avoid confusion
c
c  version : 1.4
c  date    : 23-12-2011
c  modified: Hugh Summers
c            - correct problem of degeneracy with type 1 X-parameter for
c              energy on input with a dzero theshold shift in the
c              transformation matching the same inverse transformation
c              in h9ntqd.for
c
c  version : 1.5
c  date    : 25-10-2012
c  modified: Hugh Summers
c            - increased ndtrn dimension to 600000.  Modified warnings
c              to allow ndtrn field length
c
c  version : 1.6
c  date    : 03-12-2012
c  modified: Martin O'Mullane
c            - add 12 elements of headroom to ctext for cases where
c              there are more comments in the adf04 files than ndtext.
c            - Correct degeneracy problems with a dzero shift to the
c              transition energy (evt).
c
c  version : 1.7
c  date    : 18-12-2012
c  modified: Martin O'Mullane
c            - allow orbital from n=8, increased from n=6.
c
c  version : 1.8
c  date    : 25-03-2013
c  modified: Martin O'Mullane
c            - Calling arguments of h9cvrt are changed to allow
c              specification of kappa or Druyvesteyn distributions.
c            - Extend detail of the comments.
c
c  version : 1.9
c  date    : 30-04-2015
c  modified: Hugh Summers
c            - Added counter for type 1-3 ambiguous asymptotic cases.
c            - Trapped zero cross-sections and added counter of cases.
c            - Added warning message for  zero cross-section cases
c            - added info. to comments in output adf04 type 3 dataset.
c
c  version : 1.10
c  date    : 09-11-2015
c  modified: Stuart Henderson
c            - Altered loop value for zero cross-section condition
c            - Change 0,nvt to 0,numte
c
c  version : 1.11
c  date    : 15-08-2017
c  modified: Stuart Henderson
c            - Set las=.false. at beginning of code
c
c  version : 1.12
c  date    : 25-07-2018
c  modified: Martin O'Mullane
c            - Extend length of ctext to hold the maximum length
c              of filename plus leading text.
c            - Update of xxwrto_04 requires lx_untied to be set
c              and to increase ndqdn from 8 to 30.
c
c  version : 1.13
c  date    : 27-07-2018
c  modified: Martin O'Mullane
c            - Extend number of energies in adf37 files from 200 
c              to 1500.
c
c  version : 1.14
c  date    : 04-09-2018
c  modified: Martin O'Mullane
c            - Add section to process S-lines.
c            - Internal variable string should be the same length
c              as ctext().
c
c  version : 1.15
c  date    : 18-02-2020
c  modified: Martin O'Mullane
c            - ndqdn should remain at 8.
c
c  version : 1.16
c  date    : 01-03-2020
c  modified: Martin O'Mullane
c            - Add irestyp argument to xxdata_04.
c
C-----------------------------------------------------------------------
      integer    ndlev        , ndtrn        , ntdim     , ndmet     ,
     &           ndqdn        , nvmax        , nedim     , nfdim     ,
     &           ndtext       , iunit
C-----------------------------------------------------------------------
      real*8     dzero
C-----------------------------------------------------------------------
      parameter( ndlev = 2000 , ndtrn = 600000 ,
     &           ntdim = 50   , ndmet = 4      , ndqdn = 8  ,
     &           nedim = 50   , nfdim = 1500   , nvmax = 50 ,
     &           ndtext= 2000 )
      parameter( iunit = 67   )
      parameter( dzero = 1.0d-30 )
C-----------------------------------------------------------------------
      integer     i4unit
      integer     it           , numte
      integer     iz           , iz0          , iz1       ,
     &            il           , itran        , maxlev     ,
     &            icnte        , icntp        , icntr     , icnth      ,
     &            icnti        , icntl        , icnts     ,
     &            nv           , nvt          ,
     &            iupper       , ilower       , lupper    , llower     ,
     &            iorb         , iadftyp      , itieactn  , irestyp
      integer     ifint        , itype        , itypt     , ilinr      ,
     &            iescl        , ic
      integer     dtype        , n_en         , ntemp
      integer     npl          ,
     &            icateg       , nform1       , nform2    ,
     &            nplr         , npli         , ntext     , iword      ,
     &            nwords
      integer     i            , j            , iadftyp_out            ,
     &            itran_new    , l1           , l2
      integer     nct0         , nct1         , nct2      , nct3
C-----------------------------------------------------------------------
      logical     lexist
      logical     lprn         , lcpl         , lorb      , lbeth
      logical     letyp        , lptyp        , lrtyp     , lhtyp     ,
     &            lityp        , lstyp        , lltyp     , las       ,
     &            lthshft      , lzero        ,
     &            lamb0        , lamb1        , lamb2     , lamb3     ,
     &            lx_untied
C-----------------------------------------------------------------------
      real*8      zeta         , ip, ip_ev,
     &            param1
      real*8      eupper       , elower       , wupper    , wlower     ,
     &            bwno         , evt          , a04
      real*8      kap_val      , dru_val      , blimit
      real*8      ommax
C-----------------------------------------------------------------------
      character   adf04*120    , adf37*120    , dsnout*120   ,
     &            titled*3     , header*80    , string*146
      character   str_today*10 , realname*30
C-----------------------------------------------------------------------
      integer     ia(ndlev)    , isa(ndlev)   , ila(ndlev)   ,
     &            i1a(ndtrn)   , i2a(ndtrn)   , npla(ndlev)  ,
     &            ipla(ndmet,ndlev)
      integer     ietrn(ndtrn) , iptrn(ndtrn) ,
     &            irtrn(ndtrn) , ihtrn(ndtrn) ,
     &            iitrn(ndtrn) , iltrn(ndtrn) ,
     &            istrn(ndtrn) ,
     &            ie1a(ndtrn)  , ie2a(ndtrn)  ,
     &            ip1a(ndtrn)  , ip2a(ndtrn)  ,
     &            ia1a(ndtrn)  , ia2a(ndtrn)  ,
     &            il1a(ndtrn)  , il2a(ndtrn)  ,
     &            is1a(ndtrn)  , is2a(ndtrn)
      integer     ifirst(2)    , ilast(2)
      integer     i1a_new(ndtrn)   , i2a_new(ndtrn)
C-----------------------------------------------------------------------
      real*8      xja(ndlev)   , wa(ndlev)    , aval(ndtrn) , cea(ndtrn)
      real*8      omga(nedim,ndtrn)           , aa(ndtrn)
      real*8      gammaup(ntdim)      , gammadn(ntdim),
     &            scomup(ntdim,ndtrn ), scomdn(ntdim,ndtrn),
     &            en(ntdim,nfdim)     , f(ntdim,nfdim)     ,
     &            utva(ntdim)         , utva2(ntdim)
      real*8      qdorb((ndqdn*(ndqdn+1))/2)  , qdn(ndqdn)   ,
     &            beth(ndtrn)  , bwnoa(ndmet) , prtwta(ndmet),
     &            zpla(ndmet,ndlev)           ,
     &            param2(2)    , median(ntdim) ,
     &            auga(ndtrn)  , wvla(ndlev)         ,
     &            rtions(ntdim), rtrecm(ntdim)
      real*8      scom(nedim)  , scx(nedim)   , oma(nedim)
      real*8      te(ntdim)    , tek(ntdim)   , tva(ntdim)
      real*8      beth_new(ndtrn)    , scom_new(ntdim,ndtrn)    ,
     &            aval_new(ndtrn)
C-----------------------------------------------------------------------
      logical     lss04a(ndlev,ndmet)
      logical     lqdorb((ndqdn*(ndqdn+1))/2)
      logical     ltied(ndlev), lbseta(ndmet)
C-----------------------------------------------------------------------
      character   filnam(nfdim)*120 , filout*120 , calgeb(ntdim,4)*25  ,
     &            ealgeb(ntdim)*25
      character   tcode(ndtrn)*1 , pecode(ndtrn)*1 , tecode(ndtrn)*1
      character   cstrga(ndlev)*64, cpla(ndlev)*1  , cprta(ndmet)*9
      character   ctext(ndtext+13)*146
      character   tcode_new(ndtrn)*1
C-----------------------------------------------------------------------

      lthshft = .false.
      las     = .false.

      if (numte.GT.ntdim) then
         write(i4unit(-1),1001)
     &          'Maximum number of temperatures allowed : ', ntdim
         write(i4unit(-1),1002)
         stop
      endif

C----------------------------------------------------------------------
C zero abnormal cross-section type counts and logicals
C----------------------------------------------------------------------

      nct0 = 0
      nct1 = 0
      nct2 = 0
      nct3 = 0

      lamb0 = .false.
      lamb1 = .false.
      lamb2 = .false.
      lamb3 = .false.

C-----------------------------------------------------------------------
C Read specific ion file and electron distribution file data if required.
C Note that only numerical distributions and not superposition
C descriptions are enabled.
C-----------------------------------------------------------------------

      inquire( file=adf04, exist=lexist)
      if (lexist) then
         open( unit=iunit , file=adf04  , status='old' )
      else
         write(i4unit(-1),1001)
     &           'Input adf04 dataset does not exist: '//adf04
         write(i4unit(-1),1002)
         stop
      endif

      itieactn = 1
      call xxdata_04( iunit  ,
     &                ndlev  , ndtrn  , ndmet   , ndqdn , nvmax ,
     &                titled , iz     , iz0     , iz1   , bwno  ,
     &                npl    , bwnoa  , lbseta  , prtwta, cprta ,
     &                il     , qdorb  , lqdorb  , qdn   , iorb  ,
     &                ia     , cstrga , isa     , ila   , xja   ,
     &                wa     ,
     &                cpla   , npla   , ipla    , zpla  ,
     &                nv     , scom   ,
     &                itran  , maxlev ,
     &                tcode  , i1a    , i2a     , aval  , omga  ,
     &                beth   ,
     &                iadftyp, lprn   , lcpl    , lorb  , lbeth ,
     &                letyp  , lptyp  , lrtyp   , lhtyp , lityp ,
     &                lstyp  , lltyp  , itieactn, ltied , irestyp
     &              )

      close(iunit)

      if (adf37(1:4).NE.'NULL') then

         dtype = 2
         inquire( file=adf37, exist=lexist)
         if (lexist) then
            open( unit=iunit , file=adf37  , status='old' )
         else
            write(i4unit(-1),1001)
     &          'Input adf37 dataset does not exist: '//adf37
            write(i4unit(-1),1002)
            stop
         endif

         call xxdata_37( iunit  ,
     &                   nfdim  , nvmax  ,
     &                   header , icateg , n_en   , nvt    ,
     &                   nform1 , param1 , nform2 , param2 ,
     &                   en     , f      , utva   , utva2  ,
     &                   median , filnam , filout , calgeb ,
     &                   ealgeb
     &                 )

        if (icateg.ne.2) then
          write(i4unit(-1),1001)
     &        'Only numerical distributions allowed for now'
          write(i4unit(-1),1002)
          stop
        endif

        write(i4unit(-1),1000)'Using temperatures from adf37 file'
        write(i4unit(-1),1003)

        numte = nvt
        do it = 1, nvt
           tva(it) = utva(it)
           tek(it) = utva(it) * 11605.0D0
        end do

        close(iunit)

      else

         do it = 1, numte
            tek(it) = te(it)
            tva(it) = te(it) / 11605.0D0
         end do

      endif


C-----------------------------------------------------------------------
C Acquire comments from adf04 file
C-----------------------------------------------------------------------

      open(unit=iunit, file=adf04, status='old')
      call xxcomm(iunit, ndtext, ctext, ntext)
      close(iunit)

      do i = 1, ntext

        iword  = 1
        j      = 2
        call xxcase(ctext(i),string,'lc')
        call xxword(string, ' '      , iword    , j ,
     &              ifirst(1)   , ilast(1) , nwords
     &             )
        if (nwords .ge. 2) then
          if(string(ifirst(1):ilast(1)).eq.'c'.and.
     &       string(ifirst(2):ilast(2)).eq.'autostructure')
     &       las = .TRUE.
        endif
      end do

C-----------------------------------------------------------------------
C Remap collision strengths for AUTOSTRUCTURE non-neutral type 1 data
C-----------------------------------------------------------------------

      if (las.and.iz1.ne.1.and.iadftyp.eq.1) then
        call remapx(ndtrn, nvmax, nv, itran, scom, omga)
        lthshft = .true.
        write(i4unit(-1),1000)'Using zero-threshold correction'
      endif

C-----------------------------------------------------------------------
C Sort transitions into transition/recombination types.
C-----------------------------------------------------------------------

      call bxttyp( ndlev  , ndmet  , ndtrn  , nplr  , npli  ,
     &             itran  , tcode  , i1a    , i2a   , aval  ,
     &             icnte  , icntp  , icntr  , icnth , icnti ,
     &             icntl  , icnts  ,
     &             ietrn  , iptrn  , irtrn  , ihtrn , iitrn ,
     &             iltrn  , istrn  ,
     &                               ie1a   , ie2a  , aa    ,
     &                               ip1a   , ip2a  ,
     &                               ia1a   , ia2a  , auga  ,
     &                               il1a   , il2a  , wvla  ,
     &                               is1a   , is2a  , lss04a
     &           )

C And determine Burgess-Tully electron impact excitation types

      call bfttyp(ndlev  , ndtrn  ,
     &            iz1    , il     ,
     &            ia     , cstrga , isa    , ila   , xja   , wa ,
     &            itran  , tcode  , i1a    , i2a   , aval  ,
     &            icnte  , icntp  , icntr  , icnth , icnti ,
     &            ietrn  , pecode , tecode , ie1a  , ie2a  , aa ,
     &            cea
     &          )

C-----------------------------------------------------------------------
C  Set standard parameters for the Maxwell averaging
C-----------------------------------------------------------------------

      ifint = 1
      ilinr = 1
      iescl = 1

      if (iz1.eq.1) then
        itypt = 2
      else
        itypt = 1
      endif

C-----------------------------------------------------------------------
C  Cycle through electron excitation transition in turn
C-----------------------------------------------------------------------

      do it = 1, nv
        scx(it) = scom(it)
      end do

      do ic = 1, icnte

        if (tcode(ietrn(ic)).ne.' ') then
           read(tcode(ietrn(ic)),'(i1)')itype
        else
           read(tecode(ietrn(ic)),'(i1)')itype
           if (beth(ietrn(ic)).lt.0.0D0) then
              itype = 1
           elseif (beth(ietrn(ic)).gt.0.0D0) then
              itype = 2
           endif
        endif

        blimit = beth(ic)

        call h9tran( ndlev  , ndtrn     , nedim  ,
     &               il     , ietrn(ic) , nv     ,
     &               ia     , wa        , xja    ,
     &               i1a    , i2a       , aval   ,
     &               omga   ,
     &               iupper , ilower    ,
     &               lupper , llower    ,
     &               wupper , wlower    ,
     &               eupper , elower    ,
     &               a04    , oma
     &             )

        evt = 13.6048d0 * (eupper-elower+dzero) / 109737.26d0

C Convert type 5 scaled energy to x-parameter

        if (iadftyp.eq.5) then

           do it = 1, nv
             scx(it) = 13.6048D0 * scom(it) / (evt+dzero) + 1.0D0
           end do

        endif

C----------------------------------------------------------------------
C Trap zero cross-sections
C----------------------------------------------------------------------

        lzero = .false.
        ommax = 0.0d0
        do it = 1,nv
          ommax = dmax1(ommax,oma(it))
        enddo

        if (ommax.le.dzero) then
           lzero=.true.
           if (nct0.eq.0) then
              write(i4unit(-1),1000)'adf04 dataset has zero xsects.'
              write(i4unit(-1),1003)
           endif
           do i = 1,numte
             gammaup(i)=1.0d-30
             gammadn(i)=1.0d-30
           enddo
           lamb0 = .true.
           nct0 = nct0+1
        endif

        if (.not.lzero) then
           call h9ntqd( nedim   , ntdim   , nfdim    ,
     &                  ifint   , itype   , itypt    , ilinr , iescl ,
     &                  ietrn(ic) ,
     &                  nv      , numte   ,
     &                  evt     , scx     ,
     &                  oma     , tva     ,
     &                  kap_val , dru_val , dtype    ,
     &                  n_en    , en      , f        ,
     &                  nvt     , utva    ,
     &                  lbeth   , blimit  ,
     &                  nform1  , param1  , nform2   , param2 ,
     &                  gammaup , gammadn ,
     &                  lamb1   , lamb2   , lamb3     ,
     &                  nct1    , nct2    , nct3
     &                )
        endif

        do it = 1, numte
          scomup(it,ietrn(ic)) = gammaup(it)
          scomdn(it,ietrn(ic)) = gammadn(it)
        enddo

      enddo

C----------------------------------------------------------------------
C  Cycle through ionisation S-lines
C----------------------------------------------------------------------

      do ic = 1,icnts

	 call h9trni( ndlev  , ndtrn     , nedim  , ndmet ,
     &                il     , istrn(ic) , nv     ,
     &                ia     , wa        , xja    ,
     &                i1a    , i2a       , aval   , 
     &                omga   , zpla      , bwnoa  , ipla ,
     &                iupper , ilower    , 
     &                lupper , llower    ,
     &                wupper , wlower    ,
     &                eupper , elower    ,
     &                aa     , oma       ,
     &                zeta   , ip
     &              )
	
	 ip_ev = 13.6048d0*ip/109737.26d0
	                 
         call h9qd3b( nedim   , ntdim   , nfdim  ,
     &                wupper  ,
     &                nv      , numte   ,
     &                ip_ev   , scx     ,
     &                oma     , tva     ,
     &                kap_val , dru_val , dtype  ,
     &                n_en    , en      , f      ,
     &                nvt     , utva    ,
     &                nform1  , param1  , nform2 , param2 ,
     &                rtions  , rtrecm
     &              )

	 do it = 1, numte
           scomup(it,istrn(ic)) = rtions(it)
           scomdn(it,istrn(ic)) = rtrecm(it) 
         enddo

      enddo

C----------------------------------------------------------------------
C Write output file and add on original ctext with additions
C----------------------------------------------------------------------

      open( unit=iunit , file=dsnout)

      do i = ntext,1,-1
         ctext(i+2)=ctext(i)
      end do

      ntext=ntext+2
      write(ctext(1), '("c",79("-"))')
      write(ctext(2), '(a)')'c comments from original adf04 file'

      call xxtoday(str_today)
      call xxname(realname)

      write(ctext(ntext+1),'("c",79("-"))')
      write(ctext(ntext+2), '("c")')
      write(ctext(ntext+3), '(A,I1,A,I1)')'c  Conversion of an '//
     &                'adf04 dataset from type ',iadftyp,
     &                ' to type ',iadftyp_out
      write(ctext(ntext+4), '("c")')
      write(ctext(ntext+5), '(a,l1)')
     &     "c  Parameter: Cowan ion threshold shift = ",lthshft

      if (dtype.eq.0) then
         write(ctext(ntext+6), '("c")')
      elseif (dtype.eq.1) then
         write(ctext(ntext+6), '(a,1p,e12.4)')
     &     "c  kappa distribution : ", kap_val
      elseif (dtype.eq.2) then
         call xxslen(adf37, L1, L2)
         write(ctext(ntext+6), '(a,a)')
     &     "c  adf37 distribution : ", adf37(L1:L2)
      elseif (dtype.eq.3) then
         write(ctext(ntext+6), '(a,1p,e12.4)')
     &     "c  Druyvesteyn distribution : ", dru_val
      endif

      write(ctext(ntext+7), '("c")')
      write(ctext(ntext+8), '("c  Producer : ", A)')realname
      write(ctext(ntext+9), '("c  Date     : ", A)')str_today
      write(ctext(ntext+10), '("c")')
      write(ctext(ntext+11), '("c",79("-"))')

      ntext = ntext+11

      itran_new = 0

C excitation

      do i=1,icnte
        tcode_new(itran_new+1)=tcode(ietrn(i))
        i1a_new(itran_new+1)=i1a(ietrn(i))
        i2a_new(itran_new+1)=i2a(ietrn(i))
        aval_new(itran_new+1)=aval(ietrn(i))
        do it=1,numte
          scom_new(it,itran_new+1)=scomup(it,ietrn(i))
        enddo
        if(iadftyp_out.eq.4) then
            tcode_new(itran_new+2)=tcode(ietrn(i))
            i1a_new(itran_new+2)=i2a(ietrn(i))
            i2a_new(itran_new+2)=i1a(ietrn(i))
            aval_new(itran_new+2)=0.0d0
            do it=1,numte
              scom_new(it,itran_new+2)=scomdn(it,ietrn(i))
            enddo
            itran_new = itran_new+2
        elseif(iadftyp_out.eq.3) then
            itran_new = itran_new+1
        else
            write(i4unit(-1),1000)'Output adf04 type ', iadftyp_out,
     &                            ' not allowed.'
            write(i4unit(-1),1003)
            return
        endif
      enddo


C ionisation

      do i=1,icnts
        
        tcode_new(itran_new+1) = tcode(istrn(i))
        i1a_new(itran_new+1)   = i1a(istrn(i))
        i2a_new(itran_new+1)   = i2a(istrn(i))
        aval_new(itran_new+1)  = 0.0
        
        do it=1,numte
          scom_new(it,itran_new+1)=scomup(it,istrn(i))
        enddo
        
        if(iadftyp_out.eq.4) then
            tcode_new(itran_new+2) = tcode(istrn(i))
            i1a_new(itran_new+2)   = i2a(istrn(i))
            i2a_new(itran_new+2)   = i1a(istrn(i))
            aval_new(itran_new+2)  = 0.0d0
            
            do it=1,numte
              scom_new(it,itran_new+2)=scomdn(it,istrn(i))
            enddo
            
            itran_new = itran_new+2
        
        elseif(iadftyp_out.eq.3) then
            itran_new = itran_new+1
        else
            write(i4unit(-1),1000)'Output adf04 type ', iadftyp_out,
     &                            ' not allowed.'
            write(i4unit(-1),1003)
            return
        endif
        
      enddo

C write out adf04 file

      lx_untied = .TRUE.
      call xxwrto_04( iunit  ,
     &                ndlev  , ndtrn  , ntdim , ndmet  , ndqdn ,
     &                ndtext ,
     &                iz     , iz0    , iz1   ,
     &                npl    , bwnoa  , cprta ,
     &                il     , qdorb  , lqdorb, qdn    , iorb  ,
     &                cstrga , isa    , ila   , xja    , wa    ,
     &                cpla   , npla   , ipla  , zpla   ,
     &                numte  , tek    ,
     &                itran_new       ,
     &                tcode_new       , i1a_new        , i2a_new   ,
     &                aval_new        , scom_new       ,
     &                iadftyp_out     ,
     &                lbeth  , beth   ,
     &                ntext  , ctext  ,
     &                lx_untied
     &              )


      close(iunit)

C-----------------------------------------------------------------------

 1000 format(1x,31('*'),' h9cvrt warning ',31('*'),/,
     &       a,i7,a)
 1001 format(1x,31('*'),'  h9cvrt error  ',31('*'),/,
     &       a,i7,a)
 1002 format(1x,29('*'),' program terminated ',29('*'))
 1003 format(1x,29('*'),' program continues  ',29('*'))
C-----------------------------------------------------------------------

      return
      end
