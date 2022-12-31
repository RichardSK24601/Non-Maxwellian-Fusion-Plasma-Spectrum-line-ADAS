       subroutine h9ntqd  ( nedim   , ntdim   , nfdim  ,
     &                      ifint   , itype   , itypt  , ilinr , iescl ,
     &                      itrn    ,
     &                      ne      , nt      ,
     &                      evt     , xa      ,
     &                      oma     , tva     ,
     &                      kappa   , dru_val , dist   ,
     &                      nef     , en      , f      ,
     &                      ntf     , tvf     ,
     &                      lbeth   , beth    ,
     &                      nform1  , param1  , nform2 , param2 ,
     &                      upsilon , dnsilon ,
     &                      lamb1   , lamb2   , lamb3  ,
     &                      nct1    , nct2    , nct3
     &                     )
      implicit none
C-----------------------------------------------------------------------
C
C  VERSION:   1.0
C
C  PURPOSE:  Executes quadratures over collision strengths to form
C            excitation and de-excitation effective collision strengths
C            for atoms and ions with tabulated collision stregths
C            as a function of x parameter.
C
C            Quadrature can be executed over a Maxwellian, kappa
C            distribution, Druyvesteyn or numerical distribution.
C            Linear interpolation is recommended and is default (ilinr=1)
C            Quadratic interpolation is also allowed for analytic
C            distributions.
C            1/E variable interpolation is allowed for Maxwellian only.
C
C
C  INPUT and OUTPUT:
C
C          (i*4)  nedim       = input  = max no of energies in omega file
C          (i*4)  ntdim       = input  = max no of temperatures
C          (i*4)  nfdim       = input  = max no of energies in adf37 file
C          (i*4)  ifint       = input  = indep. var. for interpolation
C                                           (1 = E)
C                                           (2 = 1/E)
C          (i*4)  itype       = input  = collision strength type, to give
C                                        high energy behaviour
C                                           (1 = dipole --> a*log(X+b) )
C                                           (2 = non-dp --> a+X/b      )
C                                           (3 = spin ch--> a/(X+b)**2 )
C          (i*4)  itypt       = input  = threshold behaviours allowed
C                                           (1(ion)     = const to 1st pt.)
C                                           (2(neutral) = 0     to 1st pt.)
C          (i*4)  ilinr       = input  = allow linear or quadratic interp
C                                           (1 = linear   )
C                                           (2 = quadratic)
C          (i*4)  iescl       = input  = allow e**2*omega + lin. interp
C                                           (1 = normal use)
C                                           (2 = e**2*omega +lin.)
C                                        iescl=2 not implimented
C          (i*4)  itrn        = input  = index of current transition
C          (i*4)  ne          = input  = number of energies in omega file
C          (i*4)  nt          = input  = number of temperatures
C          (r*8)  evt         = input  = theshold energy (eV)
C          (r*8)  xa()        = input  = tabul. x param. for coll. str.
C          (r*8)  oma()       = input  = tabul. coll. str.
C          (r*8)  tva()       = input  = temperatures (eV)
C          (r*8)  kappa       = input  = kappa value of electron dist.
C          (r*8)  dru_val     = input  = x parameter from Druyvesteyn dist.
C          (i*4)  dist        = input  = electron distribution
C                                           (0 = Maxwellian )
C                                           (1 = kappa      )
C                                           (2 = numerical  )
C                                           (3 = Druyvesteyn)
C          (i*4)  nef         = input  = number of energies in adf37 file
C          (r*8)  en(,)       = input  = adf37 energy (eV)
C          (r*8)  f(,)        = input  = adf37 distribution
C          (i*4)  ntf         = input  = number of temperatures in adf37
C          (r*8)  tvf()       = input  = temperatures (eV) from adf37
C          (l*4)  lbeth       = input  = true if limit point exists
C          (r*8)  beth        = input  = infinite energy limit point of omega
C                                        (note this is a signed vriable: 
C                                         -ve => dipole, +ve => non-dipole or
C                                         spin change)
C          (i*4)  nform1      = input  = type of threshold behaviour
C                                          (1 = cutoff       )
C                                          (2 = energy^param1)
C          (r*8)  param1      = input  = parameter of threshold form
C          (i*4)  nform2      = input  = type of high-energy behaviour
C                                          (1 => cutoff                )
C                                          (2 => energy^-param2(1)     )
C                                          (3 => exp(-param2(1)*energy))
C          (r*8)  param2()    = input  = parameter of high-energy form
C
C          (r*8)  upsilon(,)  = output = upsilon values
C          (r*8)  dnsilon(,)  = output = downsilon values
C          (l*1)  lamb1       = i/o    = .true.=> type 1 asymp. ambiguity
C                                        .false. => no type 1 ambig. cases 
C          (l*1)  lamb2       = i/o    = .true.=> type 2 asymp. ambiguity
C                                        .false. => no type 2 ambig. cases 
C          (l*1)  lamb3       = i/o    = .true.=> type 3 asymp. ambiguity
C                                        .false. => no type 3 ambig. cases 
C          (i*4)  nct1        = i/o    = type 1 ambiguous case count
C          (i*4)  nct2        = i/o    = type 2 ambiguous case count
C          (i*4)  nct3        = i/o    = type 3 ambiguous case count
C
C
C  PROGRAM:
C
C          (r*8)  xf      = program   = current en(i)/evt
C          (r*8)  omega(,)= program   = oma interpolated to distribution
C                                       function energy grid
C          (r*8)  sumi()  = program   = gamma contrib. from i    -> i+1
C          (r*8)  sumn()  = program   = gamma contrib. from ne-1  -> ne
C          (r*8)  sumu()  = program   = gamma contrib. from ne --> inf.
C          (r*8)  suml()  = program   = gamma contrib. from thres. -> 1
C          (r*8)  en()    = program   = tabul. ener. for coll. str. (ev)
C          (r*8)  fva()   = program   = indep. var. for interpolation
C          (r*8)  expi()  = program   = current exp(-(ui-ut))
C          (r*8)  expi1() = program   = current exp(-(ui1-ut))
C          (r*8)  exp1()  = program   = exp(-(u1-ut))
C          (r*8)  ui()    = program   = current eva(i)/kte
C          (r*8)  ui1()   = program   = current eva(i+1)/kte
C          (r*8)  u1()    = program   = eva(1)/kte
C          (r*8)  ut      = program   = evt/kte
C          (r*8)  uj()    = program   = ui-ut
C          (r*8)  uj1()   = program   = ui1-ut
C          (r*8)  w0      = program   = interpolation working variable
C          (r*8)  w1      = program   = interpolation working variable
C          (r*8)  w2      = program   = interpolation working variable
C          (r*8)  v0      = program   = interpolation working variable
C          (r*8)  v1      = program   = interpolation working variable
C          (r*8)  v2      = program   = interpolation working variable
C          (r*8)  y1      = program   = interpolation working variable
C          (r*8)  y2      = program   = interpolation working variable
C          (r*8)  c0      = program   = interpolation working variable
C          (r*8)  c1      = program   = interpolation working variable
C          (r*8)  c2      = program   = interpolation working variable
C          (r*8)  cc0     = program   = interpolation working variable
C          (r*8)  cc1     = program   = interpolation working variable
C          (r*8)  ww0     = program   = interpolation working variable
C          (r*8)  ww1     = program   = interpolation working variable
C          (r*8)  ww2     = program   = interpolation working variable
C          (r*8)  a1      = program   = interpolation working variable
C          (r*8)  a2      = program   = interpolation working variable
C          (r*8)  b1      = program   = interpolation working variable
C          (r*8)  b2      = program   = interpolation working variable
C
C
C  routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          eei        copase    evaluates exp(x)*E1(x)
C          ee2        copase    evaluates exp(x)*E2(x)
C          lngama               evaluates ln(gamma(a))
c          ingama               evaluates incomplete gamma P(a,x)
c          ingamq               evaluates incomplete gamma 1-P(a,x)
C
C author:  H P Summers
C          K1/1/57
C          JET ext. 4941
C
C date:    26/05/93
C
C update:  30/11/01  HP Summers - altered input to use x parameter
C
C update:  23/11/04  P Bryans - altered to evaluate non-maxwellian
C                               electron distributions
C
C update:  20/07/07  A Whiteford - Modified comments slightly to allow
C                                  for automatic generation of
C                                  documentation.
C
C update:  29/05/09  A Whiteford - Fixed bug in indexing of last
C                                  element of xa() array.
C
C  VERSION : 1.6                     
C  DATE    : 07-09-2010
C  MODIFIED: Martin O'Mullane
C            - Move upsilon and dnsilon, the only outputs, to the end
C              of the argument list.
C            - Re-order documentation to follow argument listing.
c
c  version : 1.7
c  date    : 07-10-2011
c  modified: Hugh Summers
c            - negative sign of beth parameter not taken note of in 
c              formulae involving type 1 transitions.  Corrected by
c              introducing bthp = dabs(beth) internally.  
c
c  version : 1.8
c  date    : 23-12-2011
c  modified: Hugh Summers
c            - correct problem of degeneracy with type 1 X-parameter for
c              energy on input with a dzero theshold shift in the 
c              transformation matching the same inverse transformation
c              in h9cvrt.for  
c
c  version : 1.9
c  date    : 05-12-2012
c  modified: Martin O'Mullane
c            - Remove ambiguity in the error warning.
c
c  version : 1.10
c  date    : 20-04-2013
c  modified: Martin O'Mullane
c            - Trap divide by zero possibility before determining 
c              step size.
c
c  version : 1.11
c  date    : 30-04-2015
c  modified: Hugh Summers
c            - Added counter for type 1-3 ambiguous asymptotic cases.
c            - Modified warning message output for ambiguous cases.
c
C-----------------------------------------------------------------------
      integer     max
C-----------------------------------------------------------------------
      real*8      tol   , dzero
C-----------------------------------------------------------------------
      parameter ( max = 1000)
C-----------------------------------------------------------------------
      parameter ( tol = 1.0D-10 , dzero = 1.0d-30 )
C-----------------------------------------------------------------------
      integer     nedim  , ntdim   , nfdim   , iter   , itrn   , i4unit
      integer     ne     , nt      , itype   , itypt  , ifint  , ilinr
      integer     iescl  , i       , j       , it     , dist   , ie
      integer     ntf    , nef     , nform1  , nform2
      integer     nct1   , nct2    , nct3 
C-----------------------------------------------------------------------
      logical     lbeth  , lhighe
      logical     lamb1  , lamb2   , lamb3
C-----------------------------------------------------------------------
      real*8      pi     , eei     , ee2     , lngama , ingama , ingamq
      real*8      y0     , y1      , k1      , k2     , a      , b
      real*8      v0     , v1      , v2
      real*8      w0     , w1      , w2
      real*8      ww0    , ww1     , ww2
      real*8      c0     , c1      , c2
      real*8      a0     , a1      , b0      , b1
      real*8      beth   , dru_val , param1  , bthp  
      real*8      uiu    , uiu1    , uju     , uju1
      real*8      intgl1 , intgl2  , evt
      real*8      estep  , kappa   , xf
      real*8      igam1i , igam2i  , igam3i  , igam1j , igam2j , igam3j
      real*8      igmxi1 , igmxi2  , igmxj1  , igmxj2
      real*8      intgl_tot1       , intgl_tot2
C-----------------------------------------------------------------------
      real*8      sumu(2), suml(2) , sumi(2) , sumn(2), sumx(2)
      real*8      param2(2)
      real*8      xa(nedim)   , oma(nedim)   , eva(nedim)
      real*8      fva(nedim)  , expi(ntdim)  , expi1(ntdim)
      real*8      tva(ntdim)  , tvf(ntdim)   , sum(2,ntdim)
      real*8      ut(ntdim)   , ui(ntdim)    , ui1(ntdim)
      real*8      uj(ntdim)   , uj1(ntdim)
      real*8      u1(ntdim)   , exp1(ntdim)
      real*8      upsilon(ntdim)      , dnsilon(ntdim)
      real*8      en(ntdim,nfdim)     , eni(ntdim,nfdim)
      real*8      f(ntdim,nfdim)
      real*8      fj(ntdim,nfdim)     , fi(ntdim,nfdim)
      real*8      omegai(ntdim,nfdim) , omegaj(ntdim,nfdim)
      real*8      atest1, atest2
C-----------------------------------------------------------------------

      lhighe=.false.
      
      if(lbeth) bthp=dabs(beth)

C-----------------------------------------------------------------------
C  zero gama vector
C-----------------------------------------------------------------------

      pi = acos(-1.0d0)

      if (dist.ne.2) then
        do it = 1,nt
          sum(1,it)=0.d0
          sum(2,it)=0.d0
        enddo
      else
        do it = 1,ntf
          sum(1,it)=0.d0
          sum(2,it)=0.d0
        enddo
      endif

C-----------------------------------------------------------------------
C  convert x parameter to incident energy in ev
C-----------------------------------------------------------------------

       do i = 1,ne
         eva(i)=xa(i)*(evt+dzero)
c         eva(i)=xa(i)*evt
       enddo

C-----------------------------------------------------------------------
C  analytic distributions
C-----------------------------------------------------------------------

      if (dist.ne.2) then

C-----------------------------------------------------------------------
C  set up initial values at each temperature
C-----------------------------------------------------------------------

      do it = 1,nt

        if     (dist .eq. 0) then
          ut(it) = evt/tva(it)
          ui(it) = eva(1)/tva(it)
          if ((ui(it)-ut(it)) .gt. 165.d0) then
            expi(it)=0.d0
          else
            expi(it)=dexp(-ui(it)+ut(it))
          endif
        elseif (dist .eq. 1) then
          ut(it) = 2.0d0*evt/((2.0d0*kappa-3.0d0)*tva(it))
          ui(it) = 2.0d0*eva(1)/((2.0d0*kappa-3.0d0)*tva(it))
        elseif (dist .eq. 3) then
          ut(it) = (2.d0*evt/(3.d0*tva(it)) *
     &             dexp(lngama(5.d0/(2.d0*dru_val))-
     &             lngama(3.d0/(2.d0*dru_val))))
          ui(it) = (2.d0*eva(1)/(3.d0*tva(it)) *
     &             dexp(lngama(5.d0/(2.d0*dru_val))-
     &             lngama(3.d0/(2.d0*dru_val))))
        endif
        uj(it) = ui(it)-ut(it)
        u1(it) = ui(it)
        exp1(it) = expi(it)

C-----------------------------------------------------------------------
C  set up interpolation independent variable
C-----------------------------------------------------------------------

        do i = 1 , ne
          if     (ifint .eq. 1 .or. iescl .ne. 1) then
            fva(i) = eva(i)
          elseif (ifint .eq. 2) then
            fva(i) = 1.d0/eva(i)
          else
            stop
          endif
        enddo

      enddo

C-----------------------------------------------------------------------
C  loop over energy points
C  approximate collision strength as first order power series
C-----------------------------------------------------------------------

      if (ilinr.eq.1) then

        do 80 i = 1,ne-1

          v0 = 1.d0/(fva(i)-fva(i+1))
          v1 = 1.d0/(fva(i+1)-fva(i))

          w0 = oma(i)*v0 + oma(i+1)*v1
          w1 = -fva(i+1)*oma(i)*v0 -fva(i)*oma(i+1)*v1


C-----------------------------------------------------------------------
C  temperature loop
C-----------------------------------------------------------------------


          do 70 it = 1,nt
            if     (dist .eq. 0) then
              ui1(it) = eva(i+1)/tva(it)
              expi1(it)=dexp(-ui1(it)+ut(it))
            elseif (dist .eq. 1) then
              ui1(it) = 2.d0*eva(i+1)/((2.d0*kappa-3.d0)*tva(it))
            elseif (dist .eq. 3) then
              ui1(it) = (2.d0*eva(i+1)/(3.d0*tva(it)) *
     &                  dexp(lngama(5.d0/(2.d0*dru_val))-
     &                  lngama(3.d0/(2.d0*dru_val))))
            endif
            uj1(it) = ui1(it) - ut(it)

C-----------------------------------------------------------------------
C  add in component part of the rate parameter in this interval
C-----------------------------------------------------------------------

            if (iescl .eq. 1) then

              if (dist .eq. 0) then
                sumi(1) =   expi(it)*(tva(it)*w0*ui(it) + tva(it)*w0+w1)
     &                  - expi1(it)*(tva(it)*w0*ui1(it) + tva(it)*w0+w1)
              elseif (dist .eq. 1) then
                ww0 = w0 * ((2.0d0*kappa-3.0d0)*tva(it)/2.d0)
                ww1 = w1
                sumi(1) =  (1.d0/kappa*(1.d0+ui(it))**(-kappa) *
     &                         ( ww0*ui(it)+ww1 ) +
     &                         1.d0/(kappa*(kappa-1.d0))*(1.d0+ui(it))**
     &                         (1.d0-kappa)* ww0)
     &                     - (1.d0/kappa*(1.d0+ui1(it))**(-kappa) *
     &                         ( ww0*ui1(it)+ww1 ) +
     &                         1.d0/(kappa*(kappa-1.d0))*(1.d0+ui1(it))
     &                         **(1.d0-kappa)* ww0)
                sumi(2) =  (1.d0/kappa*(1.d0+uj(it))**(-kappa) *
     &                         ( ww0*ui(it)+ww1 ) +
     &                         1.d0/(kappa*(kappa-1.d0))*(1.d0+uj(it))**
     &                         (1.d0-kappa)* ww0)
     &                     - (1.d0/kappa*(1.d0+uj1(it))**(-kappa) *
     &                         ( ww0*ui1(it)+ww1 ) +
     &                         1.d0/(kappa*(kappa-1.d0))*(1.d0+uj1(it))
     &                         **(1.d0-kappa)* ww0)

              elseif (dist.eq.3) then
                ww0 = w0 * (3.d0*tva(it) / (2.d0*
     &                dexp(lngama(5.d0/(2.d0*dru_val))-
     &                lngama(3.d0/(2.d0*dru_val)))))
                ww1 = w1
                if (ui1(it)**dru_val.lt.1.d0/dru_val+1d0) then
                  igam1i = dexp(lngama(1.d0/dru_val))*
     &                     (ingama(1.d0/dru_val,ui1(it)**dru_val) -
     &                     ingama(1.d0/dru_val,ui(it)**dru_val))
                else
                  igam1i = dexp(lngama(1.d0/dru_val))*
     &                     (ingamq(1.d0/dru_val,ui(it)**dru_val) -
     &                     ingamq(1.d0/dru_val,ui1(it)**dru_val))
                endif
                if (ui1(it)**dru_val.lt.2.d0/dru_val+1d0) then
                  igam2i = dexp(lngama(2.d0/dru_val))*
     &                     (ingama(2.d0/dru_val,ui1(it)**dru_val) -
     &                     ingama(2.d0/dru_val,ui(it)**dru_val))
                else
                  igam2i = dexp(lngama(2.d0/dru_val))*
     &                     (ingamq(2.d0/dru_val,ui(it)**dru_val) -
     &                     ingamq(2.d0/dru_val,ui1(it)**dru_val))
                endif
                if (uj1(it)**dru_val.lt.1.d0/dru_val+1d0) then
                  igam1j = dexp(lngama(1.d0/dru_val))*
     &                     (ingama(1.d0/dru_val,uj1(it)**dru_val) -
     &                     ingama(1.d0/dru_val,uj(it)**dru_val))
                else
                  igam1j = dexp(lngama(1.d0/dru_val))*
     &                     (ingamq(1.d0/dru_val,uj(it)**dru_val) -
     &                     ingamq(1.d0/dru_val,uj1(it)**dru_val))
                endif
                if (uj1(it)**dru_val.lt.2.d0/dru_val+1d0) then
                  igam2j = dexp(lngama(2.d0/dru_val))*
     &                     (ingama(2.d0/dru_val,uj1(it)**dru_val) -
     &                     ingama(2.d0/dru_val,uj(it)**dru_val))
                else
                  igam2j = dexp(lngama(2.d0/dru_val))*
     &                     (ingamq(2.d0/dru_val,uj(it)**dru_val) -
     &                     ingamq(2.d0/dru_val,uj1(it)**dru_val))
                endif
                sumi(1) = (ww0*igam2i + ww1*igam1i)
                sumi(2) = ww0*igam2j + (ww0*ut(it)+ww1)*igam1j
              endif
            endif

            if (sumi(1).lt.0) sumi(1) = 0.d0
            if (sumi(2).lt.0) sumi(2) = 0.d0
            sum(1,it)=sum(1,it)+sumi(1)
            sum(2,it)=sum(2,it)+sumi(2)


   70     continue

C-----------------------------------------------------------------------
C  move reference point from i to i+1 for next step of quadrature
C-----------------------------------------------------------------------

          do it=1,nt
            ui(it) = ui1(it)
            uj(it) = uj1(it)
            expi(it) = expi1(it)
          end do

   80   continue
   
c        write(i4unit(-1),'(a,1p15e9.2)')'sum(1,*)=',(sum(1,it),it=1,nt)
c	 write(i4unit(-1),'(a,1p15e9.2)')'sum(2,*)=',(sum(2,it),it=1,nt)

C-----------------------------------------------------------------------
C  loop over energy points.
C  approximate collision strength as power series to second order.
C  set up interpolation parameters
C-----------------------------------------------------------------------
C
      else

        do 40 i = 1,ne-2

          v0 = 1.d0/((fva(i)-fva(i+1))*(fva(i)-fva(i+2)))
          v1 = 1.d0/((fva(i+1)-fva(i))*(fva(i+1)-fva(i+2)))
          v2 = 1.d0/((fva(i+2)-fva(i))*(fva(i+2)-fva(i+1)))

          w0= oma(i)*v0 + oma(i+1)*v1 + oma(i+2)*v2
          w1=  -(fva(i+1)+fva(i+2))*oma(i)*v0
     &         -(fva(i)+fva(i+2))*oma(i+1)*v1
     &         -(fva(i)+fva(i+1))*oma(i+2)*v2
          w2=    fva(i+1)*fva(i+2)*oma(i)*v0
     &         + fva(i)*fva(i+2)*oma(i+1)*v1
     &         + fva(i)*fva(i+1)*oma(i+2)*v2
C
C-----------------------------------------------------------------------
C  temperature loop
C-----------------------------------------------------------------------
C
          do 30 it = 1,nt

            if     (dist .eq. 0) then
              ui1(it) = eva(i+1)/tva(it)
              expi1(it)=dexp(-ui1(it)+ut(it))
            elseif (dist .eq. 1) then
              ui1(it) = 2.d0*eva(i+1)/((2.d0*kappa-3.d0)*tva(it))
            elseif (dist .eq. 2) then
              ui1(it) = eva(i+1)/tva(it)
            elseif (dist .eq. 3) then
              ui1(it) = (2.d0*eva(i+1)/(3.d0*tva(it)) *
     &                  dexp(lngama(5.d0/(2.d0*dru_val))-
     &                  lngama(3.d0/(2.d0*dru_val))))
            endif
            uj1(it) = ui1(it) - ut(it)

C-----------------------------------------------------------------------
C  add in component part of the rate parameter in this interval
C-----------------------------------------------------------------------

            if     (ifint .eq. 1) then

              if     (dist .eq. 0) then
                sumi(1) = expi(it)*((tva(it)*ui(it))**(2d0)*w0 +
     &                       (2d0*tva(it)**(2d0)*w0+tva(it)*w1)*ui(it) +
     &                       (2d0*tva(it)**(2d0)*w0+tva(it)*w1+w2))
     &                  - expi1(it)*((tva(it)*ui1(it))**(2d0)*w0 +
     &                       (2d0*tva(it)**(2d0)*w0+tva(it)*w1)*ui1(it)+
     &                       (2d0*tva(it)**(2d0)*w0+tva(it)*w1+w2))
              elseif (dist .eq. 1) then
                ww0 = w0 * ((2.0d0*kappa-3.0d0)*tva(it)/2.d0)**(2.d0)
                ww1 = w1 * ((2.0d0*kappa-3.0d0)*tva(it)/2.d0)
                ww2 = w2
                sumi(1) = (1.d0/kappa*(1.d0+ui(it))**(-kappa) *
     &                        ( ww0*ui(it)**(2.d0)+ww1*ui(it)+ww2 ) +
     &                       1.d0/(kappa*(kappa-1.d0))*(1.d0+ui(it))**
     &                        (1.d0-kappa)* ( ww0*2.d0*ui(it)+ww1 ) +
     &                        1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                        (1.d0+ui(it))**(2.d0-kappa) * 2.d0*ww0)
     &                    - (1.d0/kappa*(1.d0+ui1(it))**(-kappa) *
     &                        ( ww0*ui1(it)**(2.d0)+ww1*ui1(it)+ww2 )+
     &                      1.d0/(kappa*(kappa-1.d0))*(1.d0+ui1(it))**
     &                        (1.d0-kappa)* ( ww0*2.d0*ui1(it)+ww1 ) +
     &                        1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                        (1.d0+ui1(it))**(2.d0-kappa)*2.d0*ww0)
                sumi(2) =  (1.d0/kappa*(1.d0+uj(it))**(-kappa) *
     &                       ( ww0*ui(it)**(2.d0)+ww1*ui(it)+ww2 ) +
     &                       1.d0/(kappa*(kappa-1.d0))*(1.d0+uj(it))**
     &                       (1.d0-kappa)* ( ww0*2.d0*ui(it)+ww1 ) +
     &                       1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                       (1.d0+uj(it))**(2.d0-kappa) * 2.d0*ww0)
     &                   - (1.d0/kappa*(1.d0+uj1(it))**(-kappa) *
     &                       ( ww0*ui1(it)**(2.d0)+ww1*ui1(it)+ww2 ) +
     &                      1.d0/(kappa*(kappa-1.d0))*(1.d0+uj1(it))**
     &                       (1.d0-kappa)* ( ww0*2.d0*ui1(it)+ww1 ) +
     &                       1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                       (1.d0+uj1(it))**(2.d0-kappa) * 2.d0*ww0)
              elseif (dist .eq. 3) then
                ww0 = w0 * (3.d0*tva(it) / (2.d0*
     &                dexp(lngama(5.d0/(2.d0*dru_val))-
     &                lngama(3.d0/(2.d0*dru_val)))))**2.d0
                ww1 = w1 * (3.d0*tva(it) / (2.d0*
     &                dexp(lngama(5.d0/(2.d0*dru_val))-
     &                lngama(3.d0/(2.d0*dru_val)))))
                ww2 = w2
                if (ui1(it)**dru_val.lt.1.d0/dru_val+1d0) then
                  igam1i = dexp(lngama(1.d0/dru_val))*
     &                     (ingama(1.d0/dru_val,ui1(it)**dru_val) -
     &                     ingama(1.d0/dru_val,ui(it)**dru_val))
                else
                  igam1i = dexp(lngama(1.d0/dru_val))*
     &                     (ingamq(1.d0/dru_val,ui(it)**dru_val) -
     &                     ingamq(1.d0/dru_val,ui1(it)**dru_val))
                endif
                if (ui1(it)**dru_val.lt.2.d0/dru_val+1d0) then
                  igam2i = dexp(lngama(2.d0/dru_val))*
     &                     (ingama(2.d0/dru_val,ui1(it)**dru_val) -
     &                     ingama(2.d0/dru_val,ui(it)**dru_val))
                else
                  igam2i = dexp(lngama(2.d0/dru_val))*
     &                     (ingamq(2.d0/dru_val,ui(it)**dru_val) -
     &                     ingamq(2.d0/dru_val,ui1(it)**dru_val))
                endif
                if (ui1(it)**dru_val.lt.3.d0/dru_val+1d0) then
                  igam3i = dexp(lngama(3.d0/dru_val))*
     &                     (ingama(3.d0/dru_val,ui1(it)**dru_val) -
     &                     ingama(3.d0/dru_val,ui(it)**dru_val))
                else
                  igam3i = dexp(lngama(3.d0/dru_val))*
     &                     (ingamq(3.d0/dru_val,ui(it)**dru_val) -
     &                     ingamq(3.d0/dru_val,ui1(it)**dru_val))
                endif
                if (uj1(it)**dru_val.lt.1.d0/dru_val+1d0) then
                  igam1j = dexp(lngama(1.d0/dru_val))*
     &                     (ingama(1.d0/dru_val,uj1(it)**dru_val) -
     &                     ingama(1.d0/dru_val,uj(it)**dru_val))
                else
                  igam1j = dexp(lngama(1.d0/dru_val))*
     &                     (ingamq(1.d0/dru_val,uj(it)**dru_val) -
     &                     ingamq(1.d0/dru_val,uj1(it)**dru_val))
                endif
                if (uj1(it)**dru_val.lt.2.d0/dru_val+1d0) then
                  igam2j = dexp(lngama(2.d0/dru_val))*
     &                     (ingama(2.d0/dru_val,uj1(it)**dru_val) -
     &                     ingama(2.d0/dru_val,uj(it)**dru_val))
                else
                  igam2j = dexp(lngama(2.d0/dru_val))*
     &                     (ingamq(2.d0/dru_val,uj(it)**dru_val) -
     &                     ingamq(2.d0/dru_val,uj1(it)**dru_val))
                endif
                if (uj1(it)**dru_val.lt.3.d0/dru_val+1d0) then
                  igam3j = dexp(lngama(3.d0/dru_val))*
     &                     (ingama(3.d0/dru_val,uj1(it)**dru_val) -
     &                     ingama(3.d0/dru_val,uj(it)**dru_val))
                else
                  igam3j = dexp(lngama(3.d0/dru_val))*
     &                     (ingamq(3.d0/dru_val,uj(it)**dru_val) -
     &                     ingamq(3.d0/dru_val,uj1(it)**dru_val))
                endif
                sumi(1) = ww0*igam3i + ww1*igam2i + ww2*igam1i
                sumi(2) = ww0*igam3j +
     &                    (2d0*ut(it)*ww0+ww1)*igam2j +
     &                    (ww0*ut(it)*ut(it)+ww1*ut(it)+ww2)*igam1j
              endif
            elseif (ifint .eq. 2) then
              if     (dist .eq. 0) then
                sumi(1) =   expi(it)*(w0*ee2(ui(it))
     &                        /(ui(it)*tva(it)**2.d0) +
     &                        w1*eei(ui(it))/tva(it)+w2)
     &                    - expi1(it)*(w0*ee2(ui1(it))
     &                        /(ui1(it)*tva(it)**2.d0) +
     &                        w1*eei(ui1(it))/tva(it)+w2)
              endif
            endif
            if (sumi(1).lt.0) sumi(1) = 0.d0
            if (sumi(2).lt.0) sumi(2) = 0.d0
            sum(1,it)=sum(1,it)+sumi(1)
            sum(2,it)=sum(2,it)+sumi(2)
   30     continue

C-----------------------------------------------------------------------
C  move reference from i to i+1 for next step of quadrature
C-----------------------------------------------------------------------

          do it=1,nt
            ui(it)=ui1(it)
            uj(it)=uj1(it)
            expi(it)=expi1(it)
          end do

   40   continue

C-----------------------------------------------------------------------
C  complete computation for last energy point
C-----------------------------------------------------------------------

        do 50 it = 1,nt

          if     (dist .eq. 0) then
            ui1(it) = eva(ne)/tva(it)
            expi1(it)=dexp(-ui1(it)+ut(it))
          elseif (dist .eq. 1) then
            ui1(it) = 2.d0*eva(ne)/((2.d0*kappa-3.d0)*tva(it))
          elseif (dist .eq. 3) then
            ui1(it) = (2.d0*eva(ne)/(3.d0*tva(it)) *
     &                dexp(lngama(5.d0/(2.d0*dru_val))-
     &                lngama(3.d0/(2.d0*dru_val))))
          endif
          uj1(it) = ui1(it) - ut(it)

C-----------------------------------------------------------------------
C  add in component part of the rate parameter in last interval
C-----------------------------------------------------------------------

          if     (ifint .eq. 1) then

            if     (dist .eq. 0) then
              sumn(1) =  expi(it)*((tva(it)*ui(it))**(2.d0)*w0 +
     &                     (2.d0*tva(it)**(2.d0)*w0+tva(it)*w1)*ui(it)+
     &                     (2.d0*tva(it)**(2.d0)*w0+tva(it)*w1+w2))
     &                 - expi1(it)*((tva(it)*ui1(it))**(2.d0)*w0 +
     &                     (2.d0*tva(it)**(2.d0)*w0+tva(it)*w1)*ui1(it)+
     &                     (2.d0*tva(it)**(2.d0)*w0+tva(it)*w1+w2))
            elseif (dist .eq. 1) then
              ww0 = w0 * ((2.0d0*kappa-3.0d0)*tva(it)/2.d0)**2.d0
              ww1 = w1 * ((2.0d0*kappa-3.0d0)*tva(it)/2.d0)
              ww2 = w2
              sumn(1) =   (1.d0/kappa*(1.d0+ui(it))**(-kappa) *
     &                      ( ww0*ui(it)**(2.d0)+ww1*ui(it)+ww2 ) +
     &                      1.d0/(kappa*(kappa-1.d0))*(1.d0+ui(it))**
     &                      (1.d0-kappa)* ( ww0*2.d0*ui(it)+ww1 ) +
     &                      1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                      (1.d0+ui(it))**(2.d0-kappa) * 2.d0*ww0)
     &                  - (1.d0/kappa*(1.d0+ui1(it))**(-kappa) *
     &                      ( ww0*ui1(it)**(2.d0)+ww1*ui1(it)+ww2 ) +
     &                      1.d0/(kappa*(kappa-1.d0))*(1.d0+ui1(it))**
     &                      (1.d0-kappa)* ( ww0*2.d0*ui1(it)+ww1 ) +
     &                      1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                      (1.d0+ui1(it))**(2.d0-kappa) * 2.d0*ww0)
              sumn(2) =  (1.d0/kappa*(1.d0+uj(it))**(-kappa) *
     &                     ( ww0*ui(it)**(2.d0)+ww1*ui(it)+ww2 ) +
     &                     1.d0/(kappa*(kappa-1.d0))*(1.d0+uj(it))**
     &                     (1.d0-kappa)* ( ww0*2.d0*ui(it)+ww1 ) +
     &                     1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                     (1.d0+uj(it))**(2.d0-kappa) * 2.d0*ww0)
     &                 - (1.d0/kappa*(1.d0+uj1(it))**(-kappa) *
     &                     ( ww0*ui1(it)**(2.d0)+ww1*ui1(it)+ww2 ) +
     &                     1.d0/(kappa*(kappa-1.d0))*(1.d0+uj1(it))**
     &                     (1.d0-kappa)* ( ww0*2.d0*ui1(it)+ww1 ) +
     &                     1.d0/(kappa*(kappa-1.d0)*(kappa-2.d0))*
     &                     (1.d0+uj1(it))**(2.d0-kappa) * 2.d0*ww0)
            elseif (dist .eq. 3) then
              ww0 = w0 * (3.d0*tva(it) / (2.d0*
     &              dexp(lngama(5.d0/(2.d0*dru_val))-
     &              lngama(3.d0/(2.d0*dru_val)))))**2.d0
              ww1 = w1 * (3.d0*tva(it) / (2.d0*
     &              dexp(lngama(5.d0/(2.d0*dru_val))-
     &              lngama(3.d0/(2.d0*dru_val)))))
              ww2 = w2
              if (ui1(it)**dru_val.lt.1.d0/dru_val+1d0) then
                igam1i = dexp(lngama(1.d0/dru_val))*
     &                   (ingama(1.d0/dru_val,ui1(it)**dru_val) -
     &                   ingama(1.d0/dru_val,ui(it)**dru_val))
              else
                igam1i = dexp(lngama(1.d0/dru_val))*
     &                   (ingamq(1.d0/dru_val,ui(it)**dru_val) -
     &                   ingamq(1.d0/dru_val,ui1(it)**dru_val))
              endif
              if (ui1(it)**dru_val.lt.2.d0/dru_val+1d0) then
                igam2i = dexp(lngama(2.d0/dru_val))*
     &                   (ingama(2.d0/dru_val,ui1(it)**dru_val) -
     &                   ingama(2.d0/dru_val,ui(it)**dru_val))
              else
                igam2i = dexp(lngama(2.d0/dru_val))*
     &                   (ingamq(2.d0/dru_val,ui(it)**dru_val) -
     &                   ingamq(2.d0/dru_val,ui1(it)**dru_val))
              endif
              if (ui1(it)**dru_val.lt.3.d0/dru_val+1d0) then
                igam3i = dexp(lngama(3.d0/dru_val))*
     &                   (ingama(3.d0/dru_val,ui1(it)**dru_val) -
     &                   ingama(3.d0/dru_val,ui(it)**dru_val))
              else
                igam3i = dexp(lngama(3.d0/dru_val))*
     &                   (ingamq(3.d0/dru_val,ui(it)**dru_val) -
     &                   ingamq(3.d0/dru_val,ui1(it)**dru_val))
              endif
              if (uj1(it)**dru_val.lt.1.d0/dru_val+1d0) then
                igam1j = dexp(lngama(1.d0/dru_val))*
     &                   (ingama(1.d0/dru_val,uj1(it)**dru_val) -
     &                   ingama(1.d0/dru_val,uj(it)**dru_val))
              else
                igam1j = dexp(lngama(1.d0/dru_val))*
     &                   (ingamq(1.d0/dru_val,uj(it)**dru_val) -
     &                   ingamq(1.d0/dru_val,uj1(it)**dru_val))
              endif
              if (uj1(it)**dru_val.lt.2.d0/dru_val+1d0) then
                igam2j = dexp(lngama(2.d0/dru_val))*
     &                   (ingama(2.d0/dru_val,uj1(it)**dru_val) -
     &                   ingama(2.d0/dru_val,uj(it)**dru_val))
              else
                igam2j = dexp(lngama(2.d0/dru_val))*
     &                   (ingamq(2.d0/dru_val,uj(it)**dru_val) -
     &                   ingamq(2.d0/dru_val,uj1(it)**dru_val))
              endif
              if (uj1(it)**dru_val.lt.3.d0/dru_val+1d0) then
                igam3i = dexp(lngama(3.d0/dru_val))*
     &                   (ingama(3.d0/dru_val,uj1(it)**dru_val) -
     &                   ingama(3.d0/dru_val,uj(it)**dru_val))
              else
                igam3i = dexp(lngama(3.d0/dru_val))*
     &                   (ingamq(3.d0/dru_val,uj(it)**dru_val) -
     &                   ingamq(3.d0/dru_val,uj1(it)**dru_val))
              endif
              sumn(1) = ww0*igam3i + ww1*igam2i + ww2*igam1i
              sumn(2) = ww0*igam3j +
     &                  (2d0*ut(it)*ww0+ww1)*igam2j +
     &                  (ww0*ut(it)*ut(it)+ww1*ut(it)+ww2)*igam1j
            endif
          elseif (ifint.eq.2.and.dist.eq.0) then
             sumn(1) = expi(it)*(w0*ee2(ui(it))/(ui(it)*tva(it)**2)
     &                    + w1*eei(ui(it))/tva(it)+w2)
     &               - expi1(it)*(w0*ee2(ui1(it))/(ui1(it)*tva(it)**2)
     &                    + w1*eei(ui1(it))/tva(it)+w2)
          endif
          if (sumn(1).lt.0) sumn(1) = 0.d0
          if (sumn(2).lt.0) sumn(2) = 0.d0
          sum(1,it)=sum(1,it)+sumn(1)
          sum(2,it)=sum(2,it)+sumn(2)

   50   continue
      endif

C-----------------------------------------------------------------------
C  Add in part of the rate parameter from last point to infinity.
C  Asymptotic behaviour determined by itype.
C-----------------------------------------------------------------------

      do 90 it = 1 , nt

C-----------------------------------------------------------------------
C  Maxwellian distribution
C-----------------------------------------------------------------------
        
        if (dist .eq. 0) then

          if (itype .eq. 1) then
C-----------------------------------------------------------------------
C  if infinite energy limit point exists then draw straight line in
C  Burgess-Tully space from last point to limit point, with C = e
C-----------------------------------------------------------------------
            if (lbeth) then
              c2      = -oma(ne)+bthp*dlog(ui1(it)/ut(it)-1d0+dexp(1d0))
              sumu(1) = expi1(it)*(bthp*(dlog(ui1(it)/ut(it)+
     &                  dexp(1d0))+eei(ui1(it)+dexp(1d0)*ut(it)))-c2)
            else
              if (oma(ne).le.oma(ne-1)) then
                if ((it.eq.1).and.(.not.lamb1)) then
                  write(i4unit(-1),1011)'electron excitation transition'
     &              ,itrn
                  write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                  lamb1=.true.
                endif
                sumu(1) = oma(ne)*expi1(it)*(1.d0+eei(ui1(it))/
     &                    dlog(xa(ne)))
                nct1 = nct1+1
              else
                b = 0.0D0
                do iter=1,20
                  b = -xa(ne-1)+(xa(ne)+b)**(oma(ne-1)/oma(ne))
                enddo
                a = oma(ne)/(dlog(xa(ne)+b))
                sumu(1) = expi1(it)*(oma(ne)+
     &                    a*eei(ui1(it)+evt*b/tva(it)))
              endif
            endif

          elseif (itype .eq. 2) then
            if (lbeth) then
              sumu(1) = expi1(it)*(bthp+ui1(it)*eei(ui1(it))*
     &                  (oma(ne)-bthp))
            else
              a = oma(ne)-(oma(ne)-oma(ne-1))/(1.d0-xa(ne)/xa(ne-1))
              b = (oma(ne)-oma(ne-1))/(1.d0/xa(ne)-1.d0/xa(ne-1))
              if (a.lt.0.d0.or.xa(ne).lt.-b/a) then
               if ((it.eq.1).and.(.not.lamb2)) then
                  write(i4unit(-1),1011)'electron excitation transition'
     &              ,itrn
                  write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                  lamb2 = .true.
                endif
                sumu(1) = oma(ne)*expi1(it)
		nct2 = nct2+1
              else
                sumu(1) = expi1(it)*(a+evt/tva(it)*b*eei(ui1(it)))
              endif
            endif

          elseif (itype .eq. 3) then
            if (oma(ne).ge.oma(ne-1)) then
              if ((it.eq.1).and.(.not.lamb3)) then
                write(i4unit(-1),1011)'electron excitation transition'
     &            ,itrn
                write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                lamb3 = .true.
              endif
              b = 100.d0
	      nct3 = nct3+1
            else
              b = (dsqrt(oma(ne))*xa(ne)-dsqrt(oma(ne-1))*
     &            xa(ne-1))/(dsqrt(oma(ne-1))-dsqrt(oma(ne)))
            endif
            a = oma(ne)*(xa(ne)+b)**2.d0
            sumu(1) = a*evt/tva(it)*evt/tva(it)*expi1(it)*
     &                ee2(ui1(it)+evt/tva(it)*b)

          else
            sumu(1) = 0.d0
          endif

C-----------------------------------------------------------------------
C  kappa distribution
C-----------------------------------------------------------------------
        elseif (dist .eq. 1) then

          if     (itype .eq. 1) then
            if (lbeth) then
              a = -oma(ne)+bthp*dlog(xa(ne)-1d0+dexp(1d0))
            else
              b = 0.d0
              if (oma(ne).le.oma(ne-1)) then
                if ((it.eq.1).and.(.not.lamb1)) then
                  write(i4unit(-1),1011)'electron excitation transition'
     &              ,itrn
                  write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                  lamb1=.true.
                endif
                nct1 = nct1+1
              else
                do iter=1,20
                  b = -xa(ne-1)+(xa(ne)+b)**(oma(ne-1)/oma(ne))
                enddo
              endif
              a = oma(ne)/(dlog(xa(ne)+b))
            endif
            estep = 2.0d0*(eva(ne)-eva(ne-1))/
     &              ((2.0d0*kappa-3.0d0)*tva(it))/10.d0
            uiu1 = ui1(it)
            uju1 = uj1(it)
            intgl_tot1 = 0.d0
            intgl_tot2 = 0.d0
            do ie = 1,max
              uiu  = uiu1
              uiu1 = uiu1+estep
              uju  = uju1
              uju1 = uju1+estep
              if (lbeth) then
                intgl1 = ((1.d0+uiu)**(-kappa-1.d0)*
     &                    (bthp*dlog(uiu*tva(it)*(kappa-1.5d0)/evt-
     &                     1d0+dexp(1d0))-a) +
     &                   (1.d0+uiu1)**(-kappa-1.d0)*
     &                   (bthp*dlog(uiu1*tva(it)*(kappa-1.5d0)/evt-
     &                     1d0+dexp(1d0))-a))*
     &                   estep/2.d0
                intgl2 = ((1.d0+uju)**(-kappa-1.d0)*
     &                    (bthp*dlog(uiu*tva(it)*(kappa-1.5d0)/evt-
     &                     1d0+dexp(1d0))-a) +
     &                   (1.d0+uju1)**(-kappa-1.d0)*
     &                   (bthp*dlog(uiu1*tva(it)*(kappa-1.5d0)/evt-
     &                     1d0+dexp(1d0))-a))*
     &                   estep/2.d0
              else
                intgl1 = a*((1.d0+uiu)**(-kappa-1.d0)*
     &                   dlog(uiu*tva(it)*(kappa-1.5d0)/evt+b) +
     &                   (1.d0+uiu1)**(-kappa-1.d0)*
     &                   dlog(uiu1*tva(it)*(kappa-1.5d0)/evt+b)) *
     &                   estep/2.d0
                intgl2 = a*((1.d0+uju)**(-kappa-1.d0)*
     &                   dlog(uiu*tva(it)*(kappa-1.5d0)/evt+b) +
     &                   (1.d0+uju1)**(-kappa-1.d0)*
     &                   dlog(uiu1*tva(it)*(kappa-1.5d0)/evt+b)) *
     &                   estep/2.d0
              endif
              intgl_tot1 = intgl_tot1 + intgl1
              intgl_tot2 = intgl_tot2 + intgl2
              if (intgl_tot1.gt.0.0) then
                 if (intgl1/intgl_tot1.lt.tol) goto 310
              endif
              estep=estep*1.1d0
            enddo
  310       continue
            sumu(1) = intgl_tot1
            sumu(2) = intgl_tot2

          elseif (itype .eq. 2) then
            if (lbeth) then
              b = -(xa(ne)-1d0+dexp(1d0))/dexp(1d0)*(oma(ne)-bthp)
              a = oma(ne)-b
            else
              a = oma(ne)-(oma(ne)-oma(ne-1))/(1.d0-xa(ne)/xa(ne-1))
              b = (oma(ne)-oma(ne-1))/(1.d0/xa(ne)-1.d0/xa(ne-1))
              if (a.lt.0.d0.or.xa(ne).lt.-b/a) then
                if ((it.eq.1).and.(.not.lamb2)) then
                  write(i4unit(-1),1011)'electron excitation transition'
     &              ,itrn
                  write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                  lamb2 = .true.
                endif
                sumu(1) = oma(ne)/kappa*(1.0d0+ui1(it))**(-kappa)*
     &                    dexp(evt/tva(it))
                sumu(2) = oma(ne)/kappa*(1.0d0+ui1(it)-ut(it))**(-kappa)
		nct2 = nct2+1
                goto 810
              endif
            endif
            estep = 2.0d0*(eva(ne)-eva(ne-1))/
     &              ((2.0d0*kappa-3.0d0)*tva(it))/10.d0
            uiu1 = ui1(it)
            uju1 = uj1(it)
            intgl_tot1 = 0.d0
            intgl_tot2 = 0.d0
            do ie = 1,max
              uiu  = uiu1
              uiu1 = uiu1+estep
              uju  = uju1
              uju1 = uju1+estep
              if (lbeth) then
                intgl1 = ((1.d0+uiu)**(-kappa-1.d0)*
     &                   (uiu*tva(it)*(kappa-1.5d0)/evt-1d0)/
     &                   (uiu*tva(it)*(kappa-1.5d0)/evt-1d0+dexp(1d0)) +
     &                   (1.d0+uiu1)**(-kappa-1.d0)*
     &                   (uiu1*tva(it)*(kappa-1.5d0)/evt-1d0)/
     &                   (uiu1*tva(it)*(kappa-1.5d0)/evt-1d0+dexp(1d0)))
     &                   *estep/2.d0
                intgl2 = ((1.d0+uju)**(-kappa-1.d0)*
     &                   (uiu*tva(it)*(kappa-1.5d0)/evt-1d0)/
     &                   (uiu*tva(it)*(kappa-1.5d0)/evt-1d0+dexp(1d0)) +
     &                   (1.d0+uju1)**(-kappa-1.d0)*
     &                   (uiu1*tva(it)*(kappa-1.5d0)/evt-1d0)/
     &                   (uiu1*tva(it)*(kappa-1.5d0)/evt-1d0+dexp(1d0)))
     &                   *estep/2.d0
              else
                intgl1 = ((1.d0+uiu)**(-kappa-1.d0)/uiu +
     &                   (1.d0+uiu1)**(-kappa-1.d0)/uiu1) *
     &                   (evt/(tva(it)*(kappa-1.5d0)))*estep/2.d0
                intgl2 = ((1.d0+uju)**(-kappa-1.d0)/uiu +
     &                   (1.d0+uju1)**(-kappa-1.d0)/uiu1) *
     &                   (evt/(tva(it)*(kappa-1.5d0)))*estep/2.d0
              endif
              intgl_tot1 = intgl_tot1 + intgl1
              intgl_tot2 = intgl_tot2 + intgl2
              if (intgl_tot1.gt.0.0) then
                 if (intgl1/intgl_tot1.lt.tol) goto 300
              endif
              estep=estep*1.1d0
            enddo
  300       continue
            sumu(1) = a/kappa*(1.0d0+ui1(it))**(-kappa)+b*intgl_tot1
            sumu(2) = a/kappa*(1.0d0+uj1(it))**(-kappa)+b*intgl_tot2
  810       continue

          elseif (itype .eq. 3) then
            if (oma(ne).ge.oma(ne-1)) then
              if ((it.eq.1).and.(.not.lamb3)) then
                write(i4unit(-1),1011)'electron excitation transition'
     &            ,itrn
                write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                lamb3 = .true.
              endif
              b = 100.d0
	      nct3 = nct3+1
            else
              b = (dsqrt(oma(ne))*xa(ne)-dsqrt(oma(ne-1))
     &            *xa(ne-1))/(dsqrt(oma(ne-1))-dsqrt(oma(ne)))
            endif
            a = oma(ne)*(xa(ne)+b)**2.d0
            estep = 2.0d0*(eva(ne)-eva(ne-1))/
     &              ((2.0d0*kappa-3.0d0)*tva(it))/10.d0
            uiu1 = ui1(it)
            uju1 = uj1(it)
            intgl_tot1 = 0.d0
            intgl_tot2 = 0.d0
            do ie = 1,max
              uiu  = uiu1
              uiu1 = uiu1+estep
              uju  = uju1
              uju1 = uju1+estep
              intgl1 = ((1.d0+uiu)**(-kappa-1.d0)/
     &                 (uiu*tva(it)*(kappa-1.5d0)/evt+b)**2d0 +
     &                 (1.d0+uiu1)**(-kappa-1.d0)/
     &                 (uiu1*tva(it)*(kappa-1.5d0)/evt+b)**2d0) *
     &                 estep/2.d0
              intgl2 = ((1.d0+uju)**(-kappa-1.d0)/
     &                 (uiu*tva(it)*(kappa-1.5d0)/evt+b)**2d0 +
     &                 (1.d0+uju1)**(-kappa-1.d0)/
     &                 (uiu1*tva(it)*(kappa-1.5d0)/evt+b)**2d0) *
     &                 estep/2.d0
              intgl_tot1 = intgl_tot1 + intgl1
              intgl_tot2 = intgl_tot2 + intgl2
              if (intgl_tot1.gt.0.0) then
                 if (intgl1/intgl_tot1.lt.tol) goto 320
              endif
              estep=estep*1.1d0
            enddo
  320       continue
            sumu(1) = a*intgl_tot1
            sumu(2) = a*intgl_tot2
          else
            sumu(1) = 0.d0
            sumu(2) = 0.d0
          endif

C-----------------------------------------------------------------------
C  Druyvesteyn distribution
C-----------------------------------------------------------------------
        elseif (dist .eq. 3) then

          if     (itype .eq. 1) then
            if (lbeth) then
              a = -oma(ne)+bthp*dlog(xa(ne)-1d0+dexp(1d0))
            else
              b = 0.d0
              if (oma(ne).le.oma(ne-1)) then
                if ((it.eq.1).and.(.not.lamb1)) then
                  write(i4unit(-1),1011)'electron excitation transition'
     &              ,itrn
                  write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                  lamb1=.true.
                endif
                nct1 = nct1+1
              else
                do iter=1,20
                  b = -xa(ne-1)+(xa(ne)+b)**(oma(ne-1)/oma(ne))
                enddo
              endif
              a = oma(ne)/(dlog(xa(ne)+b))
            endif
            estep = (2.d0*(eva(ne)-eva(ne-1))/(3.d0*tva(it)) *
     &              dexp(lngama(5.d0/(2.d0*dru_val))-
     &              lngama(3.d0/(2.d0*dru_val))))/10d0
            uiu1 = ui1(it)
            uju1 = uj1(it)
            intgl_tot1 = 0.d0
            intgl_tot2 = 0.d0
            do ie = 1,max
              uiu  = uiu1
              uiu1 = uiu1+estep
              uju  = uju1
              uju1 = uju1+estep
              if (lbeth) then
                intgl1 = (dexp(-uiu**dru_val)*(bthp*dlog(1.5d0*tva(it)*
     &                    uiu*dexp(lngama(3d0/(2d0*dru_val))-
     &                    lngama(5d0/(2d0*dru_val)))/evt-1d0+dexp(1d0))
     &                    -a)
     &                 + dexp(-uiu1**dru_val)*(bthp*dlog(1.5d0*tva(it)*
     &                    uiu1*dexp(lngama(3d0/(2d0*dru_val))-
     &                    lngama(5d0/(2d0*dru_val)))/evt-1d0+dexp(1d0))
     &                    -a))
     &                  *estep/2.d0
                intgl2 = (dexp(-uju**dru_val)*(bthp*dlog(1.5d0*tva(it)*
     &                    uiu*dexp(lngama(3d0/(2d0*dru_val))-
     &                    lngama(5d0/(2d0*dru_val)))/evt-1d0+dexp(1d0))
     &                    -a)
     &                 + dexp(-uju1**dru_val)*(bthp*dlog(1.5d0*tva(it)*
     &                    uiu1*dexp(lngama(3d0/(2d0*dru_val))-
     &                    lngama(5d0/(2d0*dru_val)))/evt-1d0+dexp(1d0))
     &                    -a))
     &                  *estep/2.d0
              else
                intgl1 = a*( dexp(-uiu**dru_val)*dlog(1.5d0*tva(it)*uiu*
     &                       dexp(lngama(3d0/(2d0*dru_val))-
     &                       lngama(5d0/(2d0*dru_val)))/evt+b) +
     &                     dexp(-uiu1**dru_val)*dlog(1.5d0*tva(it)*uiu1*
     &                       dexp(lngama(3d0/(2d0*dru_val))-
     &                       lngama(5d0/(2d0*dru_val)))/evt+b) ) *
     &                   estep/2.d0
                intgl2 = a*( dexp(-uju**dru_val)*dlog(1.5d0*tva(it)*uiu*
     &                       dexp(lngama(3d0/(2d0*dru_val))-
     &                       lngama(5d0/(2d0*dru_val)))/evt+b) +
     &                     dexp(-uju1**dru_val)*dlog(1.5d0*tva(it)*uiu1*
     &                       dexp(lngama(3d0/(2d0*dru_val))-
     &                       lngama(5d0/(2d0*dru_val)))/evt+b) ) *
     &                   estep/2.d0
              endif
              intgl_tot1 = intgl_tot1 + intgl1
              intgl_tot2 = intgl_tot2 + intgl2
              if (intgl_tot1.GT.0.0) then
                 if (intgl1/intgl_tot1.lt.tol) goto 330
              endif
              estep=estep*1.1d0
            enddo
  330       continue
            sumu(1) = dru_val*intgl_tot1
            sumu(2) = dru_val*intgl_tot2

          elseif (itype .eq. 2) then
            if (lbeth) then
              b = -(xa(ne)-1d0+dexp(1d0))/dexp(1d0)*(oma(ne)-bthp)
              a = oma(ne)-b
            else
              a = oma(ne)-(oma(ne)-oma(ne-1))/(1.d0-xa(ne)/xa(ne-1))
              b = (oma(ne)-oma(ne-1))/(1.d0/xa(ne)-1.d0/xa(ne-1))
              if (a.lt.0.d0.or.xa(ne).lt.-b/a) then
                if ((it.eq.1).and.(.not.lamb2)) then
                  write(i4unit(-1),1011)'electron excitation transition'
     &              ,itrn
                  write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                  lamb2 = .true.
                endif
                sumu(1) = oma(ne)*
     &                    ingamq(1.d0/dru_val,ui1(it)**dru_val)
     &                    *dexp(lngama(1.d0/dru_val)+evt/tva(it))
                sumu(2) = oma(ne)*
     &                    ingamq(1.d0/dru_val,(ui1(it)-ut(it))**dru_val)
     &                    *dexp(lngama(1.d0/dru_val))
		nct2 = nct2+1
                goto 820
              endif
            endif
            estep = (2.d0*(eva(ne)-eva(ne-1))/(3.d0*tva(it)) *
     &              dexp(lngama(5.d0/(2.d0*dru_val))-
     &              lngama(3.d0/(2.d0*dru_val))))/10d0
            uiu1 = ui1(it)
            uju1 = uj1(it)
            intgl_tot1 = 0.d0
            intgl_tot2 = 0.d0
            do ie = 1,max
              uiu  = uiu1
              uiu1 = uiu1+estep
              uju  = uju1
              uju1 = uju1+estep
              if (lbeth) then
                intgl1 = (dexp(-uiu**dru_val)*
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0)/
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0+dexp(1d0)) +
     &                   dexp(-uiu1**dru_val)*
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0)/
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0+dexp(1d0)))
     &                   *estep/2.d0
                intgl2 =(dexp(-uju**dru_val)*
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0)/
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0+dexp(1d0)) +
     &                   dexp(-uju1**dru_val)*
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0)/
     &                     (uiu*(1.5d0*tva(it))/evt *
     &                     dexp(lngama(3.d0/(2.d0*dru_val))-
     &                     lngama(5.d0/(2.d0*dru_val)))-1d0+dexp(1d0)))
     &                   *estep/2.d0
              else
                intgl1 = (dexp(-uiu**dru_val)/uiu +
     &                   dexp(-uiu1**dru_val)/uiu1) *
     &                   (evt/(1.5d0*tva(it)) *
     &                   dexp(lngama(5.d0/(2.d0*dru_val))-
     &                   lngama(3.d0/(2.d0*dru_val))))*estep/2.d0
                intgl2 = (dexp(-uju**dru_val)/uiu +
     &                   dexp(-uju1**dru_val)/uiu1) *
     &                   (evt/(1.5d0*tva(it)) *
     &                   dexp(lngama(5.d0/(2.d0*dru_val))-
     &                   lngama(3.d0/(2.d0*dru_val))))*estep/2.d0
              endif
              intgl_tot1 = intgl_tot1 + intgl1
              intgl_tot2 = intgl_tot2 + intgl2
              if (intgl_tot1.GT.0.0) then
                if (intgl1/intgl_tot1.lt.tol) goto 340
              endif
              estep=estep*1.1d0
            enddo
  340       continue
            sumu(1) = (a*ingamq(1.d0/dru_val,ui1(it)**dru_val)*
     &                dexp(lngama(1.d0/dru_val)) + b*intgl_tot1)
            sumu(2) = a*ingamq(1.d0/dru_val,
     &                (ui1(it)-ut(it))**dru_val)
     &                *dexp(lngama(1.d0/dru_val)) + b*intgl_tot2
  820       continue

          elseif (itype .eq. 3) then
            if (oma(ne).ge.oma(ne-1)) then
              if ((it.eq.1).and.(.not.lamb3)) then
                write(i4unit(-1),1011)'electron excitation transition'
     &            ,itrn
                write(i4unit(-1),1013)
     &              'some x-sects ambiguous at high energy for type ',
     &              itype
                lamb3 = .true.
              endif
              b = 100.d0
	      nct3 = nct3+1
            else
              b = (dsqrt(oma(ne))*xa(ne)-dsqrt(oma(ne-1))
     &            *xa(ne-1))/(dsqrt(oma(ne-1))-dsqrt(oma(ne)))
            endif
            a = oma(ne)*(xa(ne)+b)**2.d0
            estep = (2.d0*(eva(ne)-eva(ne-1))/(3.d0*tva(it)) *
     &              dexp(lngama(5.d0/(2.d0*dru_val))-
     &              lngama(3.d0/(2.d0*dru_val))))/10d0
            uiu1 = ui1(it)
            uju1 = uj1(it)
            intgl_tot1 = 0.d0
            intgl_tot2 = 0.d0
            do ie = 1,max
              uiu  = uiu1
              uiu1 = uiu1+estep
              uju  = uju1
              uju1 = uju1+estep
              intgl1 = (dexp(-uiu**dru_val)/
     &                   (uiu*tva(it)*1.5d0*
     &                   dexp(lngama(3.d0/(2.d0*dru_val))-
     &                   lngama(5.d0/(2.d0*dru_val)))/evt+b)**2d0 +
     &                 dexp(-uiu1**dru_val)/
     &                   (uiu1*tva(it)*1.5d0*
     &                   dexp(lngama(3.d0/(2.d0*dru_val))-
     &                   lngama(5.d0/(2.d0*dru_val)))/evt+b)**2d0) *
     &                 estep/2.d0
              intgl2 = (dexp(-uju**dru_val)/
     &                   (uiu*tva(it)*1.5d0*
     &                   dexp(lngama(3.d0/(2.d0*dru_val))-
     &                   lngama(5.d0/(2.d0*dru_val)))/evt+b)**2d0 +
     &                 dexp(-uju1**dru_val)/
     &                   (uiu1*tva(it)*1.5d0*
     &                   dexp(lngama(3.d0/(2.d0*dru_val))-
     &                   lngama(5.d0/(2.d0*dru_val)))/evt+b)**2d0) *
     &                 estep/2.d0
              intgl_tot1 = intgl_tot1 + intgl1
              intgl_tot2 = intgl_tot2 + intgl2
              if (intgl_tot1.gt.0.0) then
                 if (intgl1/intgl_tot1.lt.tol) goto 350
              endif
              estep=estep*1.1d0
            enddo
  350       continue
            sumu(1) = a*intgl_tot1
            sumu(2) = a*intgl_tot2
          else
            sumu(1) = 0.d0
            sumu(2) = 0.d0
          endif
        endif
        if (sumu(1).lt.0d0) sumu(1) = 0.d0
        if (sumu(2).lt.0d0) sumu(2) = 0.d0
	
c	write(i4unit(-1),'(a,i5,1p2e9.2)')
c     &              'it,sumu(1),sumu(2)=',it,sumu(1),sumu(2)
         sum(1,it)=sum(1,it)+sumu(1)
         sum(2,it)=sum(2,it)+sumu(2)

C-----------------------------------------------------------------------
C  Add in part of the rate parameter from threshold to first point.
C  Threshold behaviour determined by itypt
C-----------------------------------------------------------------------

        if     (dist .eq. 0) then
          if     (itypt .eq. 1) then
            suml(1) = oma(1)*(1.d0-exp1(it))
          elseif (itypt .eq. 2) then
            suml(1) = 0.d0
          else
            suml(1) = 0.d0
          endif
        elseif (dist .eq. 1) then
          if     (itypt .eq. 1) then
            suml(1) = oma(1)/kappa*((1d0+ut(it))**(-kappa)-
     &                              (1d0+u1(it))**(-kappa))
            suml(2) = oma(1)/kappa*(1d0-(1d0+u1(it)-ut(it))**(-kappa))
          elseif (itypt .eq. 2) then
            suml(1) = 0.d0
            suml(2) = 0.d0
          else
            suml(1) = 0.d0
            suml(2) = 0.d0
          endif
        elseif (dist .eq. 3) then
          if     (itypt .eq. 1) then
            if (u1(it)**dru_val.lt.1.d0/dru_val+1d0) then
              suml(1) = oma(1)*dexp(lngama(1.d0/dru_val))*
     &                (ingama(1.d0/dru_val,u1(it)**dru_val)-
     &                ingama(1.d0/dru_val,ut(it)**dru_val))
            else
              suml(1) = oma(1)*dexp(lngama(1.d0/dru_val))*
     &                (ingamq(1.d0/dru_val,ut(it)**dru_val)-
     &                ingamq(1.d0/dru_val,u1(it)**dru_val))
            endif
            suml(2) = oma(1)*ingama(1.d0/dru_val,(u1(it)-ut(it))
     &                **dru_val)*dexp(lngama(1.d0/dru_val))
          elseif (itypt .eq. 2) then
            suml(1) = 0.d0
            suml(2) = 0.d0
          else
            suml(1) = 0.d0
            suml(2) = 0.d0
          endif
        endif
        if (suml(1).lt.0) suml(1) = 0.d0
        if (suml(2).lt.0) suml(2) = 0.d0
        sum(1,it)=sum(1,it)+suml(1)
        sum(2,it)=sum(2,it)+suml(2)

   90 continue

C-----------------------------------------------------------------------
C  numerical distribution
C-----------------------------------------------------------------------

      else

C-----------------------------------------------------------------------
C  multiply distribution function by energy^(-1/2) for interpolation
C-----------------------------------------------------------------------

        do it = 1,ntf
          if (en(it,1).lt.evt.and.en(it,nef).gt.evt) then
            eni(it,1) = evt
            fj(it,1) = f(it,1)/dsqrt(en(it,1))
            do i = 1,nef-1
              eni(it,i+1) = en(it,i+1)-en(it,1)+evt
              fj(it,i+1) = f(it,i+1)/dsqrt(en(it,i+1))

C-----------------------------------------------------------------------
C  determine distribution function at shifted energies
C-----------------------------------------------------------------------

              v1 = 1.d0/(en(it,i+1)-en(it,i))
              a0 = v1*dlog(fj(it,i)/fj(it,i+1))
              a1 = fj(it,i)*(fj(it,i)/fj(it,i+1))**(en(it,i)*v1)
              do j = 1,nef
                if (eni(it,j).ge.en(it,i).and.eni(it,j).lt.en(it,i+1))
     &            then
                  fi(it,j) = a1*dexp(-a0*(en(it,j)-en(it,1)+evt))
                elseif (eni(it,j).gt.en(it,nef)) then
                  if (nform2.eq.1) then
                    fi(it,j) = 0d0
                  elseif (nform2.eq.2) then
                    fi(it,j) = f(it,nef)/dsqrt(eni(it,j))*
     &                          (en(it,nef)/eni(it,j))**param2(1)
                  elseif (nform2.eq.3) then
                    fi(it,j) = f(it,nef)*dexp(param2(1)/tvf(it)*
     &                          (en(it,nef)-eni(it,j)))/dsqrt(eni(it,j))
                  endif
                endif
                if (fi(it,j).lt.0d0) fi(it,j)=0.d0
              enddo
            enddo
          else
            do i = 1,nef
              eni(it,i) = en(it,i)
              fi(it,i) = f(it,i)/dsqrt(en(it,i))
              fj(it,i) = fi(it,i)
            enddo
          endif
        enddo

C-----------------------------------------------------------------------
C  set up interpolation parameters for asymptotic limits of omega,
C  dependent on transition type
C-----------------------------------------------------------------------

        if     (itype.eq.1) then
          if (lbeth) then
            a = -oma(ne)+bthp*dlog(xa(ne)-1d0+dexp(1d0))
          else
            if (oma(ne).le.oma(ne-1)) then
              lhighe=.true.
              if ((it.eq.1).and.(.not.lamb1)) then		        
                write(i4unit(-1),1011)'electron excitation transition'
     &            ,itrn
                write(i4unit(-1),1013)
     &            'some x-sects ambiguous at high energy for type ',    
     &            itype 					        
                lamb1=.true.					        
              endif
              a = (oma(ne)-oma(ne-1))/(dlog(xa(ne))-dlog(xa(ne-1)))
              b = oma(ne)-a*dlog(xa(ne))
              nct1 = nct1+1
            else
              b = 0.0d0
              do iter=1,20
                b = -xa(ne-1)+(xa(ne)+b)**(oma(ne-1)/oma(ne))
              enddo
              a = oma(ne)/(dlog(xa(ne)+b))
            endif
          endif
        elseif (itype.eq.2) then
          if (lbeth) then
            b = -(xa(ne)-1d0+dexp(1d0))/dexp(1d0)*(oma(ne)-bthp)
            a = oma(ne)-b
          else
			atest1 = oma(ne)-(oma(ne)-oma(ne-1))
			atest2 = (1.d0-xa(ne)/xa(ne-1))
			a=atest1/atest2
C           a = oma(ne)-(oma(ne)-oma(ne-1))/(1.d0-xa(ne)/xa(ne-1))
            b = (oma(ne)-oma(ne-1))/(1.d0/xa(ne)-1.d0/xa(ne-1))
            if (a.lt.0.d0.or.xa(ne).lt.-b/a) then
              lhighe=.true.
              if ((it.eq.1).and.(.not.lamb2)) then		        
                write(i4unit(-1),1011)'electron excitation transition'
     &            ,itrn
                write(i4unit(-1),1013)
     &            'some x-sects ambiguous at high energy for type ',    
     &            itype 					        
                lamb2 = .true.					        
              endif
              nct2 = nct2+1
            endif
          endif

        elseif (itype.eq.3) then
          if (oma(ne).ge.oma(ne-1)) then
            if ((it.eq.1).and.(.not.lamb3)) then		      
              write(i4unit(-1),1011)'electron excitation transition'
     &          ,itrn
              write(i4unit(-1),1013)
     &      	  'some x-sects ambiguous at high energy for type ',  
     &      	  itype 					      
              lamb3 = .true.					      
            endif
            b = 100.d0
	    nct3 = nct3+1
          else
            b = (dsqrt(oma(ne))*xa(ne)-dsqrt(oma(ne-1))*
     &          xa(ne-1))/(dsqrt(oma(ne-1))-dsqrt(oma(ne)))
          endif
          a = oma(ne)*(xa(ne)+b)**2.d0

        endif

C-----------------------------------------------------------------------
C  convert collision strength energy grid to distribution
C  function energy grid
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C  linear interpolation
C-----------------------------------------------------------------------

        if     (ilinr .eq. 1) then

          do i = 1,ne-1
            c1 = -xa(i+1)*oma(i)/(xa(i)-xa(i+1))
     &           -xa(i)*oma(i+1)/(xa(i+1)-xa(i))
            c0 = oma(i)/(xa(i)-xa(i+1)) + oma(i+1)/(xa(i+1)-xa(i))
            do it = 1,ntf
              do j = 1,nef
                xf = eni(it,j)/evt
                if ((xf.ge.1.d0).and.(((xf.lt.xa(1)).and.(i.eq.1))
     &             .or.((xf.ge.xa(i)).and.(xf.lt.xa(i+1)))))
     &             omegai(it,j) = c0*xf + c1
                if (xf.ge.xa(ne)) then
                  if (itype.eq.1) then
                    if (lbeth) then
                      omegai(it,j) = bthp*dlog(xf-1d0+dexp(1d0))-a
                    else
                      if (lhighe) then
                        omegai(it,j) = a*dlog(xf)+b
                      else
                        omegai(it,j) = a*dlog(xf+b)
                      endif
                    endif
                  elseif (itype.eq.2) then
                    if (lbeth) then
                      omegai(it,j) = b*(xf-1d0)/(xf-1d0+dexp(1d0))+a
                    else
                      if (lhighe) then
                        omegai(it,j) = oma(ne)
                      else
                        omegai(it,j) = a+b/xf
                      endif
                    endif
                  elseif (itype.eq.3) then
                    omegai(it,j) = a/(xf+b)**2d0
                  endif
                endif
                xf = en(it,j)/evt+1d0
                if ((xf.ge.1.d0).and.(((xf.lt.xa(1)).and.(i.eq.1))
     &             .or.((xf.ge.xa(i)).and.(xf.lt.xa(i+1)))))
     &             omegaj(it,j) = c0*xf + c1
                if (xf.ge.xa(ne)) then
                  if (itype.eq.1) then
                    if (lbeth) then
                      omegaj(it,j) = bthp*dlog(xf-1d0+dexp(1d0))-a
                    else
                      if (lhighe) then
                        omegaj(it,j) = a*dlog(xf)+b
                      else
                        omegaj(it,j) = a*dlog(xf+b)
                      endif
                    endif
                  elseif (itype.eq.2) then
                    if (lbeth) then
                      omegaj(it,j) = b*(xf-1d0)/(xf-1d0+dexp(1d0))+a
                    else
                      if (lhighe) then
                        omegaj(it,j) = oma(ne)
                      else
                        omegaj(it,j) = a+b/xf
                      endif
                    endif
                  elseif (itype.eq.3) then
                    omegaj(it,j) = a/(xf+b)**2d0
                  endif
                endif
              enddo
            enddo
          enddo

        endif

C-----------------------------------------------------------------------
C  loop over energy points
C  approximate collision strength as first order power series
C-----------------------------------------------------------------------

        do 380 it = 1,ntf
          do 370 i = 1,nef-1

C-----------------------------------------------------------------------
C  SET UP INTERPOLATION PARAMETERS
C  omega = w0*energy + w1
C  fj = a1*exp(-a0*energy)
C  fi = b1*exp(-b0*energy)
C-----------------------------------------------------------------------

            v0 = 1.d0/(en(it,i)-en(it,i+1))
            v1 = 1.d0/(en(it,i+1)-en(it,i))
            a0 = v1*dlog(fi(it,i)/fi(it,i+1))
            a1 = fi(it,i)*(fi(it,i)/fi(it,i+1))**(eni(it,i)*v1)
            b0 = v1*dlog(fj(it,i)/fj(it,i+1))
            b1 = fj(it,i)*(fj(it,i)/fj(it,i+1))**(en(it,i)*v1)
            w0 = omegai(it,i)*v0 + omegai(it,i+1)*v1
            w1 = -eni(it,i+1)*omegai(it,i)*v0
     &           -eni(it,i)*omegai(it,i+1)*v1
            y0 = omegaj(it,i)*v0 + omegaj(it,i+1)*v1
            y1 = -en(it,i+1)*omegaj(it,i)*v0
     &           -en(it,i)*omegaj(it,i+1)*v1

C-----------------------------------------------------------------------
C  add in component part of the rate parameter in this interval
C-----------------------------------------------------------------------

            if (eni(it,i).ge.evt) then
              sumi(1) = dexp(evt/tvf(it)) * dsqrt(tvf(it))*a1/a0 * (
     &                  dexp(-a0*eni(it,i))*(w1+w0*eni(it,i)+w0/a0)
     &                - dexp(-a0*eni(it,i+1))*(w1+w0*eni(it,i+1)+w0/a0))
              if (sumi(1).lt.0d0) sumi(1) = 0.d0
              sum(1,it)=sum(1,it)+sumi(1)
            endif
            sumi(2) = dsqrt(tvf(it))*b1/b0 * (
     &                  dexp(-b0*en(it,i))*(y1+y0*en(it,i)+y0/b0)
     &                - dexp(-b0*en(it,i+1))*(y1+y0*en(it,i+1)+y0/b0))
            if (sumi(2).lt.0d0) sumi(2) = 0.d0
            sum(2,it)=sum(2,it)+sumi(2)
  370     enddo

C-----------------------------------------------------------------------
C  Add in part of the rate parameter from last point to infinity
C  Take last point as maximum of adf37 energy and adf04 energy
C-----------------------------------------------------------------------

          estep = (eva(ne)-eva(ne-1))/10d0
          uiu1 = dmax1(eva(ne),eni(it,nef))
          uju1 = dmax1(eva(ne),en(it,nef))
          intgl_tot1 = 0.d0
          intgl_tot2 = 0.d0
          do ie = 1,max
            uiu  = uiu1
            uiu1 = uiu1+estep
            uju  = uju1
            uju1 = uju1+estep
            if (itype.eq.1) then
              if (nform2.eq.2) then
                if (lbeth) then
                  intgl1 = f(it,nef)*en(it,nef)**param2(1)*
     &                   ((bthp*dlog(uiu1/evt-1d0+dexp(1d0))-a)*uiu1**
     &                        (-param2(1)-.5d0) +
     &                     (bthp*dlog(uiu/evt-1d0+dexp(1d0))-a)*uiu**
     &                        (-param2(1)-.5d0))*estep/2d0
                  intgl2 = f(it,nef)*en(it,nef)**param2(1)*
     &                    ((bthp*dlog(uju1/evt+dexp(1d0))-a)*uju1**
     &                        (-param2(1)-.5d0) +
     &                     (bthp*dlog(uju/evt+dexp(1d0))-a)*uju**
     &                        (-param2(1)-.5d0))*estep/2d0
                else
                  if (lhighe) then
                    intgl1 = f(it,nef)*en(it,nef)**param2(1)*
     &                       ((a*dlog(uiu1/evt)+b)*uiu1**
     &                          (-param2(1)-.5d0) +
     &                        (a*dlog(uiu/evt)+b)*uiu**
     &                          (-param2(1)-.5d0))*estep/2d0
                    intgl2 = f(it,nef)*en(it,nef)**param2(1)*
     &                       ((a*dlog(uju1/evt+1d0)+b)*uju1**
     &                          (-param2(1)-.5d0) +
     &                        (a*dlog(uju/evt+1d0)+b)*uju**
     &                          (-param2(1)-.5d0))*estep/2d0
                  else
                    intgl1 = f(it,nef)*en(it,nef)**param2(1)*
     &                       ((a*dlog(uiu1/evt+b))*uiu1**
     &                          (-param2(1)-.5d0) +
     &                        (a*dlog(uiu/evt+b))*uiu**
     &                          (-param2(1)-.5d0))*estep/2d0
                    intgl2 = f(it,nef)*en(it,nef)**param2(1)*
     &                       ((a*dlog(uju1/evt+1d0+b))*uju1**
     &                          (-param2(1)-.5d0) +
     &                        (a*dlog(uju/evt+1d0+b))*uju**
     &                          (-param2(1)-.5d0))*estep/2d0
                  endif
                endif
              elseif (nform2.eq.3) then
                if (lbeth) then
                  intgl1 = f(it,nef)*dexp(en(it,nef)*param2(1)/
     &                     tvf(it))*estep/2d0*
     &                    ((bthp*dlog(uiu1/evt-1d0+dexp(1d0))-a)*
     &                        dexp(-uiu1*param2(1)/tvf(it))/dsqrt(uiu1)+
     &                     (bthp*dlog(uiu/evt-1d0+dexp(1d0))-a)*
     &                        dexp(-uiu*param2(1)/tvf(it))/dsqrt(uiu))
                  intgl2 = f(it,nef)*dexp(en(it,nef)*param2(1)/
     &                     tvf(it))*estep/2d0*
     &                    ((bthp*dlog(uju1/evt+dexp(1d0))-a)*
     &                        dexp(-uju1*param2(1)/tvf(it))/dsqrt(uju1)+
     &                     (bthp*dlog(uju/evt+dexp(1d0))-a)*
     &                        dexp(-uju*param2(1)/tvf(it))/dsqrt(uju))
                else
                  if (lhighe) then
                    intgl1 = f(it,nef)*dexp(en(it,nef)*param2(1)/
     &                       tvf(it))*estep/2d0*
     &                       ((a*dlog(uiu1/evt)+b)*dexp(-uiu1*
     &                          param2(1)/tvf(it))/dsqrt(uiu1) +
     &                        (a*dlog(uiu/evt)+b)*dexp(-uiu*
     &                          param2(1)/tvf(it))/dsqrt(uiu))
                    intgl2 = f(it,nef)*dexp(en(it,nef)*param2(1)/
     &                       tvf(it))*estep/2d0*
     &                       ((a*dlog(uju1/evt+1d0)+b)*dexp(-uju1*
     &                          param2(1)/tvf(it))/dsqrt(uju1) +
     &                        (a*dlog(uju/evt+1d0)+b)*dexp(-uju*
     &                          param2(1)/tvf(it))/dsqrt(uju))
                  else
                    intgl1 = f(it,nef)*dexp(en(it,nef)*param2(1)/
     &                       tvf(it))*
     &                       ((a*dlog(uiu1/evt+b))*dexp(-uiu1*
     &                          param2(1)/tvf(it)) +
     &                        (a*dlog(uiu/evt+b))*dexp(-uiu*
     &                          param2(1)/tvf(it)))*estep/2d0
                    intgl2 = f(it,nef)*dexp(en(it,nef)*param2(1)/
     &                       tvf(it))*
     &                       ((a*dlog(uju1/evt+1d0+b))*dexp(-uju1*
     &                          param2(1)/tvf(it)) +
     &                        (a*dlog(uju/evt+1d0+b))*dexp(-uju*
     &                          param2(1)/tvf(it)))*estep/2d0
                  endif
                endif
              else
                intgl1 = 0d0
                intgl2 = 0d0
              endif
            elseif (itype.eq.2) then
              if (nform2.eq.2) then
                if (lbeth) then
                  intgl1 = f(it,nef)*en(it,nef)**param2(1)*estep/2d0*(
     &                     (b*(uiu1/evt-1d0)/(uiu1/evt-1d0+dexp(1d0))+a)
     &                       *uiu1**(-param2(1)-.5d0) +
     &                     (b*(uiu/evt-1d0)/(uiu/evt-1d0+dexp(1d0))+a)
     &                       *uiu**(-param2(1)-.5d0) )
                  intgl2 = f(it,nef)*en(it,nef)**param2(1)*estep/2d0*(
     &                     (b*(uju1/evt)/(uju1/evt+dexp(1d0))+a)
     &                       *uju1**(-param2(1)-.5d0) +
     &                     (b*(uju/evt)/(uju/evt+dexp(1d0))+a)
     &                       *uju**(-param2(1)-.5d0) )
                else
                  if (lhighe) then
                    intgl1 = f(it,nef)*en(it,nef)**param2(1)*oma(ne)*
     &                       (uiu1**(-param2(1)-.5d0) +
     &                        uiu**(-param2(1)-.5d0))*estep/2d0
                    intgl2 = f(it,nef)*en(it,nef)**param2(1)*oma(ne)*
     &                       (uju1**(-param2(1)-.5d0) +
     &                        uju**(-param2(1)-.5d0))*estep/2d0
                  else
                    intgl1 = f(it,nef)*en(it,nef)**param2(1)*estep/2d0*
     &                       ((a+b*evt/uiu1)*uiu1**(-param2(1)-.5d0) +
     &                        (a+b*evt/uiu)*uiu**(-param2(1)-.5d0))
                    intgl2 = f(it,nef)*en(it,nef)**param2(1)*estep/2d0*
     &                       ((a+b/(uju1/evt+1d0))*
     &                          uju1**(-param2(1)-.5d0) +
     &                        (a+b/(uju/evt+1d0))*
     &                          uju**(-param2(1)-.5d0))
                  endif
                endif
              elseif (nform2.eq.3) then
                if (lbeth) then
                  intgl1 = f(it,nef)*dexp(en(it,nef)*param2(1)
     &                     /tvf(it))*estep/2d0*(
     &                     (b*(uiu1/evt-1d0)/(uiu1/evt-1d0+dexp(1d0))+a)
     &                     *dexp(-uiu1*param2(1)/tvf(it))/dsqrt(uiu1)+
     &                     (b*(uiu/evt-1d0)/(uiu/evt-1d0+dexp(1d0))+a)
     &                     *dexp(-uiu*param2(1)/tvf(it))/dsqrt(uiu) )
                  intgl2 = f(it,nef)*dexp(en(it,nef)*param2(1)
     &                     /tvf(it))*estep/2d0*(
     &                     (b*(uju1/evt)/(uju1/evt+dexp(1d0))+a)
     &                     *dexp(-uju1*param2(1)/tvf(it))/dsqrt(uju1)+
     &                     (b*(uju/evt)/(uju/evt+dexp(1d0))+a)
     &                     *dexp(-uju*param2(1)/tvf(it))/dsqrt(uju) )
                else
                  if (lhighe) then
                    intgl1 = f(it,nef)*dexp(en(it,nef)*param2(1)
     &                       /tvf(it))*oma(ne)*estep/2d0*
     &                       (dexp(-uiu1*param2(1)/tvf(it))/dsqrt(uiu1)+
     &                       dexp(-uiu*param2(1)/tvf(it))/dsqrt(uiu))
                    intgl2 = f(it,nef)*dexp(en(it,nef)*param2(1)
     &                       /tvf(it))*oma(ne)*estep/2d0*
     &                       (dexp(-uju1*param2(1)/tvf(it))/dsqrt(uju1)+
     &                       dexp(-uju*param2(1)/tvf(it))/dsqrt(uju))
                  else
                    intgl1 = f(it,nef)*dexp(en(it,nef)*param2(1)
     &                       /tvf(it))*estep/2d0*
     &                       ((a+b*evt/uiu1)*
     &                       dexp(-uiu1*param2(1)/tvf(it))/dsqrt(uiu1)+
     &                       (a+b*evt/uiu)*
     &                       dexp(-uiu*param2(1)/tvf(it))/dsqrt(uiu))
                    intgl2 = f(it,nef)*dexp(en(it,nef)*param2(1)
     &                       /tvf(it))*estep/2d0*
     &                       ((a+b/(uju1/evt+1d0))*
     &                       dexp(-uju1*param2(1)/tvf(it))/dsqrt(uju1)+
     &                       (a+b/(uju/evt+1d0))*
     &                       dexp(-uju*param2(1)/tvf(it))/dsqrt(uju))
                  endif
                endif
              else
                intgl1 = 0d0
                intgl2 = 0d0
              endif
            elseif (itype.eq.3) then
              if (nform2.eq.2) then
                intgl1 = f(it,nef)*en(it,nef)**param2(1)*a*estep/2d0*
     &                   (uiu1**(-param2(1)-.5d0)/(uiu1/evt+b)**2d0 +
     &                   uiu**(-param2(1)-.5d0)/(uiu/evt+b)**2d0)
                intgl2 = f(it,nef)*en(it,nef)**param2(1)*a*estep/2d0*
     &                   (uju1**(-param2(1)-.5d0)/(uju1/evt+1d0+b)**2d0+
     &                   uju**(-param2(1)-.5d0)/(uju/evt+1d0+b)**2d0)
              elseif (nform2.eq.3) then
                intgl1 = f(it,nef)*dexp(en(it,nef)*param2(1)/tvf(it))*
     &                   a*(dexp(-uiu1*param2(1)/tvf(it))/dsqrt(uiu1)/
     &                   (uiu1/evt+b)**2d0 +
     &                   dexp(-uiu*param2(1)/tvf(it))/dsqrt(uiu)/
     &                   (uiu/evt+b)**2d0)*estep/2d0
                intgl2 = f(it,nef)*dexp(en(it,nef)*param2(1)/tvf(it))*
     &                   a*(dexp(-uju1*param2(1)/tvf(it))/dsqrt(uju1)/
     &                   (uju1/evt+1d0+b)**2d0 +
     &                   dexp(-uju*param2(1)/tvf(it))/dsqrt(uju)/
     &                   (uju/evt+1d0+b)**2d0)*estep/2d0
              else
                intgl1 = 0d0
                intgl2 = 0d0
              endif
            else
              intgl1 = 0d0
              intgl2 = 0d0
            endif

            intgl_tot1 = intgl_tot1 + intgl1
            intgl_tot2 = intgl_tot2 + intgl2
            if (intgl_tot1.gt.0.0) then
               if (intgl1/intgl_tot1.lt.tol) goto 850
            endif
            estep=estep*1.1d0
          enddo
  850     continue

          sumu(1) = dexp(evt/tvf(it))*dsqrt(tvf(it))*intgl_tot1
          sumu(2) = dsqrt(tvf(it))*intgl_tot2
          if (sumu(1).lt.0d0) sumu(1) = 0.d0
          if (sumu(2).lt.0d0) sumu(2) = 0.d0
          sum(1,it)=sum(1,it)+sumu(1)
          sum(2,it)=sum(2,it)+sumu(2)

C-----------------------------------------------------------------------
C  Add in part of the rate parameter from threshold to first point.
C  Threshold behaviour determined by itypt
C-----------------------------------------------------------------------

          if (itypt.eq.1.and.eva(1).gt.evt.and.eni(it,nef).lt.evt) then
            if (nform2.eq.2) then
              suml(1) = dexp(evt/tvf(it))*oma(1)*f(it,nef)*
     &                  en(it,nef)**param2(1)*dsqrt(tvf(it))/
     &                  (.5d0-param2(1))*
     &                  (eva(1)**(.5d0-param2(1))-evt**(.5d0-param2(1)))
            elseif (nform2.eq.3) then
              if (param2(1)*eva(1)/tvf(it).lt.1.5d0) then
                suml(1) = dexp(evt/tvf(it))*oma(1) *
     &                    tvf(it)/dsqrt(param2(1))*f(it,nef)*
     &                    dexp(param2(1)*en(it,nef)/tvf(it)+
     &                    lngama(.5d0))*(ingama(.5d0,param2(1)*eva(1)/
     &                    tvf(it))-ingama(.5d0,param2(1)*evt/tvf(it)))
              else
                suml(1) = dexp(evt/tvf(it))*oma(1) *
     &                    tvf(it)/dsqrt(param2(1))*f(it,nef)*
     &                    dexp(param2(1)*en(it,nef)/tvf(it)+
     &                    lngama(.5d0))*(ingamq(.5d0,param2(1)*evt/
     &                    tvf(it))-ingamq(.5d0,param2(1)*eva(1)/
     &                    tvf(it)))
              endif
            endif
          else
            suml(1) = 0d0
          endif
          if (suml(1).lt.0d0) suml(1) = 0.d0
          sum(1,it)=sum(1,it)+suml(1)

C-----------------------------------------------------------------------
C  add in part of rate parameter from last adf37 energy to last adf04
C  energy (if last 37 energy < last 04 energy)
C-----------------------------------------------------------------------

          do 180 i = 1,ne-1

c-----------------------------------------------------------------------
c  set up interpolation parameters
c-----------------------------------------------------------------------

            v0 = 1.d0/(eva(i)-eva(i+1))
            v1 = 1.d0/(eva(i+1)-eva(i))
            w0 = oma(i)*v0 + oma(i+1)*v1
            w1 = -eva(i+1)*oma(i)*v0-eva(i)*oma(i+1)*v1
            j=0
  800       continue
            j=j+1
            if (eva(i)+evt.lt.eva(j).and.j.lt.ne-1) goto 800
            y0 = oma(j)*v0 + oma(j+1)*v1
            y1 = -eva(i+1)*oma(j)*v0-eva(i)*oma(j+1)*v1

C-----------------------------------------------------------------------
C  add in component part of the rate parameter in this interval
C-----------------------------------------------------------------------

            if (eni(it,nef).lt.eva(i+1)) then
              if (nform2.eq.2) then
                sumx(1) = dexp(evt/tvf(it)) * en(it,nef)**param2(1)*
     &                    f(it,nef)*dsqrt(tvf(it))*
     &                    (w0/(1.5d0-param2(1))*(eva(i+1)**
     &                    (1.5d0-param2(1))-dmax1(eva(i),eni(it,nef))**
     &                    (1.5d0-param2(1)))
     &                    + w1/(.5d0-param2(1))*(eva(i+1)**
     &                    (0.5d0-param2(1))-eva(i+1)**
     &                    (0.5d0-param2(1))))
              elseif (nform2.eq.3) then
                if (param2(1)/tvf(it)*eva(i+1).lt.2.5d0) then
                  igmxi1 = dexp(lngama(1.5d0))*
     &                     (ingama(1.5d0,param2(1)/tvf(it)*eva(i+1)) -
     &                     ingama(1.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),eni(it,nef))))
                else
                  igmxi1 = dexp(lngama(1.5d0))*
     &                     (-ingamq(1.5d0,param2(1)/tvf(it)*eva(i+1)) +
     &                     ingamq(1.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),eni(it,nef))))
                endif
                if (param2(1)/tvf(it)*eva(i+1).lt.1.5d0) then
                  igmxi2 = dexp(lngama(0.5d0))*
     &                     (ingama(0.5d0,param2(1)/tvf(it)*eva(i+1)) -
     &                     ingama(0.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),eni(it,nef))))
                else
                  igmxi2 = dexp(lngama(0.5d0))*
     &                     (-ingamq(0.5d0,param2(1)/tvf(it)*eva(i+1)) +
     &                     ingamq(0.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),eni(it,nef))))
                endif
                sumx(1) = dexp(evt/tvf(it)) *
     &                    tvf(it)/param2(1)*f(it,nef)*
     &                    dexp(param2(1)*en(it,nef)/tvf(it))*
     &                    (w0*tvf(it)/dsqrt(param2(1))*igmxi1 +
     &                    w1*dsqrt(param2(1))*igmxi2 )
              else
                sumx(1) = 0.d0
              endif
              if (sumx(1).lt.0d0) sumx(1) = 0.d0
              sum(1,it)=sum(1,it)+sumx(1)

            endif
            if (en(it,nef)+evt.lt.eva(i+1)) then
              if (nform2.eq.2) then
                sumx(2) = en(it,nef)**param2(1)*
     &                    f(it,nef)*dsqrt(tvf(it))*
     &                    (y0/(1.5d0-param2(1))*(eva(i+1)**
     &                    (1.5d0-param2(1))-dmax1(eva(i),en(it,nef))
     &                    **(1.5d0-param2(1))) + y1/
     &                    (.5d0-param2(1))*(eva(i+1)**
     &                    (0.5d0-param2(1))-eva(i+1)**
     &                    (0.5d0-param2(1))))
              elseif (nform2.eq.3) then
                if (param2(1)/tvf(it)*(eva(i+1)-evt).lt.2.5d0) then
                  igmxj1 = dexp(lngama(1.5d0))*
     &                     (ingama(1.5d0,param2(1)/tvf(it)*
     &                     eva(i+1)) -
     &                     ingama(1.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),en(it,nef))))
                else
                  igmxj1 = dexp(lngama(1.5d0))*
     &                     (-ingamq(1.5d0,param2(1)/tvf(it)*
     &                     eva(i+1)) +
     &                     ingamq(1.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),en(it,nef))))
                endif
                if (param2(1)/tvf(it)*(eva(i+1)-evt).lt.1.5d0) then
                  igmxj2 = dexp(lngama(0.5d0))*
     &                     (ingama(0.5d0,param2(1)/tvf(it)*
     &                     eva(i+1)) -
     &                     ingama(0.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),en(it,nef))))
                else
                  igmxj2 = dexp(lngama(0.5d0))*
     &                     (-ingamq(0.5d0,param2(1)/tvf(it)*
     &                     eva(i+1)) +
     &                     ingamq(0.5d0,param2(1)/tvf(it)*
     &                     dmax1(eva(i),en(it,nef))))
                endif
                sumx(2) = tvf(it)/param2(1)*f(it,nef)*
     &                    dexp(param2(1)*en(it,nef)/tvf(it))*
     &                    (y0*tvf(it)/dsqrt(param2(1))*igmxj1 +
     &                    y1*dsqrt(param2(1))*igmxj2 )
              else
                sumx(2) = 0.d0
              endif
              if (sumx(2).lt.0d0) sumx(2) = 0.d0
              sum(2,it)=sum(2,it)+sumx(2)
            endif

C-----------------------------------------------------------------------
C  add in part of rate parameter from first adf04 energy to first adf37
C  energy (if first 37 energy > first 04 energy)
C-----------------------------------------------------------------------

            if (nform1.eq.2) then
              if (eni(it,1).gt.eva(i)) then
                sumx(1) = dexp(evt/tvf(it))*
     &                    f(it,1)*en(it,1)**(-param1)*dsqrt(tvf(it))*
     &                    (w0/(param1+1.5d0)*
     &                      (dmin1(eva(i+1),eni(it,1))-eva(i)) +
     &                    w1/(param1+.5d0)*
     &                      (dmin1(eva(i+1),eni(it,1))-eva(i)))
              else
                sumx(1) = 0d0
              endif
              if (sumx(1).lt.0d0) sumx(1) = 0.d0
              sum(1,it)=sum(1,it)+sumx(1)
              if (en(it,1)+evt.gt.eva(i)) then
                sumx(2) = f(it,1)*en(it,1)**(-param1)*dsqrt(tvf(it))*
     &                    (y0/(param1+1.5d0)*
     &                      (dmin1(eva(i+1),en(it,1)+evt)-eva(i)) +
     &                    y1/(param1+.5d0)*
     &                      (dmin1(eva(i+1),en(it,1)+evt)-eva(i)))
              else
                sumx(2) = 0d0
              endif
              if (sumx(2).lt.0d0) sumx(2) = 0.d0
              sum(2,it)=sum(2,it)+sumx(2)
            endif


  180     enddo

  380   enddo
      endif

C-----------------------------------------------------------------------
C  convert sum into upsilon.
C  output results for upsilon (excitation and de-excitation).
C-----------------------------------------------------------------------

      if (dist.eq.1) then
        k1 = lngama(kappa+1.0d0)-lngama(kappa-0.5d0)
        k2 = dsqrt(2.0d0/(2.0d0*kappa-3.0d0))*dexp(k1)
      endif

      do 95 it = 1 , nt
        if     (dist .eq. 0) then
          upsilon(it) = sum(1,it)
          dnsilon(it) = upsilon(it)
        elseif (dist .eq. 1) then
          upsilon(it) = dexp(evt/tva(it))*k2*sum(1,it)
          dnsilon(it) = k2*sum(2,it)
        elseif (dist .eq. 3) then
          upsilon(it) = dexp(evt/tva(it))*
     &                  sum(1,it) * dsqrt(pi/6.d0) *
     &                  dexp(0.5d0*lngama(2.5d0/dru_val)-
     &                       1.5d0*lngama(1.5d0/dru_val))
          dnsilon(it) = sum(2,it) * dsqrt(pi/6.d0) *
     &                  dexp(0.5d0*lngama(2.5d0/dru_val)-
     &                       1.5d0*lngama(1.5d0/dru_val))
        endif
        if (upsilon(it).lt.1d-30) upsilon(it) = 1.d-30
        if (dnsilon(it).lt.1d-30) dnsilon(it) = 1.d-30
   95 continue

      if (dist.eq.2) then
        do 900 it = 1 , ntf
          upsilon(it) = dsqrt(pi)/2.d0*sum(1,it)
          dnsilon(it) = dsqrt(pi)/2.d0*sum(2,it)
  900   continue
      endif
      
c      write(i4unit(-1),*)'upsilon=',(upsilon(it),it=1,nt)
c      write(i4unit(-1),*)'dnsilon=',(dnsilon(it),it=1,nt)

      return
C-----------------------------------------------------------------------
 1011 format(1x,31('*'),' h9ntqd warning ',31('*')//
     &       1x,'Fault in input adf04 file: ',a,i7)
 1013 format(1x,a, i4)
C-----------------------------------------------------------------------
      end
