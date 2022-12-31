       subroutine h9qd3b( nedim   , ntdim   , nfdim  ,
     &                    wupper  ,
     &                    ne      , nt      ,
     &                    evt     , xa      ,
     &                    omega   , tva     ,
     &                    kap_val , dru_val , dist   ,
     &                    nef     , en      , f      ,
     &                    nte     , tenum   ,
     &                    nform1  , param1  , nform2 , param2 ,
     &                    alpha   , q
     &                  )

      implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 subroutine: h9qd3b *********************
C
C  purpose:  To execute quadratures over ionisation collision strengths
C            to determine the ionisation and 3-body recombination
C            coefficients. Free electron distribution function may be
C            Maxwellian, Kappa, Druyvesteyn, or numeric from adf37 file.
C
C  calling program: adas809
C
C  input : (i*4)  nedim    = max no of energy points that can be read in
C  input : (i*4)  ntdim    = max no of temperatures
C          (i*4)  nfdim    = max no of energies in adf37 file
C  input : (r*8)  tva()   = input effective temperatures (eV)
C  input : (r*8)  tenum()  = effective temperatures from adf37 file (eV)
C  input : (i*4)  ne       = no of x parameter values in adf04 type 1
C  input : (i*4)  nt       = no of input temperatures
C  input : (r*8)  wupper   = statistical weight of upper level
C  input : (r*8)  wlower   = statistical weight of lower level
C  input : (r*8)  xa()     = x parameter from adf04 type1
C  input : (i*4)  dist     = distribution type
C                              0 => Maxwellian
C                              1 => kappa
C                              2 => numeric
C                              3 => Druyvesteyn
C  input : (r*8)  kap_val  = value of kappa parameter
C  input : (r*8)  dru_val  = value of x parameter in Druyvesteyn dist
C  input : (r*8)  evt      = ionisation potential
C  input : (r*8)  omega()  = collision strength from adf04 type1
C  input : (i*4)  nef      = no of energy points in adf37
C  input : (r*8)  en(,)    = energy points of distribution tabulation
C  input : (r*8)  f(,)     = distribution function tabulation
C  input : (i*4)  nte      = no of temperatures in adf37
C  input : (i*4)  nform1   = type of threshold behaviour
C                              1 => cutoff
C                              2 => energy^param1
C  input : (r*8)  param1() = parameter of threshold form
C  input : (i*4)  nform2   = type of high-energy behaviour
C                              1 => cutoff
C                              2 => energy^-param2(1)
C                              3 => exp(-param2(1)*energy)
C  input : (r*8)  param2   = parameter of high-energy form
C
C  output: (r*8)  alpha    = 3-body recombination coefficient
C  output: (r*8)  q        = ionisation coefficient
C
C  routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          eei        ADAS      evaluates  exp(x)*E1(x)
C          lngama     ADAS      evaluates ln(gamma(x))
C
C  author: Paul Bryans, University of Strathclyde
C
C  date:   22/11/04
C
c  modification history:
c
c  version : 1.1
c  date    : 22/11/04                     
c  modified: Paul Bryans
c               - first release
c
c  version : 1.2
c  date    : 26/11/04
c  modified: Allan Whiteford
c               - made some arrays locally dimensioned to prevent
c                 g77 compiled code failing when it tried to make
c                 large automatic arrays.
c
c  version : 1.3
c  date    : 26/11/04
c  modified: Paul Bryans
c               - Moved final exponential coefficient of Maxwellian
c                 and Druyvesteyn ionisation rate inside integrals to
c                 eliminate overflow problem when this becomes large.
c
c  version : 1.4
c  date    : 06/12/04
c  modified: Paul Bryans
c               - Added in line which had been accidentally removed
c                 in last version.
c               - Changed dimensioning and limits slightly.
c               - Removed some commented out code.
c               - Changes by Paul Bryans, these comments by
c                 Allan Whiteford
c
c  version : 1.5
c  date    : 25-10-2016
c  modified: Martin O'Mullane
c             - Re-order arguments to match order in h9ntqd.
c             - Tidy code to ADAS standards. 
c
C-----------------------------------------------------------------------
      integer     npts  , limit     , limit1    , qdne  , qdnt  , qdnf
C-----------------------------------------------------------------------
      real*8      tol   , dzero
      real*8      fine  , a0        , ryd       , c
C-----------------------------------------------------------------------
      parameter ( npts = 10 , limit = 20 , limit1 = 20  ,
     &            qdne = 50 , qdnt  = 50 , qdnf   = 1500 )
C-----------------------------------------------------------------------
      parameter ( tol  = 1.0D-10    , dzero = 1.0d-30 )
      parameter ( fine = 7.2974D-03 , a0    = 5.2918D-09 ,
     &            ryd  = 13.6048D0  , c     = 2.9979D+10)
C-----------------------------------------------------------------------
      integer    it, nt, ne, ie1, ie2, dist, nedim, ntdim, i, nef
      integer    nfdim, i4unit
      integer    ntemp, nform1, nform2, ie3, ie4, nte, ie
C-----------------------------------------------------------------------
      logical    lfirst
C-----------------------------------------------------------------------
      integer    limit2(ntdim)
C-----------------------------------------------------------------------
      real*8     evt     , pi      , sum2    , kek   , ex      , 
     &           wupper  , kap_val , c0      , c1    , param1  , 
     &           sumi    , sumu    , d1      , d0    , energy2 
      real*8     eei     , lngama
C-----------------------------------------------------------------------
      real*8     omega(nedim) , f(ntdim,nfdim) , en(ntdim,nfdim)
      real*8     alpha(ntdim), q(ntdim)
      real*8     dru_val, tva(ntdim), xa(nedim)
      real*8     oma(ntdim,qdnf+2*limit1)
      real*8     ea(nedim), ut(ntdim), ui(ntdim), ui1(ntdim)
      real*8     sum(qdnt)
      real*8     tenum(qdnt) , te(qdnt)
      real*8     param2(2), c2, sumx, sumf
      real*8     omg(qdne*npts+limit), estep
      real*8     sum1(qdne*npts) ,energy(qdne*npts+limit)
      real*8     fn(qdnt,qdnf+2*limit1), enrg(qdnt,qdnf+2*limit1)
      real*8     dist1(qdnt,qdnf+2*limit1)
      real*8     dist2(qdnt,qdnf+2*limit1,qdnf+2*limit1)
      real*8     int1(qdnf*npts+limit,qdnf*npts+limit)
	  real*8	 temp1, temp2, temp3
C-----------------------------------------------------------------------
      external   eei  , lngama 
C-----------------------------------------------------------------------


      if (qdne.ne.nedim) then
        write(i4unit(-1),1000) 'QDNE does not match NEDIM'
        stop
      endif
      if (qdnf.ne.nfdim) then
        write(i4unit(-1),1000) 'QDNF does not match NFDIM'
        stop
      endif
      if (qdnt.ne.ntdim) then
        write(i4unit(-1),1000) 'QDNT does not match NTDIM'
        stop
      endif

C-----------------------------------------------------------------------
C Setup
C-----------------------------------------------------------------------

      pi=acos(-1.d0)

C-----------------------------------------------------------------------
C convert x parameter to incident energy in ev
C-----------------------------------------------------------------------

      do i = 1,ne
        ea(i) = xa(i)*(evt+dzero)
      end do

C-----------------------------------------------------------------------
C  For numerical distribution, set number of temperatures to that in
C  adf37. For analytic distribution, set to nt.
C-----------------------------------------------------------------------

      if (dist.eq.2) then
        ntemp = nte
      else
        ntemp = nt
      endif

C-----------------------------------------------------------------------
C  Initialise temperature and energy arrays.
C-----------------------------------------------------------------------

      do it = 1,ntemp
        sum(it) = 0d0
        if (dist.eq.2) then
          te(it) = tenum(it)
          ui(it) = en(it,1)/te(it)
        else
          te(it) = tva(it)
          ui(it) = ea(1)/te(it)
        endif
        ut(it) = evt/te(it)
      enddo

C-----------------------------------------------------------------------
C  Calculate ionisation coefficient for Maxwellian distribution.
C-----------------------------------------------------------------------
C  Set omega = c0*ln(energy) + c1
C-----------------------------------------------------------------------

      if (dist.eq.0) then

        do 30 ie1 = 1,ne-1
          do 40 it = 1,ntemp
            ui1(it) = ea(ie1+1)/te(it)
            c1 = dlog(ui(it))*omega(ie1+1)/(dlog(ui(it))-dlog(ui1(it)))+
     &           dlog(ui1(it))*omega(ie1)/(dlog(ui1(it))-dlog(ui(it)))
            c0 = omega(ie1)/(dlog(ui(it))-dlog(ui1(it))) +
     &           omega(ie1+1)/(dlog(ui1(it))-dlog(ui(it)))
            if (ie1.eq.1) then
              sumi = (c0*(dlog(ut(it))+eei(ut(it)))+c1)
     &             - dexp((evt-ea(ie1+1))/te(it))*
     &               (c0*(dlog(ui1(it))+eei(ui1(it)))+c1)
            elseif (ie1.eq.ne-1) then
              sumi = dexp((evt-ea(ie1))/te(it))*
     &               (c0*(dlog(ui(it))+eei(ui(it)))+c1)
            else
              sumi = dexp((evt-ea(ie1))/te(it))*
     &                (c0*(dlog(ui(it))+eei(ui(it)))+c1)
     &             - dexp((evt-ea(ie1+1))/te(it))*
     &                (c0*(dlog(ui1(it))+eei(ui1(it)))+c1)
            endif
            sum(it) = sum(it) + sumi
   40     enddo
          do 50 it=1,ntemp
            ui(it) = ui1(it)
   50     enddo
   30   enddo

        do it=1,ntemp
          alpha(it) = sum(it)/dsqrt(te(it))
     &                       *2d0*fine*c*a0*a0/wupper*dsqrt(pi*ryd)
          q(it) = alpha(it)
        enddo

C-----------------------------------------------------------------------
C  Calculate ionisation coefficient for kappa distribution.
C-----------------------------------------------------------------------
C  omega = c0*ln(energy)+c1 cannot be integrated analytically over this
C  distribution. Increase number of omega tabulations by factor npts
C  and use linear interpolation on this finer grid.
C-----------------------------------------------------------------------

      elseif (dist.eq.1) then

        do 70 ie1 = 1,ne-1
          c1 = dlog(ea(ie1))*omega(ie1+1)/
     &         (dlog(ea(ie1))-dlog(ea(ie1+1)))+
     &         dlog(ea(ie1+1))*omega(ie1)/
     &         (dlog(ea(ie1+1))-dlog(ea(ie1)))
          c0 = omega(ie1)/(dlog(ea(ie1))-dlog(ea(ie1+1))) +
     &         omega(ie1+1)/(dlog(ea(ie1+1))-dlog(ea(ie1)))
          energy(1) = ea(ie1)
          omg(1) = omega(ie1)
          do 80 ie2 = 1,npts
            energy(ie2+1) = (ea(ie1+1)-ea(ie1))*ie2/npts + ea(ie1)
            omg(ie2+1) = c0*dlog(energy(ie2+1)) + c1
            do 90 it = 1,ntemp
              kek = te(it)*(kap_val-1.5d0)
              sumi = (energy(ie2+1)-energy(ie2))/2.d0*(omg(ie2)*
     &               (1.d0+energy(ie2)/kek)**(-kap_val-1.d0)+omg(ie2+1)*
     &               (1.d0+energy(ie2+1)/kek)**(-kap_val-1.d0))
              sum(it) = sum(it) + sumi*2d0*fine*c*a0*a0/wupper
     &                  *dexp(lngama(kap_val+1.d0)-lngama(kap_val-.5d0))
     &                  *kek**(-1.5d0)*dsqrt(pi*ryd)
   90       enddo
   80     enddo
   70   enddo

C-----------------------------------------------------------------------
C  Add contributions above the last tabulated omega value. Continue this
C  until contribution is less than tolerance (tol), or until a maximum
C  number of points (limit) is reached.
C-----------------------------------------------------------------------
 
        energy(1) = ea(ne)
        omg(1) = omega(ne)
        estep = (ea(ne)-ea(ne-1))/npts
        do 100 it = 1,ntemp
          kek = te(it)*(kap_val-1.5d0)
          do 110 ie = 1,limit
            energy(ie+1) = energy(ie)+estep
            omg(ie+1) = c0*dlog(energy(ie+1)) + c1
            sumi = (energy(ie+1)-energy(ie))/2.d0*(omg(ie)*
     &             (1.d0+energy(ie)/kek)**(-kap_val-1.d0)+omg(ie+1)*
     &             (1.d0+energy(ie+1)/kek)**(-kap_val-1.d0))
            sum(it) = sum(it) + sumi*2d0*fine*c*a0*a0/wupper
     &                *dexp(lngama(kap_val+1.d0)-lngama(kap_val-.5d0))
     &                *kek**(-1.5d0)*dsqrt(pi*ryd)
            if (sumi/sum(it).lt.tol) goto 120
  110     enddo
  120     continue
  100   enddo

C-----------------------------------------------------------------------
C  Calculate ionisation coefficient for numerical distribution.
C-----------------------------------------------------------------------

      elseif (dist.eq.2) then

C-----------------------------------------------------------------------
C  calculate contribution to integral between lower and upper energy
C  points in adf37 file
C-----------------------------------------------------------------------

        do 130 it = 1,ntemp
          do 140 ie1 = 1,nef-1
            if (en(it,ie1).ge.evt) then
              do 150 ie2 = 1,ne-1

C-----------------------------------------------------------------------
C  interpolate/extrapolate omega to adf37 energy points using
C  omega = c0*ln(energy)+c1
C-----------------------------------------------------------------------

                if ((ea(ie2).lt.en(it,ie1).and.ea(ie2+1).gt.en(it,ie1))
     &            .or.(ea(ie2).gt.en(it,ie1).and.ie2.eq.1)
     &            .or.(ea(ie2).lt.en(it,ie1).and.ie2.eq.ne-1)) then
                  c1 = dlog(ea(ie2))*omega(ie2+1)/
     &                 (dlog(ea(ie2))-dlog(ea(ie2+1)))+
     &                 dlog(ea(ie2+1))*omega(ie2)/
     &                 (dlog(ea(ie2+1))-dlog(ea(ie2)))
                  c0 = omega(ie2)/(dlog(ea(ie2))-dlog(ea(ie2+1))) +
     &                 omega(ie2+1)/(dlog(ea(ie2+1))-dlog(ea(ie2)))

C-----------------------------------------------------------------------
C  execute quadrature with distribution set to
C  f = d0*sqrt(energy)*exp(-d1*energy)
C-----------------------------------------------------------------------

                  if (f(it,ie1)/dsqrt(en(it,ie1)).ge.
     &                f(it,ie1+1)/dsqrt(en(it,ie1+1))) then
                    d1 = dlog(f(it,ie1)/f(it,ie1+1)*dsqrt(en(it,ie1+1)/
     &                   en(it,ie1)))/(en(it,ie1+1)-en(it,ie1))
                    d0 = f(it,ie1)/dsqrt(en(it,ie1))*(f(it,ie1)/
     &                   f(it,ie1+1)*dsqrt(en(it,ie1+1)/en(it,ie1)))**
     &                   (en(it,ie1)/(en(it,ie1+1)-en(it,ie1)))
                    sumi = (dexp(-d1*en(it,ie1))*(c0*(dlog(en(it,ie1))+
     &                     eei(d1*en(it,ie1)))+c1)
     &                     - dexp(-d1*en(it,ie1+1))*
     &                     (c0*(dlog(en(it,ie1+1))+
     &                     eei(d1*en(it,ie1+1)))+c1))/d1*d0

C-----------------------------------------------------------------------
C  if f/sqrt(energy) is not decreasing then execute quadrature with
C  f = d0*energy+d1
C-----------------------------------------------------------------------

                  else
                    d1 = (f(it,ie1)*en(it,ie1+1)/dsqrt(en(it,ie1))-
     &                   f(it,ie1+1)*en(it,ie1)/dsqrt(en(it,ie1+1)))
     &                   /(en(it,ie1+1)-en(it,ie1))
                    d0 = (f(it,ie1+1)/dsqrt(en(it,ie1+1))-
     &                   f(it,ie1)/dsqrt(en(it,ie1)))/
     &                   (en(it,ie1+1)-en(it,ie1))
                    sumi =d0*c0*(en(it,ie1+1)**(2.d0)*dlog(en(it,ie1+1))
     &                      - en(it,ie1)**(2.d0)*dlog(en(it,ie1)) )/2.d0
     &                   + d1*c0*( en(it,ie1+1)*dlog(en(it,ie1+1))
     &                           - en(it,ie1)*dlog(en(it,ie1)) )
     &                   + (d0/2.d0*(c1-c0/2.d0))*(en(it,ie1+1)**(2.d0)-
     &                      en(it,ie1)**(2.d0))
     &                   + (d1*(c1-c0)*(en(it,ie1+1)-en(it,ie1)))
                  endif
                  if (sumi.lt.0.d0) sumi = 0.d0
                  sum(it) = sum(it) + sumi
     &                      *fine*c*a0*a0/wupper*dsqrt(ryd)*pi
                endif
  150         enddo
            endif
  140     enddo
  130   enddo

C-----------------------------------------------------------------------
C  Add contribution from energies between highest adf37 energy and
C  highest adf04 type1 energy.
C  High energy behaviour of distribution:
C    nform2=1: f = 0
C    nform2=2: f = c2*energy^(-param2(1))
C    nform2=3: f = c2*sqrt(energy)*exp(-param2(1)*energy)
C    nform2=4: not set
C  Contributions from ie1 to ie1+1 in adf04 type1 energy space are
C  represented by sumx.
C  For first extra interval, also take into account contribution from
C  last adf37 energy to the first adf04 energy as sumf.
C-----------------------------------------------------------------------

        do 160 it = 1,ntemp
          lfirst = .true.
          do 170 ie1 = 1,ne-1
            sumx = 0.d0
            if (ea(ie1).gt.en(it,nef)) then
              c1 = dlog(ea(ie1))*omega(ie1+1)/
     &             (dlog(ea(ie1))-dlog(ea(ie1+1)))+
     &             dlog(ea(ie1+1))*omega(ie1)/
     &             (dlog(ea(ie1+1))-dlog(ea(ie1)))
              c0 = omega(ie1)/(dlog(ea(ie1))-dlog(ea(ie1+1))) +
     &             omega(ie1+1)/(dlog(ea(ie1+1))-dlog(ea(ie1)))
              if     (nform2.eq.1) then
                sumx = 0.d0
              elseif (nform2.eq.2) then
                c2 = f(it,nef)*en(it,nef)**param2(1)
                sumx = c2/(.5d0-param2(1))*(
     &                 ea(ie1+1)**(.5d0-param2(1))*
     &                 (c1+c0*(dlog(ea(ie1+1))+1.d0/(param2(1)-.5d0))) -
     &                 ea(ie1)**(.5d0-param2(1))*
     &                 (c1+c0*(dlog(ea(ie1))+1.d0/(param2(1)-.5d0))) )
                if (lfirst) then
                  sumf = c2/(.5d0-param2(1))*(
     &                   ea(ie1)**(.5d0-param2(1))*
     &                   (c1+c0*(dlog(ea(ie1))+1.d0/(param2(1)-.5d0))) -
     &                   dmax1(en(it,nef),evt)**(.5d0-param2(1))*
     &                   (c1+c0*(dlog(dmax1(en(it,nef),evt))+
     &                   1.d0/(param2(1)-.5d0))) )
                  if (sumf.gt.0.d0) sumx = sumx + sumf
                endif
              elseif (nform2.eq.3) then
                sumx = f(it,nef)/param2(1)*te(it)/dsqrt(en(it,nef))*
     &                 (dexp(-param2(1)/te(it)*(ea(ie1)-en(it,nef)))
     &          *(c0*(dlog(ea(ie1))+eei(param2(1)/te(it)*ea(ie1)))+c1) -
     &                 dexp(-param2(1)/te(it)*ea(ie1+1))
     &                 *(c0*(dlog(ea(ie1+1))+eei(param2(1)/te(it)*
     &                 ea(ie1+1)))+c1) )
                if (lfirst) then
                  sumf = f(it,nef)/param2(1)*te(it)/dsqrt(en(it,nef))*
     &                   (dexp(-param2(1)/te(it)*
     &                   (dmax1(en(it,nef),evt)-en(it,nef)))*
     &                   (c0*(dlog(dmax1(en(it,nef),evt))+
     &              eei(param2(1)/te(it)*dmax1(en(it,nef),evt)))+c1) -
     &                   dexp(-param2(1)/te(it)*ea(ie1))*
     &            (c0*(dlog(ea(ie1))+eei(param2(1)/te(it)*ea(ie1)))+c1))
                  if (sumf.gt.0.d0) sumx = sumx + sumf
                endif
              elseif (nform2.eq.4) then
                sumx = 0.d0
              endif
              lfirst = .false.

C-----------------------------------------------------------------------
C  Add contributions from adf04 energies that are lower than lowest
C  adf37 energy.
C  Threshold behaviour of distribution:
C    nform1=1: f = 0
C    nform1=2: f = c2*energy^param1
C  Contributions from ie1 to ie1+1 in adf04 type1 energy space are
C  represented by sumx.

C-----------------------------------------------------------------------
            elseif (ea(ie1+1).lt.en(it,1)) then
              c1 = dlog(ea(ie1))*omega(ie1+1)/
     &             (dlog(ea(ie1))-dlog(ea(ie1+1)))+
     &             dlog(ea(ie1+1))*omega(ie1)/
     &             (dlog(ea(ie1+1))-dlog(ea(ie1)))
              c0 = omega(ie1)/(dlog(ea(ie1))-dlog(ea(ie1+1))) +
     &             omega(ie1+1)/(dlog(ea(ie1+1))-dlog(ea(ie1)))
              if     (nform1.eq.1) then
                sumx = 0.d0
              elseif (nform1.eq.2) then
                c2 = f(it,1)/(en(it,1)**param1)
                sumx = c2/(.5d0+param1)*(
     &                 ea(ie1+1)**(.5d0+param1)*(c1+c0*
     &                 (dlog(ea(ie1+1))-1.d0/(param1+.5d0))) -
     &                 ea(ie1)**(.5d0+param1)*(c1+c0*
     &                 (dlog(ea(ie1))-1.d0/(param1+.5d0))))
              endif
            elseif (ea(ie1).lt.en(it,1)) then
              c1 = dlog(ea(ie1))*omega(ie1+1)/
     &             (dlog(ea(ie1))-dlog(ea(ie1+1)))+
     &             dlog(ea(ie1+1))*omega(ie1)/
     &             (dlog(ea(ie1+1))-dlog(ea(ie1)))
              c0 = omega(ie1)/(dlog(ea(ie1))-dlog(ea(ie1+1))) +
     &             omega(ie1+1)/(dlog(ea(ie1+1))-dlog(ea(ie1)))
              if     (nform1.eq.1) then
                sumx = 0.d0
              elseif (nform1.eq.2) then
                c2 = f(it,1)/(en(it,1)**param1)
                sumx = c2/(.5d0+param1)*(
     &                 en(it,1)**(.5d0+param1)*(c1+c0*
     &                 (dlog(en(it,1))-1.d0/(param1+.5d0))) -
     &                 ea(ie1)**(.5d0+param1)*(c1+c0*
     &                 (dlog(ea(ie1))-1.d0/(param1+.5d0))))
              endif
            else
              sumx = 0.d0
            endif
            if (sumx.lt.0.d0) sumx = 0.d0
            sum(it) = sum(it) + sumx
     &                  *fine*c*a0*a0/wupper*dsqrt(ryd)*pi

  170     enddo

C-----------------------------------------------------------------------
C  Add contributions from highest energy point (either adf04 or adf37)
C  to infinity.
C-----------------------------------------------------------------------

          c1 = dlog(ea(ne-1))*omega(ne)/
     &         (dlog(ea(ne-1))-dlog(ea(ne)))+
     &         dlog(ea(ne))*omega(ne-1)/
     &         (dlog(ea(ne))-dlog(ea(ne-1)))
          c0 = omega(ne-1)/(dlog(ea(ne-1))-dlog(ea(ne))) +
     &         omega(ne)/(dlog(ea(ne))-dlog(ea(ne-1)))
          if     (nform2.eq.1) then
            sumu = 0.d0
          elseif (nform2.eq.2) then
            c2 = f(it,nef)*en(it,nef)**param2(1)
            sumu = c2/(param2(1)-.5d0)*(
     &             dmax1(ea(ne),en(it,nef))**(.5d0-param2(1))*
     &             (c1+c0*(dlog(dmax1(ea(ne),en(it,nef)))+
     &             1.d0/(param2(1)-.5d0))))
          elseif (nform2.eq.3) then
            sumu = f(it,nef)/param2(1)*(dexp(-param2(1)*
     &             (dmax1(ea(ne),en(it,nef))-en(it,nef)))*
     &             (c0*(dlog(dmax1(ea(ne),en(it,nef)))+
     &             eei(param2(1)*dmax1(ea(ne),en(it,nef))))+c1))/
     &             dsqrt(en(it,nef))
          elseif (nform2.eq.4) then
            sumu = 0.d0
          endif
          if (sumu.lt.0.d0) sumu = 0.d0
          sum(it) = sum(it) + sumu
     &                *fine*c*a0*a0/wupper*dsqrt(ryd)*pi


  160   enddo

C-----------------------------------------------------------------------
C  Calculate ionisation coefficient for Druyvesteyn distribution.
C-----------------------------------------------------------------------
C  omega = c0*ln(energy)+c1 cannot be integrated analytically over this
C  distribution. Increase number of omega tabulations by factor npts
C  and use linear interpolation on this finer grid.
C-----------------------------------------------------------------------

      elseif (dist.eq.3) then

        do 180 ie1 = 1,ne-1
          c1 = dlog(ea(ie1))*omega(ie1+1)/
     &         (dlog(ea(ie1))-dlog(ea(ie1+1)))+
     &         dlog(ea(ie1+1))*omega(ie1)/
     &         (dlog(ea(ie1+1))-dlog(ea(ie1)))
          c0 = omega(ie1)/(dlog(ea(ie1))-dlog(ea(ie1+1))) +
     &         omega(ie1+1)/(dlog(ea(ie1+1))-dlog(ea(ie1)))
          energy(1) = ea(ie1)
          omg(1) = omega(ie1)
          do 190 ie2 = 1,npts
            energy(ie2+1) = (ea(ie1+1)-ea(ie1))*ie2/npts + ea(ie1)
            omg(ie2+1) = c0*dlog(energy(ie2+1)) + c1
            do 200 it = 1,ntemp
              ex = 3.d0/2.d0*te(it)
              sumi = (energy(ie2+1)-energy(ie2))/2.d0*
     &               ( omg(ie2)*dexp(evt/te(it)-
     &                (energy(ie2)/ex*dexp(lngama(5.d0/(2.d0*
     &                dru_val))-lngama(3.d0/(2.d0*dru_val))))**dru_val)
     &               + omg(ie2+1)*dexp(evt/te(it)-
     &                (energy(ie2+1)/ex*dexp(lngama(5.d0/(2.d0*
     &                dru_val))-lngama(3.d0/(2.d0*dru_val))))**dru_val))
              sum(it) = sum(it) + sumi*fine*c*a0*a0/wupper*pi*dsqrt(ryd)
     &                  *ex**(-1.5d0)*dru_val*dexp(1.5d0*
     &                  lngama(5.d0/(2.d0*dru_val))
     &                  -2.5d0*lngama(3.d0/(2.d0*dru_val)))
  200       enddo
  190     enddo
  180   enddo

C-----------------------------------------------------------------------
C  Add contributions above the last tabulated omega value. Continue this
C  until contribution is less than tolerance (tol), or until a maximum
C  number of points (limit) is reached.
C-----------------------------------------------------------------------

        energy(1) = ea(ne)
        omg(1) = omega(ne)
        estep = (ea(ne)-ea(ne-1))/npts
        do 210 it = 1,ntemp
          ex = 3.d0/2.d0*te(it)
          do 220 ie = 1,limit
            energy(ie+1) = energy(ie)+estep
            omg(ie+1) = c0*dlog(energy(ie+1)) + c1
            sumi = (energy(ie+1)-energy(ie))/2.d0*
     &             ( omg(ie)*dexp(evt/te(it)-
     &              (energy(ie)/ex*dexp(lngama(5.d0/(2.d0*
     &              dru_val))-lngama(3.d0/(2.d0*dru_val))))**dru_val)
     &             + omg(ie+1)*dexp(evt/te(it)-
     &              (energy(ie+1)/ex*dexp(lngama(5.d0/(2.d0*
     &              dru_val))-lngama(3.d0/(2.d0*dru_val))))**dru_val))
            sum(it) = sum(it) + sumi*fine*c*a0*a0/wupper*pi*dsqrt(ryd)
     &                *ex**(-1.5d0)*dru_val*dexp(1.5d0*
     &                lngama(5.d0/(2.d0*dru_val))
     &                -2.5d0*lngama(3.d0/(2.d0*dru_val)))
            if (sumi/sum(it).lt.tol) goto 230
  220     enddo
  230     continue
  210   enddo

      endif

C-----------------------------------------------------------------------
C  three-body recombination integral
C-----------------------------------------------------------------------
C  numerical distributions
C-----------------------------------------------------------------------

      if (dist.eq.2) then

        do 240 it = 1,ntemp

C-----------------------------------------------------------------------
C  Add distribution function points between ionisation potential and
C  first tabulation. Use low-energy behaviour. f(en) -> fn(enrg)
C-----------------------------------------------------------------------
		  temp3=3
          if (en(it,1).gt.1d-30.and.nform1.eq.2) then
			temp3=2
            estep = (en(it,1)-1d-30)/limit1
            do 450 ie = 1,limit1
              enrg(it,ie) = 1d-30+estep*(ie-1)
              fn(it,ie) = f(it,1)*(enrg(it,ie)/en(it,1))**param1
  450       enddo
            limit2(it) = limit1
          else
			temp3=1
            limit2(it) = 0
          endif

C-----------------------------------------------------------------------
C  Increase number of points of distribution function tabulation by
C  limit using high energy behaviour. f(en) -> fn(enrg)
C-----------------------------------------------------------------------

          estep = en(it,nef)-en(it,nef-1)
          do 420 ie = 1,nef
            enrg(it,ie+limit2(it)) = en(it,ie)
            fn(it,ie+limit2(it)) = f(it,ie)
  420     enddo
          do 430 ie = 1,limit2(it)
            enrg(it,nef+ie+limit2(it)) = en(it,nef)+estep*ie
			temp1=enrg(it,nef+ie+limit2(it))
			temp2=nef+ie+limit2(it)
            if (nform2.eq.1) then
              fn(it,nef+ie+limit2(it)) = 0d0
            elseif(nform2.eq.2) then
              fn(it,nef+ie+limit2(it)) = f(it,nef)*(en(it,nef)/
     &                           enrg(it,nef+ie+limit2(it)))**param2(1)
            elseif(nform2.eq.3) then
              fn(it,nef+ie+limit2(it)) = f(it,nef)*dexp(param2(1)/
     &                 te(it)*(en(it,nef)-enrg(it,nef+ie+limit2(it))))
            endif
            if (fn(it,nef+ie+limit2(it)).eq.0d0) then
              limit2(it) = limit2(it)+ie-1
              goto 440
            endif
  430     enddo
C          limit2(it) = limit2(it)+limit1
  440     continue

          do 250 ie1 = 1,nef+limit2(it)
            dist1(it,ie1) = fn(it,ie1)/dsqrt(enrg(it,ie1))
            do 260 ie2 = 1,nef+limit2(it)

C-----------------------------------------------------------------------
C  if projectile electron energy less than first energy point then use
C  low energy behaviour of distribution function
C-----------------------------------------------------------------------

              if ((enrg(it,ie1)-enrg(it,ie2)-evt).lt.enrg(it,1)) then
                if (nform1.eq.1) dist2(it,ie1,ie2) = 0.d0
                if (nform1.eq.2) dist2(it,ie1,ie2) = fn(it,1)*
     &            ((enrg(it,ie1)-enrg(it,ie2)-evt)/enrg(it,1))**param1

C-----------------------------------------------------------------------
C  interpolate distribution function to projectile electron energies
C  assuming exponential behaviour
C-----------------------------------------------------------------------

              else
                do 270 ie3 = 1,nef-1+limit2(it)
                  if ((enrg(it,ie1)-enrg(it,ie2)-evt).ge.enrg(it,ie3)
     &                .and. (enrg(it,ie1)-enrg(it,ie2)-evt).lt.
     &                enrg(it,ie3+1)) then
                    if (fn(it,ie3)/dsqrt(enrg(it,ie3)).ge.
     &                  fn(it,ie3+1)/dsqrt(enrg(it,ie3+1))) then
                      d1 = dlog(fn(it,ie3)/fn(it,ie3+1)*
     &                     dsqrt(enrg(it,ie3+1)/
     &                     enrg(it,ie3)))/(enrg(it,ie3+1)-enrg(it,ie3))
                      d0 = fn(it,ie3)/dsqrt(enrg(it,ie3))*(fn(it,ie3)/
     &                     fn(it,ie3+1)*dsqrt(enrg(it,ie3+1)/
     &                     enrg(it,ie3)))**
     &                     (enrg(it,ie3)/(enrg(it,ie3+1)-enrg(it,ie3)))
                      dist2(it,ie1,ie2) =
     &                  d0*dexp(-d1*(enrg(it,ie1)-enrg(it,ie2)-evt))
                    else
                      d1 =(fn(it,ie3)*enrg(it,ie3+1)/dsqrt(enrg(it,ie3))
     &                     -fn(it,ie3+1)*enrg(it,ie3)/
     &                     dsqrt(enrg(it,ie3+1)))/
     &                     (enrg(it,ie3+1)-enrg(it,ie3))
                      d0 = (fn(it,ie3+1)/dsqrt(enrg(it,ie3+1))-
     &                     fn(it,ie3)/dsqrt(enrg(it,ie3)))/
     &                     (enrg(it,ie3+1)-enrg(it,ie3))
                      dist2(it,ie1,ie2) =
     &                  d0*(enrg(it,ie1)-enrg(it,ie2)-evt)+d1
                    endif
                  endif
  270           enddo
              endif
  260       enddo

C-----------------------------------------------------------------------
C  interpolate omega to dist1 energies assuming c0*ln(energy)+c1
C  behaviour
C-----------------------------------------------------------------------

            do 280 ie4 = 1,ne-1!+limit2(it)
              c1 = (dlog(ea(ie4))*omega(ie4+1)-dlog(ea(ie4+1))*
     &             omega(ie4))/(dlog(ea(ie4))-dlog(ea(ie4+1)))
              c0 = (omega(ie4)-omega(ie4+1))/
     &             (dlog(ea(ie4))-dlog(ea(ie4+1)))
              if((enrg(it,ie1).ge.ea(ie4).and.enrg(it,ie1).lt.ea(ie4+1))
     &          .or.(enrg(it,ie1).lt.ea(1).and.ie4.eq.1)
     &          .or.(enrg(it,ie1).ge.xa(ne).and.ie4.eq.ne-1))
     &          oma(it,ie1) = c0*dlog(enrg(it,ie1)) + c1
  280       enddo
  250     enddo
  240   enddo


C-----------------------------------------------------------------------
C  temperature loop
C-----------------------------------------------------------------------

        do 290 it = 1,ntemp

C-----------------------------------------------------------------------
C  if within limits of integration then calculate integrand at each
C  incident/projectile energy combination
C-----------------------------------------------------------------------

          do 300 ie1 = 1,nef+limit2(it)
            do 310 ie2 = 1,nef+limit2(it)
              if ((enrg(it,ie2).ge.0.d0)             .and.
     &            (enrg(it,ie2).le.enrg(it,ie1)-evt) .and.
     &            (enrg(it,ie1).gt.evt)              )then
                int1(ie1,ie2) = dist1(it,ie2)*dist2(it,ie1,ie2)/
     &            ((enrg(it,ie2)+evt)**(2.d0)*
     &            (1.d0/evt-1.d0/enrg(it,ie1)))*oma(it,ie1)
              else
                int1(ie1,ie2) = 0.d0
              endif
  310       enddo

C-----------------------------------------------------------------------
C  sum the contributions from projectile energies
C-----------------------------------------------------------------------
            sum1(ie1) = 0.d0
            do 320 ie2 = 1,nef-1+limit2(it)
              sum1(ie1) = sum1(ie1) + (enrg(it,ie2+1)-enrg(it,ie2))/2d0*
     &                                (int1(ie1,ie2)+int1(ie1,ie2+1))
  320       enddo
  300     enddo

C-----------------------------------------------------------------------
C  sum the contributions from incident energies
C-----------------------------------------------------------------------

          sum2 = 0.d0
          do 330 ie1 = 1,nef-1+limit2(it)
            sum2 = sum2 + (enrg(it,ie1+1)-enrg(it,ie1))*
     &                    (sum1(ie1)+sum1(ie1+1))/2.d0
  330     enddo

C-----------------------------------------------------------------------
C  convert integrals to rate coefficients
C-----------------------------------------------------------------------

          alpha(it) = sum2*te(it)
     &                *pi*fine*c*a0*a0/wupper*dsqrt(pi*ryd*te(it))/2.d0
          q(it) = dexp(evt/te(it))*sum(it)

  290   enddo

C-----------------------------------------------------------------------
C  kappa and Druyvesteyn distributions
C-----------------------------------------------------------------------

      elseif(dist.ne.0) then

C-----------------------------------------------------------------------
C  Increase number of points of collision strength tabulation by factor
C  npts using linear interpolation. omega(ea) -> omg(energy)
C-----------------------------------------------------------------------

        energy(1) = ea(1)
        omg(1) = omega(1)
        do 340 ie1 = 1,ne-1
          c1 = dlog(ea(ie1))*omega(ie1+1)/
     &         (dlog(ea(ie1))-dlog(ea(ie1+1)))+
     &         dlog(ea(ie1+1))*omega(ie1)/
     &         (dlog(ea(ie1+1))-dlog(ea(ie1)))
          c0 = omega(ie1)/(dlog(ea(ie1))-dlog(ea(ie1+1))) +
     &         omega(ie1+1)/(dlog(ea(ie1+1))-dlog(ea(ie1)))
          do 350 ie2 = 1,npts
            energy(npts*(ie1-1)+ie2+1) = ea(ie1) +
     &                                   (ea(ie1+1)-ea(ie1))*ie2/npts
            omg(npts*(ie1-1)+ie2+1) = dlog(energy(npts*(ie1-1)+ie2+1))
     &                                *c0 + c1
  350     enddo
  340   enddo

C-----------------------------------------------------------------------
C  add points beyond the last tabulated omega by extrapolation
C-----------------------------------------------------------------------

        do ie = 1,limit
          energy(npts*(ne-1)+1+ie) = ea(ne)+estep*(ie-1)
          omg(npts*(ne-1)+1+ie) = dlog(energy(npts*(ne-1)+1+ie))*c0 + c1
        enddo

C-----------------------------------------------------------------------
C  calculate 3-body recombination integral
C-----------------------------------------------------------------------

        do 370 it = 1,ntemp
          kek = te(it)*(kap_val-1.5d0)
          ex = 3.d0/2.d0*te(it)

C-----------------------------------------------------------------------
C  Calculate integrand for every combination of impacting electron
C  energy (energy) and projectile electron energy (energy2).
C  Integration limits: evt < energy  < infinity
C                        0 < energy2 < energy-evt
C-----------------------------------------------------------------------

          do 380 ie1 = 1,(ne-1)*npts+1+limit
            do 390 ie2 = 1,(ne-1)*npts+1+limit
              energy2 = energy(ie2) - evt
              if ((energy2.ge.0.d0)           .and.
     &            (energy2.le.energy(ie1)-evt).and.
     &            (energy(ie1).gt.evt)       )then
                if (dist.eq.1) then
                  int1(ie1,ie2) = ((1.d0+(energy(ie1)-energy2-evt)
     &              /kek)*(1.d0+energy2/kek))**(-kap_val-1.d0)*
     &              omg(ie1)/((energy2+evt)**(2.d0)*
     &              (1.d0/evt-1.d0/energy(ie1)))
                elseif (dist.eq.3) then
                  int1(ie1,ie2) = dexp(-(energy2**dru_val+
     &              (energy(ie1)-energy2-evt)**dru_val)*
     &              (dexp(lngama(5.d0/(2.d0*dru_val))-
     &              lngama(3.d0/(2.d0*dru_val)))/ex)**dru_val)
     &              *omg(ie1)/
     &              ((energy2+evt)**(2.d0)*(1.d0/evt-1.d0/energy(ie1)))
                endif
              else
                int1(ie1,ie2) = 0.d0
              endif
  390       enddo

C-----------------------------------------------------------------------
C  for each incident electron energy, sum the contributions from the
C  projectile electron energies
C-----------------------------------------------------------------------

            sum1(ie1) = 0.d0
            do 400 ie2 = 1,(ne-1)*npts+limit
              sum1(ie1) = sum1(ie1) + (energy(ie2+1)-energy(ie2))/2.d0*
     &                                (int1(ie1,ie2)+int1(ie1,ie2+1))
  400       enddo
  380     enddo

C-----------------------------------------------------------------------
C  now sum the contributions from incident electron energies
C-----------------------------------------------------------------------

          sum2 = 0.d0
          do ie1 = 1,(ne-1)*npts+limit
            sum2 = sum2 + (energy(ie1+1)-energy(ie1))*
     &                    (sum1(ie1)+sum1(ie1+1))/2.d0
          enddo

C-----------------------------------------------------------------------
C  convert integral into 3-body recombination rate coefficient by
C  multiplying by relevant factor
C-----------------------------------------------------------------------

          if (dist.eq.1) then
            alpha(it) = sum2*2.d0*fine*c*dsqrt(ryd*pi)/wupper*a0*a0/
     &                  kek**(3.d0)*te(it)**(1.5d0)*dexp((2.d0*
     &                  (lngama(kap_val+1.d0)-lngama(kap_val-.5d0))))
            q(it) = dexp(evt/te(it))*sum(it)
          elseif (dist.eq.3) then
            alpha(it) = sum2/2.d0*fine*c*dsqrt(ryd*pi)/wupper*a0*a0*pi*
     &                  dru_val**(2.d0)/ex**(3.d0)*te(it)**(1.5d0)*
     &                  dexp(3.d0*lngama(5.d0/(2.d0*dru_val))
     &                  -5.d0*lngama(3.d0/(2.d0*dru_val)))
            q(it) = sum(it)
          endif
  370   enddo

      endif

      do it = 1,ntemp
        if (alpha(it).lt.1d-30) alpha(it) = 1.d-30
        if (q(it).lt.1d-30) q(it) = 1.d-30
      enddo

C-----------------------------------------------------------------------
      return
C-----------------------------------------------------------------------
1000  format(1x,31('*'),'  h9qd3b error  ',31('*')//
     &       1x,a)
C-----------------------------------------------------------------------
      end
