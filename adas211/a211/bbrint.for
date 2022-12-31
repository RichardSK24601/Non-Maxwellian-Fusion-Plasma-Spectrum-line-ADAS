      subroutine bbrint( ndgnt  ,
     &		  	 gaunt  , vve    , z1     , v    ,
     &		         maxe   , temp   , dparam , dist , 
     &			 f      , rrcint
     &			)
      
      implicit none
C-----------------------------------------------------------------------
C
C    ****************** fortran77 subroutine: bbrint ****************** 
C
C  purpose: To calculate radiative recombination coefficient when 
C           electron distribution is not Maxwellian.
C
C  calling program: adas211
C
C  input : (i*4)  ndgnt   = max no of vve gaunt pairs
C  input : (r*8)  gaunt() = the bound-free gaunt factor
C  input : (r*8)  vve()   = v**2*e 
C                           where e=(free electron energy)/z1**2 (ryd)
C  input : (r*8)  z1      = parent ion charge
C  input : (r*8)  v       = effective principal quantum number
C  input : (i*4)  maxe    = number of vve gaunt pairs
C  input : (r*8)  temp    = effective temperature (kelvin),
C                           2/3 * mean energy of distribution
C  input : (r*8)  dparam  = parameter describing distribution function:
C                             kappa dist.       => kappa
C                             Druyvesteyn dist. => x
C  input : (i*4)  dist    = non-Maxwellian distribution type:
C                             1 => kappa distribution
C                             2 => numerical distribution
C                             3 => Druyvesteyn distribution
C  input : (r*8)  f()     = numerical distribution function at vve
C
C  output: (r*8)  rrcint  = radiative recombination coefficient (cm3 sec-1)
C
C  local : (r*8)  ryd     = Rydberg constant (eV)
C  local : (r*8)  te      = effective temperature (eV)
C  local : (r*8)  ip      = ionisation potential (eV)
C  local : (r*8)  kek     = kappa * characteristic energy of kappa dist.
C  local : (r*8)  ex      = characteristic energy of Druyvesteyn dist.
C  local : (r*8)  alpha   = fine structure constant
C  local : (r*8)  c       = speed of light in vacuum (cm sec-1)
C  local : (r*8)  a0      = Bohr radius (cm)
C  local : (r*8)  energy()= free electron energy (eV)
C  local : (r*8)  int1    = integrand at energy(i)
C  local : (r*8)  int2    = integrand at energy(i+1)
C  local : (r*8)  de      = energy difference from i to i+1
C
C  routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          lngama               evaluates ln(gamma(x))
C
C  author: Paul Bryans, University of Strathclyde
C
C  date:   23/01/04
C
C  update: 26/01/04 - Paul Bryans
C          added Druyvesteyn distribution (dist = 3)
C
C  update: 02/12/04 - Paul Bryans
C          added numerical distribution (dist = 2)
C
C  update: 02/02/05 - Allan Whiteford
C          Declared i4unit as an integer.
C
C  update: 20/07/07 - Allan Whiteford
C          Removed commet stating that Druyvesteyn and numerical
C          distributions can't be handled.
C
C-----------------------------------------------------------------------
      integer     maxe          , dist   , i   , ndgnt , i4unit
C-----------------------------------------------------------------------
      real*8      gaunt(ndgnt)  , v      , pi  , temp  , rrcint ,
     &            vve(ndgnt)    , te     , ip  , int1  , dparam ,
     &            energy(ndgnt) , c      , a0  , kek   , int2   ,     
     &            sum(ndgnt)    , de     , z1  , ryd   , alpha  ,
     &            f(ndgnt)      , lngama , ex
C-----------------------------------------------------------------------
      parameter ( ryd = 13.606d0  , alpha = 7.2974d-3 ,
     &            c   = 2.9979d10 , a0    = 5.2918d-9  )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  set values of physical constants
C-----------------------------------------------------------------------
C
      pi=acos(-1.0d0)      
      te=temp/11605.4d0
      ip=(z1/v)**2d0*ryd
      energy(1)=vve(1)*ip
      kek=(dparam-1.5d0)*te
      ex=1.5d0*te
C
C-----------------------------------------------------------------------
C  initialise rrcint to zero
C-----------------------------------------------------------------------
C      
      rrcint=0.0d0
C
C-----------------------------------------------------------------------
C  loop over energy points and add to integral after each iteration
C-----------------------------------------------------------------------
C       
      do i=1,maxe-1
	energy(i+1)=vve(i+1)*ip
	de=energy(i+1)-energy(i)
	if (dist.eq.1) then
	  int1=gaunt(i)/(energy(i)+ip)*(1.0d0+energy(i)/kek)
     &	        **(-dparam-1.0d0)
	  int2=gaunt(i+1)/(energy(i+1)+ip)*(1.0d0+energy(i+1)/kek)
     &	        **(-dparam-1.0d0)
	  sum(i)=de/2.0d0*(int1+int2)*
     &           2.0d0/dsqrt(pi)*(ryd/kek)**(1.5d0)*
     &           dexp(lngama(dparam+1)-lngama(dparam-.5d0))
	elseif (dist.eq.2) then
	  int1 = gaunt(i)/(energy(i)+ip)*f(i)
	  int2 = gaunt(i+1)/(energy(i+1)+ip)*f(i+1)
	  sum(i) = de/2d0*(int1+int2)*ryd**1.5d0
	elseif (dist.eq.3) then
	  int1=gaunt(i)/(energy(i)+ip)*dexp(-(energy(i)/ex*dexp
     &	       (lngama(5.0d0/(2.0d0*dparam))-
     &          lngama(3.0d0/(2.0d0*dparam))))**dparam)
	  int2=gaunt(i+1)/(energy(i+1)+ip)*dexp(-(energy(i+1)/ex*dexp
     &	       (lngama(5.0d0/(2.0d0*dparam))-
     &          lngama(3.0d0/(2.0d0*dparam))))**dparam)  
   	  sum(i)=de/2.0d0*(int1+int2)*
     &           dparam*(ryd/ex)**(1.5d0)*dexp(3.0d0/2.0d0*
     &           lngama(5.0d0/(2.0d0*dparam))-5.0d0/2.0d0*
     &           lngama(3.0d0/(2.0d0*dparam)))
	else
          write(i4unit(-1),1000)'(dist .ne. 1, 2 or 3) -'
	  write(i4unit(-1),1001)
     &    'invalid electron distribution type ',
     &    '(must equal 1, 2 or 3)'
	  write(i4unit(-1),1002)
	endif
	rrcint=rrcint+sum(i)
      enddo
C
C-----------------------------------------------------------------------
C  convert integral to radiative recombination coefficient
C-----------------------------------------------------------------------
C      
      rrcint=rrcint*32.0d0*alpha**4d0*c*pi*a0**2d0*z1**4d0/
     &       (3.0d0**(1.50d0)*v**3d0)
C
C-----------------------------------------------------------------------
C      
 1000 format(1x,32('*'),' bbrint error ',32('*')//
     &       1x,'Fault in distribution type: ',a)
 1001 format(1x,a) 
 1002 format(/1x,29('*'),' program terminated ',29('*'))
C 
C-----------------------------------------------------------------------
      return
      end
      
