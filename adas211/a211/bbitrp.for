      subroutine bbitrp( ninmx  , noutmx ,
     &                   nein   , neout  , te     ,
     &		         nform1 , param1 , nform2 , param2 ,
     &			 ein    , fin    , eout   , fout
     &		        )
      
      implicit none
C-----------------------------------------------------------------------
C
C    ****************** fortran77 subroutine: bbitrp ****************** 
C
C  purpose: To interpolate/extrapolate numerical distribution from 
C           fin(ein) to fout(eout).
C           A f=sqrt(E)*exp(-E) fit is chosen for interpolation
C           Extrapolation uses limit behaviour from nform1 and nform2
C
C  calling program: adas211
C
C  input : (i*4)  ninmx    = max no of input energies
C  input : (i*4)  noutmx   = max no of output energies
C  input : (i*4)  nein     = no of input energies
C  input : (i*4)  neout    = no of output energies
C  input : (r*8)  te       = temperature
C  input : (i*4)  nform1   = type of threshold behaviour
C                              1 => cutoff
C                              2 => energy^param1
C  input : (r*8)  param1   = parameter of threshold form
C  input : (i*4)  nform2   = type of high-energy behaviour
C                              1 => cutoff
C                              2 => energy^-param2(1)
C                              3 => exp(-param2(1)*energy)
C                              4 => exp(-param2(1)*energy^param2(2))
C  input : (r*8)  param2() = parameter of high-energy form
C
C  input : (r*8)  ein()    = input energy of distribution
C  input : (r*8)  fin()    = value of distribution at ein
C  input:  (r*8)  eout()   = output energy
C
C  output: (r*8)  fout()   = (value of distribution at eout)/sqrt(eout)
C
C
C  author: Paul Bryans, University of Strathclyde
C
C  date:   30/11/04
C
C  update: 
C
C-----------------------------------------------------------------------
      integer i      , j         , neout        ,
     &        ninmx  , noutmx    , nein         , nform1       , nform2 
C-----------------------------------------------------------------------
      real*8  param1 , param2(2) , ein(ninmx)   , fin(ninmx)   ,
     &        a0     , a1        , eout(noutmx) , fout(noutmx) ,
     &        te
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C  determine interpolation parameters
C  f(E)=a1*sqrt(E)*exp(-a0*E)
C-----------------------------------------------------------------------
C
      do 1 i = 1,nein-1
	a0 = dlog(fin(i)/fin(i+1)*dsqrt(ein(i+1)/ein(i)))/
     &	     (ein(i+1)-ein(i))
        a1 = fin(i)/dsqrt(ein(i))*(fin(i)/fin(i+1)*
     &	     dsqrt(ein(i+1)/ein(i)))**(ein(i)/(ein(i+1)-ein(i)))
	do 2 j = 1,neout
	  if (eout(j).ge.ein(i).and.eout(j).lt.ein(i+1))
     &	    fout(j) = a1*dexp(-a0*eout(j))
    2	continue
    1 continue
C
C-----------------------------------------------------------------------
C  extrapolate if eout < lowest ein
C-----------------------------------------------------------------------
C    
      do 3 j = 1,neout
	if (eout(j).lt.ein(1)) then
	  if (nform1.eq.2) then
	    fout(j) = fin(1)*eout(j)**(param1-.5d0)*ein(1)**(-param1)
	  else
	    fout(j) = 0d0
	  endif
C
C-----------------------------------------------------------------------
C  extrapolate if eout > highest ein
C-----------------------------------------------------------------------
C 	
	elseif (eout(j).ge.ein(nein)) then
	  if (nform2.eq.2) then
	    fout(j) = fin(nein)*(ein(nein)/eout(j))**param2(1)
	  elseif (nform2.eq.3) then
	    fout(j) = fin(nein)*dexp(param2(1)/te*(ein(nein)-eout(j)))/
     &	              dsqrt(eout(j))
	  else
	    fout(j) = 0d0
	  endif
	endif
    3 continue	  
C-----------------------------------------------------------------------	
      return
      end	
            
