      subroutine bbspln( ndtem   , ntmax   ,
     &                   nblock  , maxt    ,
     &                   tin     , tout    , 
     &                   rrcin   , rrcout
     &                 )
      implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 subroutine: bbspln *********************
C
C  purpose:
C 	  1) performs cubic spline on log(temp) versus log(rad.rec.coeff.)
C 	     input data. ('tin' versus 'rrcin' , nblock data pairs)
C
C 	  2) interpolates 'maxt'  rrcout values using above splines at
C 	     temperatures read in from adf08 file for tabular output.

C
C  calling program: adas211
C
C
C  subroutine:
C
C
C  input : (i*4)  ndtem   = maximum number of adf08 temperatures
C  input : (i*4)  ntmax   = maximum number of adf37 temperatures
C  input : (i*4)  nblock  = input data file: number of rrc/temperature
C 			    pairs read for the transition being assessed
C  input : (i*4)  maxt    = number of adf08 temperature values at
C 			    which interpolated rrc values are required
C 			    for tabular output.
C  input : (r*8)  tin()   = adf37 temperatures (kelvin)
C  input : (r*8)  tout()  = adf08 entered temperatures (kelvin)
C  input : (r*8)  rrcin() = rrc values at 'tin()'.
C
C  output: (r*8)  rrcout()= spline interpolated rrc values at 'tout()'
C
C  local : (i*4)  nin	  = parameter = max. no. of input temp/rrc pairs
C 				        must be >= 'nv'
C  local : (i*4)  nout    = parameter = max. no. of output temp/rrc pairs
C 				        must be >= 'maxt' & 'npspl'
C  local : (i*4)  iarr    = array subscript used for temp/rrc pairs
C  local : (i*4)  iopt    = defines the boundary derivatives for the
C 			    spline routine 'xxspln', see 'xxspln'.
C 			    (valid values = <0, 0, 1, 2, 3, 4).
C  local : (l*4)  lsetx   = .true.  => set up spline parameters relating
C 				       to 'xin' axis.
C 			    .false. => do not set up spline parameters
C 				       relating to 'xin' axis.
C 				       (they were set in a previous call)
C 			    (value set to .false. by 'xxsple')
C  local : (r*8)  xin()   = log( 'tin()' )
C  local : (r*8)  yin()   = log( 'rrcin()' )
C  local : (r*8)  xout()  = log(temperatures at which splines required)
C  local : (r*8)  yout()  = log(output spline interpolated rrc values)
C  local : (r*8)  df()    = spline interpolated derivatives
C
C
C routines:
C 	   routine    source	brief description
C 	   ------------------------------------------------------------
C 	   xxspln     adas	spline subroutine
C 	   r8fun1     adas	real*8 function: ( x -> x )
C
C author:  Paul Bryans (University of Strathclyde)
C
C date:    01/12/04
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      integer    nin           , nout        , nblock  , ndtem
      integer    maxt          , iarr        , iopt    , ntmax
C-----------------------------------------------------------------------
      parameter( nin = 51      , nout = 100  )
C-----------------------------------------------------------------------
      logical    lsetx
C-----------------------------------------------------------------------
      real*8     rrcin(ntmax)  , tin(ntmax)  , r8fun1
      real*8     rrcout(ndtem) , tout(ndtem) , df(nin)
      real*8     xin(nin)      , xout(nout)
      real*8     yin(nin)      , yout(nout)
C-----------------------------------------------------------------------
      external   r8fun1
C-----------------------------------------------------------------------
      if (nin.lt.nblock)
     &             stop ' bbspln error: nin < nblock - increase nin'
      if (nout.lt.maxt)
     &             stop ' bbspln error: nout < maxt - increase ntout'
C
C-----------------------------------------------------------------------
C  set up spline boundary conditions - switch on extrapolation.
C-----------------------------------------------------------------------
C
      iopt  = 0
C
C-----------------------------------------------------------------------
C  interpolate splined temp/rrc pairs
C  converts from adf37 temperatures to adf08 temperatures.
C-----------------------------------------------------------------------
C
      do 1 iarr = 1,nblock
	xin(iarr) = dlog(tin(iarr))
    1	yin(iarr) = dlog(rrcin(iarr))

      do 2 iarr = 1,maxt
    2	xout(iarr) = dlog(tout(iarr))

      lsetx = .true.
      call xxspln( lsetx  , iopt   , r8fun1 ,
     &             nblock , xin    , yin    ,
     &   	   maxt   , xout   , yout   ,
     &  	   df
     &  	 )

      do 3 iarr=1,maxt
    3	rrcout(iarr) = dexp(yout(iarr))
C
C-----------------------------------------------------------------------
C
      return
      end


