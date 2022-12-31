      function ingamq(a,x)
	
      implicit none
C-----------------------------------------------------------------------
C
C    ****************** fortran77 function: ingamq ****************** 
C
C  purpose: evaluates incomplete gamma function, 1-P(a,x)
C
C  calling program: various
C
C  input : (r*8)  x      = function argument
C  input : (r*8)  a      = function argument
C
C  output: (r*8)  ingamq = function name
C
C  author: Paul Bryans, University of Strathclyde
C
C  date:   20/02/04
C
C  update: 
C
C-----------------------------------------------------------------------
      real*8  ingamq , a      , x
      real*8  ap     , sum    , del    , tol
      real*8  lngama , b      , c      , d
      real*8  h      , an     , fpmin
C-----------------------------------------------------------------------
      integer itmax  , i      , i4unit
C-----------------------------------------------------------------------
      parameter ( itmax=100, tol=3.d-12, fpmin=1.d-30 )
C-----------------------------------------------------------------------
C
      if (x.lt.a+1.d0) then
      
        ap=a
	sum=1.d0/a
	del=sum
	do 10 i=1,itmax
	  ap=ap+1.d0
	  del=del*x/ap
	  sum=sum+del
	  if (abs(del).lt.abs(sum)*tol) goto 20
 10	enddo
	write(i4unit(-1),'(a)') 'tolerance is too low: increase itmax'
 20	ingamq=1d0-sum*dexp(-x+a*dlog(x)-lngama(a))
	
      else
      
        b=x+1.d0-a
	c=1.d0/fpmin
	d=1.d0/b
	h=d
	do 30 i=1,itmax
	  an=-i*(i-a)
	  b=b+2.d0
	  d=an*d+b
	  if (abs(d).lt.fpmin) d=fpmin
	  c=b+an/c
	  if (abs(c).lt.fpmin) c=fpmin
	  d=1.d0/d
	  del=d*c
	  h=h*del
	  if (abs(del-1.d0).lt.tol) goto 40
 30	enddo
	write(i4unit(-1),'(a)') 'tolerance is too low: increase itmax'
 40	ingamq=dexp(-x+a*dlog(x)-lngama(a))*h
	
      endif
      
      return
      end
