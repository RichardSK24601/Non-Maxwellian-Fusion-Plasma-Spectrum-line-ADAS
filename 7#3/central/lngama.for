      function lngama(x)
	
      implicit none
C-----------------------------------------------------------------------
C
C    ****************** fortran77 function: lngama ****************** 
C
C  purpose: Returns the natural logarithm of the Gamma function of x
C
C  calling program: various
C
C  input : (r*8)  x      = function argument
C
C  output: (r*8)  lngama = function name
C
C  author: Paul Bryans, University of Strathclyde
C
C  date:   26/01/04
C
C  update: 
C
C-----------------------------------------------------------------------
      real*8  lngama , ser , stp , tmp , x , y , z , cof(6)
C-----------------------------------------------------------------------
      integer j
C-----------------------------------------------------------------------      
      save    cof  , stp
C-----------------------------------------------------------------------
      data    cof  , stp /76.18009172947146d0  , -86.50532032941677d-0 ,
     &                    24.01409824083091d0  , -1.231739572450155d0  ,
     &                    .1208650973866179d-2 , -.5395239384953d-5    ,
     &                    2.5066282746310005d0/
C-----------------------------------------------------------------------
C
        y=x
	z=y
	tmp=y+5.5d0
	tmp=(y+0.5d0)*log(tmp)-tmp
	ser=1.000000000190015d0
	
	do j=1,6
	  z=z+1.d0
	  ser=ser+cof(j)/z
        enddo 
	
	lngama=tmp+log(stp*ser/y)
	
	return
	end
   
