       subroutine remapx( ndtrn , nvmax , nv , itran , scef  , scom  )

       implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 subroutine: remapx *********************
C
C  PURPOSE:  Remaps threshold unit parameter X'=X+a/(b+X) where a=3, b=1
C            corresponds to Cowan's remapping, and corresponding Omegas.
C
C  calling program: main
C
C  subroutine:
C
C  input : (i*4)  ndtrn   = max. number of transitions that can be read.
C  input : (i*4)  nvmax   = max. number of temperatures that can be read.
C  i/oput: (i*4)  nv      = input data file: number of gamma/temperature
C                           pairs for a given transition.
C  i/oput: (r*8)  scef()  = input data file: electron temperatures (k)
C                           (note: te=tp=th is assumed)
C  input : (i*4)  itran   = input data file: number of transitions
C  i/oput: (r*8)  scom(,) = transition omegaa values 
C                           1st dimension - temperature 'scef()'
C                           2nd dimension - transition number
C
C AUTHOR:  N R Badnell, University of Strathclyde
C          badnell@phys.strath.ac.uk
C
C DATE:    27/10/04
C
C UPDATE:
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      integer   ndtrn , nvmax , nv , itran
      integer   i , j , n , nv0, iw, i0
C-----------------------------------------------------------------------
      real*8    scef(nvmax) , scom(nvmax,ndtrn)
      real*8    a , b, x, xp, x1, x2, w1, w2
C-----------------------------------------------------------------------
      parameter( a = 3.0d0 , b = 1.0d0 )
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C

C-----------------------------------------------------------------------
C Remap X-values.
C-----------------------------------------------------------------------

      iw=-1
      nv0=nv
      nv=0
      do n=1,nv0
        xp=scef(n)
        x=(xp-b)**2+4.0d0*(b*xp-a)
        if(x.ge.0.0d0)then
          x=(sqrt(x)+xp-b)*0.5d0
          if(x.ge.1.0d0)then
            if(nv.eq.0)then
              iw=n-1
              x1=scef(iw)
              x2=scef(iw+1)
            endif
            nv=nv+1
            scef(nv)=x
          endif
        endif
      enddo

C-----------------------------------------------------------------------
C Determine interpolate weights for X=1 (if necessary).
C-----------------------------------------------------------------------

      if(scef(1).gt.1.01d0)then
        i0=1
        do n=nv,1,-1
          scef(n+1)=scef(n)
        enddo
        scef(1)=1.0d0
        xp=1.0d0+a/(b+1.0d0)
        w1=(x2-xp)/(x2-x1)
        w2=(xp-x1)/(x2-x1)
      else
        i0=0
      endif

!      write(*,*)'Reducing number of X-values from ', nv0, ' to ', nv+i0
!      write(*,*) (n,scef(n),n=1,nv+i0)

C-----------------------------------------------------------------------
C Remap omegas,discarding nv0-nv X-values.
C-----------------------------------------------------------------------

      j=nv0-nv
      do i=1,itran
        if(i0.gt.0)scom(1,i)=w1*scom(iw,i)+w2*scom(iw+1,i)
!        write(*,*)i,scom(iw,i),scom(1,i),scom(iw+1,i)
        do n=1,nv
          scom(n+i0,i)=scom(n+j,i)
        enddo
      enddo
      nv=nv+i0

C-----------------------------------------------------------------------

      return
      end
