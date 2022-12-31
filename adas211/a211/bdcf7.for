       subroutine bdcf7(n,l,qd,jsn,z0,nshell,nc,numel,alfa,zl,z1,z2,z3,
     &                  zs,x0,x1,x2,xmax,h,f,c)
       implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 program: bcdf7.for *********************
C
C  Purpose: Calculates a numerical radial wave function in a distorted
C           Coulomb potential described by shell screening.
C
C           The code permits a search for screening parameters given
C           the eigenenergy of for the energy given the screening
C           parameters
C           (original by A. Burgess, DAMTP, University of Cambridge)
C
C
C  Subroutine:
C
C  input : (i*4)  n        = principal quantum number
C  input : (i*4)  l        = orbital quantum number
C  input : (i*4)  qd       = quantum defect for valence electron
C  input : (i*4)  jsn      = -1 => Jucys potential form adopted
C                          = 0  => Slater potential form adopted
C  input : (i*4)  z0       = nuclear charge
C  input : (i*4)  nshell   = number of screening shells
C  input : (i*4)  nc()     = principal quantum number of screening shell
C                            1st dim: index of screening shells
C  input : (i*4)  numel()  = number of electrons in screening shell
C  i/o   : (r*8)  alfa()   = screening parameters
C                            1st dim: screening shell index.
C  input : (r*8)  zl()     =
C  input : (r*8)  z1       =
C  input : (r*8)  z2       =
C  input : (r*8)  z3       =
C  input : (r*8)  zs()     =
C  input : (r*8)  x0       =
C  input : (r*8)  x1       =
C  input : (r*8)  x2       =
C  input : (r*8)  xmax     =
C  input : (r*8)  h        =
C
C  output: (r*8)  f()      = wave function tabulation
C  output: (r*8)  c        =
C
C
C  Routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          zeffl      adas
C          bdcf3      adas
C          bdcf4      adas
C          fcf6       adas
C          i4unit     adas      fetch unit number for output of messages
C
C
C  Author:  William Osborn (Tessella support services plc)
C
C  Date:    4th july 1996
C
C  Update: M O'Mullane 19-12-01 removed junk from > column 72.
C
C
C  Update: HP Summers  24/05/04  restructure and addded standard warning
C
C
C  Unix-idl port:
C
C  Version: 1.1                          Date: 04-07-96
C  Modified: William Osborn
C               - first version.
C
C  Version: 1.2                          Date: 19-12-01
C  Modified: Martin O'Mullane
C               - removed junk from > column 72.
C
C  Version: 1.3                          Date: 25-05-2004
C  Modified: H P Summers
C               - restructure.
C
C  Version: 1.4                          Date: 17-05-2007
C  Modified: Allan Whiteford
C               - Updated comments as part of subroutine documentation
C                 procedure.
C
C-----------------------------------------------------------------------
       integer  i4unit
       integer  n       , l        , jsn    , nshell
       integer  i       , i1       , i2     , imax    ,
     &          i11     , j        , jh
C-----------------------------------------------------------------------
       integer  nc(10),numel(10)
C-----------------------------------------------------------------------
       real*8   zeffl   , bdcf4
       real*8   qd      , z0       , z1     , z2     , z3     ,
     &          x0      , x1       , x2     , xmax   , h      ,
     &          c
       real*8   za      , t        , en     , el     , w2     ,
     &          h2      , hh       , h1     , w1     , c0     ,
     &          c1      , c2       , a1     , a2     , a3     ,
     &          del     , e        , fi1    , phi    ,  r1    ,
     &          r2      , t3       , w0     , x
C-----------------------------------------------------------------------
       real*8   f(1000) , zl(1000) , zs(100)   , amp(20)   , alfa(10)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
       za=z0+z1
       if(nshell)22,22,20
   20  do 21 j=1,nshell
       t=numel(j)
       if(j-jsn)24,23,24
   23  t=numel(j)-1
   24  za=za+t
   21  continue
   22  continue
       i1=x1/h+0.5d0
       x1=i1
       x1=x1*h
       imax=xmax/h+0.5d0
       en=n
       t=za/(en-qd)
       e=-t*t
       t=x1
       call fcf6(f,c,a2,amp,phi,del,x0,n,l,e,jsn,z0,nshell,nc,numel,
     1 alfa,zl,z1,z2,z3,zs,x1,h,x2,h2)
       x1=t
       i2=x2/h+0.5d0
       x2=i2
       x2=x2*h
       if(imax-i1)1,2,2
    1  write(i4unit(-1),2000)'failed in bdcf7, xmax too small'
       write(i4unit(-1),2001)
       return
    2  el=l
       w2=-el*(el+1.0d0)
       jh=h/h2+0.5d0
       hh=h2*h2
       h1=hh/12.0d0
       x=x2+h+h2
       t=1.0d0/x
       w1=e-2.0d0*za*t+w2*t*t
       r1=bdcf4(e,n,l,za,x)
       x=x2+h
       t=1.0d0/x
       w0=e-2.0d0*za*t+w2*t*t
       r2=bdcf4(e,n,l,za,x)
       c0=r2*(1.0d0+(h1-hh)*w0)-r1*(1.0d0+h1*w1)
       c1=r2*(1.0d0+h1*w0)
       fi1=f(i1)
       i=i2+1
    3  i=i-1
       do 4 j=1,jh
       x=x-h2
       c1=c1+c0
       t=1.0d0/x
       t3=2.0d0*zeffl(jsn,z0,nshell,nc,numel,alfa,x,zl,h,x1,z1,z2,z3)/x
       c2=e-t3+w2*t*t
       r2=c1/(1.0d0+h1*c2)
       c0=c0-hh*c2*r2
    4  continue
       if(i-imax)5,5,6
    5  f(i)=r2
    6  if(i-i1)7,7,3
    7  c1=f(i1)/fi1
       i11=i1-1
       do 8 i=1,i11
       f(i)=c1*f(i)
    8  continue
       if(i2-imax)9,10,10
    9  x=x2+h
       call bdcf3(f,e,n,l,za,x,xmax,h)
   10  t=0.0d0
       a3=0.0d0
       i=0
   11  i=i+2
       a1=a3
       a2=f(i-1)*f(i-1)
       a3=f(i)*f(i)
       t=t+a1+4.0d0*a2+a3
       if(i-imax)11,12,12
   12  c0=1.0d0/dsqrt(-e)
       a1=1.0d0-2.0d0*za*c0
       t3=0.5d0*c0/xmax
       a2=a3
   15  a1=a1-1.0d0
       if(a1)17,17,16
   16  a2=a2*a1*t3
       a3=a3+a2
       go to 15
   17  t=t*h/3.0d0+0.5d0*a3*c0
       c2=1.0d0/dsqrt(t)
       do 13 i=1,imax
       f(i)=c2*f(i)
   13  continue
       c=c*c1*c2

       return

C-----------------------------------------------------------------------
 2000 format(1x,32('*'),' bdcf7 warning  ',32('*')//
     &       2x,a )
 2001 format(/1x,29('*'),' program continues ',29('*'))
C-----------------------------------------------------------------------
      end
