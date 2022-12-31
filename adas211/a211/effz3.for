CX UNIX PORT - SCCS info: Module @(#)effz3.for	1.2 Date 03/15/02
CX
       subroutine effz3( jealfa , n     , l      , e     , qd    ,
     &                   jsn    , z0    , nshell , nc    , numel ,
     &                   alfa   , jalf1 , jalf2  ,
     &                   x0     , x1    , x2     , d     , m0
     &                 )
      implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 program: effz3.for *********************
C
C  Purpose: Searches for the effective potential for a single electron
C           distorted wave function for a specified screening or a
C           specified energy.
C           (original by A. Burgess, DAMTP, University of Cambridge)
C
C
C  Subroutine:
C
C  input : (i*4)  jealfa   = <0 => search for energy e
C                          = >0 => search for screening parameter alfa
C  input : (i*4)  n        = principal quantum number
C  input : (i*4)  l        = orbital quantum number
C  i/o   : (i*4)  e        = energy (Ryd) for electron.
C                             (NB -ve for a bound state)
C  input : (i*4)  qd       = quantum defect for valence electron
C  input : (i*4)  jsn      = -1 => Jucys potential form adopted
C                          = 0  => Slater potential form adopted
C  input : (i*4)  z0       = nuclear charge
C  input : (i*4)  nshell   = number of screening shells
C  input : (i*4)  nc()     = principal quantum number of screening shell
C                            1st dim: index of screening shells
C  input : (i*4)  numel()  = number of electrons in screening shell
C  i/o   : (r*8)  alfa()   = screening parameters
C                            1st dim: initial (1) and final (2) states
C                            2nd dim: screening shell index.
C  input : (i*4)  jalf1    = first screening shell for optimising
C  input : (i*4)  jalf2    = last screening shell for optimising
C  output: (i*4)  x0       = inner turning point
C  output: (i*4)  x1       = outer turning point
C  output: (i*4)  x2       = range for active electron wave function
C  input : (i*4)  d        = earch accuracy setting
C  output: (i*4)  m0       = number of nodes in wave function
C
C
C
C  Routines:
C          routine    source    brief description
C          -------------------------------------------------------------
C          zeff       adas
C          zser       adas
C          fcf6       adas
C          i4unit     adas      fetch unit number for output of messages
C
C
C  Author:  H. P. Summers, University of Strathclyde
C           ja7.08
C           tel. 0141-548-4196
C
C  Date:   24/02/03
C
C  Update: HP Summers  24/05/04  restructure and addded standard warning
C
C-----------------------------------------------------------------------
       integer   jmax      , i1
C-----------------------------------------------------------------------
       real*8    s1
C-----------------------------------------------------------------------
       parameter ( jmax = 60 , i1 = 256 )
       parameter ( s1 = 1.2d0 )
C-----------------------------------------------------------------------
       integer   i4unit
       integer   jealfa    , n      , l      , jalf1   , jalf2    ,
     &           jsn       , nshell , m0
       integer   m         , m1     , ja     , j       , i
C-----------------------------------------------------------------------
       integer   nc(10)    , numel(10)
C-----------------------------------------------------------------------
       real*8    zeff
       real*8    e         , qd     , z0     ,
     &           x0        , x1     , x2     , d
       real*8    s2        , ti1    , h      , ell     , alf      ,
     &           z1        , z2     , z3     , t1      , t        ,
     &           xa        , xb     , ya     , yb      , del      ,
     &           ymin      , a2     , c      , h2      , phi      ,
     &           qd1       , x
C-----------------------------------------------------------------------
       real*8    f(1000)   , zl(1000)    , zs(100)     , amp(20)
       real*8    alfa(10)  , tal(10)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
       s2=1.0d0/s1
       ti1=i1
       h=x1/ti1
       ell=l
       ell=ell*(ell+1.0d0)
       if(e+1.0d-40)1,3,3
    1  if(jealfa)2,2,3
    2  alf=-e
       go to 4
    3  alf=1.0d0
       do 61 ja=jalf1,jalf2
       tal(ja)=alfa(ja)
   61  continue
    4  continue
       m=-1
       m1=1
       z1=0.0d0
       z2=0.0d0
       z3=0.0d0
       j=0
    5  j=j+1
       if(j-jmax)7,7,6
    6  write(i4unit(-1),2000)'Failed in effz3'
       write(i4unit(-1),2001)
       go to 46
    7  continue
       if(e+1.0d-40)8,14,14
    8  if(j-1)15,15,9
    9  t1=zeff(jsn,z0,nshell,nc,numel,alfa,x1)
       t=t1*t1+ell*e
       if(t)10,12,11
   10  t=0.0d0
       go to 12
   11  t=dsqrt(t)
   12  x1=0.5d0*(x1+(t1-t)/e)
       h=x1/ti1
       if(jealfa)13,13,14
   13  e=-alf
       go to 19
   14  do 60 ja=jalf1,jalf2
       alfa(ja)=alf*tal(ja)
   60  continue
   15  continue
       call zser(jsn,z0,nshell,nc,numel,alfa,zs)
   19  x=0.0d0
       do 20 i=1,i1
       x=x+h
       zl(i)=zeff(jsn,z0,nshell,nc,numel,alfa,x)
   20  continue
       call fcf6(f,c,a2,amp,phi,qd1,x0,n,l,e,jsn,z0,nshell,nc,numel,
     1 alfa,zl,z1,z2,z3,zs,x1,h,x2,h2)
       if(e+1.0d-40)21,22,22
   21  del=qd1
       go to 23
   22  del=qd-qd1
   23  continue
       if(m)24,27,37
   24  xa=alf
       ya=del
       if(del)25,46,26
   25  alf=s1*alf
       m=0
       go to 5
   26  alf=s2*alf
       m=0
       go to 5
   27  xb=alf
       yb=del
       if(ya*yb)31,46,28
   28  if(yb)29,46,30
   29  alf=s1*alf
       ya=yb
       xa=xb
       go to 5
   30  alf=s2*alf
       ya=yb
       xa=xb
       go to 5
   31  if(yb)32,32,33
   32  t=xa
       xa=xb
       xb=t
       t=ya
       ya=yb
       yb=t
   33  if(ya+yb)34,35,35
   34  ymin=yb
       go to 36
   35  ymin=-ya
   36  ymin=d*ymin
       alf=0.5d0*(xa+xb)
       m=1
       go to 5
   37  if(del)38,46,39
   38  t=xa
       xa=alf
       ya=del
       go to 40
   39  t=xb
       xb=alf
       yb=del
   40  if(dabs(alf-t)-d*alf)41,41,43
   41  if(dabs(del)-d)42,42,43
   42  alf=(xa*yb-xb*ya)/(yb-ya)
       go to 46
   43  m1=-m1
       if(m1)45,45,44
   44  alf=(xa*yb-xb*ya)/(yb-ya)
       go to 5
   45  alf=0.5d0*(xa+xb)
       go to 5
   46  m0=phi+0.5d0
       if(e+1.0d-40)47,49,49
   47  if(jealfa)48,48,49
   48  e=-alf
       go to 50
   49  do 62 ja=jalf1,jalf2
       alfa(ja)=alf*tal(ja)
   62  continue

   50  return

C-----------------------------------------------------------------------
 2000 format(1x,32('*'),' dfdip0 warning  ',32('*')//
     &       2x,a )
 2001 format(/1x,29('*'),' program continues ',29('*'))
C-----------------------------------------------------------------------
      end
