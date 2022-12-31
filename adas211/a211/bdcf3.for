       subroutine bdcf3( f    , e    , n    , l    , z     ,
     &                   x0   , x1   , h
     &                 )                              
       implicit none
C-----------------------------------------------------------------------
C
C  ****************** fortran77 routine: bcdf3 *************************
C
C  Purpose:  Tabulates asymptotically decaying bound coulomb function 
C
C
C  Subroutine:
C
C  input : (r*8)  e      = energy (Ryd) : must be <0 for a bound state
C  input : (r*8)  n      = principal quantum number
C  input : (r*8)  l      = orbital anular momentum quantum number
C  input : (r*8)  z      = ion charge +1
C  input : (r*8)  x0     = inner turning point of potential
C  input : (r*8)  x1     = outer turning point of potential
C  input : (r*8)  h      = interval for tabulation
C
C  output: (r*8)  f()    = Coulomb function
C                          1st dim: tabulation index
C 
C  Routines:
C          none
C
C  Author:  William Osborn (Tessella support services plc)
C
C  Date:     4 July 1996 
C
C  Update: MG O'Mullane  19/12/01  Removed junk from > column 72
C
C  Update: HP Summers    21/05/04  made implicit none and restructured  
C
C
C  Version  : 1.1                          Date: 04-07-96
C  Modified : William Osborn
C                - First version.
C
C  Version  : 1.2                          Date: 19-12-01
C  Modified : Martin O'Mullane
C                - removed junk from > column 72.
C
C  Version  : 1.3                          Date: 21-15-04
C  Modified : Hugh Summers
C                - Restructured as above.

C-----------------------------------------------------------------------
       integer    n     , l
       integer    j     , j0    , j1    , k     , k0
C-----------------------------------------------------------------------
       real*8     gama6
       real*8     e     , z     , x0    , x1    , h
       real*8     el    , gnu   , t     , t1    , t2    , t3    ,
     &            a     , b     , c     , x     , r     , s     ,
     &            th    , tk   
C-----------------------------------------------------------------------
       real*8     f(1000)                                                
C-----------------------------------------------------------------------
       el=l                                                             
       gnu=-z/dsqrt(-e)                                                 
       t1=gnu+el+1.0d0                                                    
       t2=gnu-el
                                                               
       if(t2.le.0.5d0)then                                                  
           t2=1.0
       endif
                                                                  
       t3=gama6(t1)*gama6(t2)                                           
       b=dsqrt(-z/t3)/gnu                                               
       b=b*(-1.0d0)**(n-l+1)                                              
       x=x0                                                             
       j0=0.5d0+x0/h                                                      
       j1=0.5d0+x1/h 
                                                            
       do j=j0,j1                                                     
         r=-z*x                                                           
         t=2.0d0*r/gnu                                                      
         t3=1.0d0/t                                                         
         s=1.0d0                                                            
         a=1.0d0                                                            
         k0=gnu+gnu+1.5d0+t                                                 
         th=k0                                                            
         th=gnu+gnu+1.0d0+t-th
                                                       
         do k=1,k0                                                      
           tk=k                                                             
           a=a*(tk-gnu+el)*(gnu-tk+el+1.0d0)*t3/tk                            
           t1=dabs (a)                                                      
           t2=dabs (s)                                                      
           if(t2.gt.t1*1.0d7)go to 2                                             
           s=s+a
         enddo
                                                                     
         t1=(gnu-el)*(gnu+el+1.0d0)                                         
         t2=-2.0d0*gnu*t1                                                   
         c=-0.5d0+0.125d0*t3*((2.0d0*th-1.0d0)+
     &      (th*th-1.5d0*th+0.25d0-2.0d0*t1)*t3)     
         s=s+c*a                                                          
    2    f(j)=(t**gnu)*s*dexp(-(0.5d0*t))*b                                 
         x=x+h
         
       enddo                                                            
       return                                                           
      end                                                               
