C UNIX-IDL - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/ass2.for,v 1.3 2004/07/06 11:17:11 whitefor Exp $ Date $Date: 2004/07/06 11:17:11 $
C
       subroutine ass2( x1    , h     , x    , f0    , f1    ,
     &                  g0    , g1    , ei   , ej    , tkij  ,
     &                  li    , lj    , z    , ni    , rem
     &                )    
       implicit none
c-----------------------------------------------------------------------
c
c  ****************** fortran77 subroutine: ass2   *********************
c
c  purpose:  calculate asymptotic part of integrals over wave functions
c
c
c  subroutine:
c
c  input : (r*8)  x1
c  input : (r*8)  h
c  input : (r*8)  x
c  input : (r*8)  f0
c  input : (r*8)  f1
c  input : (r*8)  g0
c  input : (r*8)  g1
c  input : (r*8)  ei
c  input : (r*8)  ej
c  input : (r*8)  tkij
c  input : (i*4)  li
c  input : (i*4)  lj
c  input : (r*8)  z
c  input : (i*4)  ni
c
c  output: (r*8)  rem
c
c  routines:
c          none
c
c  author:  h. p. summers, university of strathclyde
c           ja7.08
c           tel. 0141-548-4196
c
c  date:    06/06/02
c
c  update: 
C 
C VERSION: 1.1                          DATE: 04-07-96
C MODIFIED: WILLIAM OSBORN
C               - FIRST VERSION.
C
C VERSION: 1.2                          DATE: 19-12-01
C MODIFIED: Martin O'MULLANE
C               - Removed junk from > column 72.c
C
C VERSION: 1.3                          DATE: 18-03-03
C MODIFIED: Hugh Summers
C               -Re-written and documented.
C
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
       integer    li    , lj    , ni   , i
c-----------------------------------------------------------------------
       real*8     x1    , h     , h1   , hh   , x      , f0    , f1    ,
     &            g0    , g1    , ei   , ej   , tkij   ,
     &            z     , z2    , rem 
       real*8     eli   , elj   , t    , wi    , wi0   , wi1   , wi2   , 
     &            wj    , wj0   , wj1  , wj2   , 
     &            ci0   , cj0   , ci1   , cj1  ,
     &            ci2   , cj2   , a1   , a2    , a3  
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
       hh=h*h                                                           
       h1=0.0833333333*hh                                               
       rem=0.0                                                          
       x=x1                                                             
       eli=li                                                           
       elj=lj                                                           
       z2=-2.0*z                                                        
       wi=-eli*(eli+1.0)                                                
       wj=-elj*(elj+1.0)                                                
       wi2=wi+wi                                                        
       wj2=wj+wj                                                        
       t=1.0/x                                                          
       wi0=ei+(z2+wi*t)*t                                               
       wj0=ej+(z2+wj*t)*t                                               
       t=1.0/(x-h)                                                      
       wi1=ei+(z2+wi*t)*t                                               
       wj1=ej+(z2+wj*t)*t                                               
       ci0=f1*(1.0+(h1-hh)*wi0)-f0*(1.0+h1*wi1)                         
       cj0=g1*(1.0+(h1-hh)*wj0)-g0*(1.0+h1*wj1)                         
       ci1=f1*(1.0+h1*wi0)                                              
       cj1=g1*(1.0+h1*wj0)                                              
       a3=f1*g1*x**ni 
                                                         
    1  continue
    
       do i=1,2                                                       
         x=x+h                                                            
         ci1=ci1+ci0                                                      
         cj1=cj1+cj0                                                      
         t=1.0/x                                                          
         ci2=ei+(z2+wi*t)*t                                               
         cj2=ej+(z2+wj*t)*t                                               
         f0=f1                                                            
         g0=g1                                                            
         f1=ci1/(1.0+h1*ci2)                                              
         g1=cj1/(1.0+h1*cj2)                                              
         ci0=ci0-hh*ci2*f1                                                
         cj0=cj0-hh*cj2*g1                                                
       enddo
                                                             
       a1=a3                                                            
       a2=f0*g0*(x-h)**ni                                               
       a3=f1*g1*x**ni                                                   
       rem=rem+a1+4.0*a2+a3
                                                    
       if((tkij*x.ge.10.0).and.(cj2*x*x+wj2.ge.60.0).and.
     &    (ci2*x*x+wi2.ge.60.0))then                                             
           rem=0.333333333*h*rem
       else
           go to 1
       endif
                                                
       return
                                                                  
      end                                                               
