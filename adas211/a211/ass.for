C UNIX-IDL - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/ass.for,v 1.3 2004/07/06 11:17:02 whitefor Exp $ Date $Date: 2004/07/06 11:17:02 $
C
       subroutine ass(a10,a1,a20,a2,phi1,phi2,x,n,e1,e2,nmax,rem)       
       implicit none
c-----------------------------------------------------------------------
c
c  ****************** fortran77 subroutine: ass   **********************
c
c  purpose:  calculates asymptotic integrals using power series modules 
c
c
c  subroutine:
c
c  input : (r*8)  a10
c  input : (r*8)  a1
c  input : (r*8)  a20
c  input : (r*8)  a2
c  input : (r*8)  phi1
c  input : (r*8)  phi2
c  input : (r*8)  x
c  input : (i*4)  n
c  input : (r*8)  e1
c  input : (r*8)  e2
c  input : (i*4)  nmax
c
c  output: (r*8)  rem
c
c
c  routines:
c          routine    source    brief description
c          -------------------------------------------------------------
c          dnaq       adas      
c          dnprod     adas      
c          i4unit     adas      fetch unit number for output of messages
c
c  author:  h. p. summers, university of strathclyde
c           ja7.08
c           tel. 0141-548-4196
c
c  date:    06/06/02
c
c  update:  

C VERSION: 1.1                          DATE: 04-07-96
C MODIFIED: WILLIAM OSBORN
C               - FIRST VERSION.
C VERSION: 1.2                          DATE: 19-12-01
C MODIFIED: Martin O'MULLANE
C               - Removed junk from > column 72.
C
C VERSION: 1.3                          DATE: 18-03-03
C MODIFIED: Hugh Summers
C               - re-written and documented
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
       integer  nmax   , n       
       integer  i      , imax   , j       , jmax
c-----------------------------------------------------------------------
       real*8   a10    , a20    , phi1    , phi2   , x      ,
     &          e1     , e2     , rem     , den
       real*8   b0     , c0     , d0      , e0     ,
     &          f0     , g0     , h0      , p0     , q0     ,   
     &          r      , r1     , r2      , r3     , r4     , s      ,
     &          t      , t1     , t2      , u
       real*8   c1     , c2     , s1      , s2   
c-----------------------------------------------------------------------
       real*8   a1(20)    , a2(20)    , b(20)    , c(20)    , d(20)    ,
     &          e(20)     , f(20)     , g(20)    , h(20)    , p(20)    ,
     &          q(20)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

       call dnaq(a10,a1,b0,b,-2.0d0,nmax,2)                             
       call dnaq(a20,a2,c0,c,-2.0d0,nmax,2)
                                    
       d0=b0+c0                                                         
       c0=b0-c0
                                                                
       do i=1,nmax                                                    
         d(i)=b(i)+c(i)                                                   
         c(i)=b(i)-c(i)                                                   
       enddo
                                                             
       call dnaq(c0,c,b0,b,-1.0d0,nmax,3)                               
       call dnaq(d0,d,c0,c,-1.0d0,nmax,3)                               
       call dnprod(a10,a1,a20,a2,d0,d,nmax)
                                    
       e0=x**n                                                          
       den=n                                                            
       e(1)=den*e0/x
                                                           
       do i=2,nmax                                                    
         den=den-1.0                                                      
         e(i)=den*e(i-1)/x                                                
       enddo
                                                                
       call dnprod(d0,d,e0,e,f0,f,nmax)                                 
       call dnprod(b0,b,f0,f,d0,d,nmax)                                 
       call dnprod(c0,c,f0,f,e0,e,nmax)
                                        
       g0=d0                                                            
       h0=e0                                                            
       imax=nmax-1
                                                             
       do i=1,imax                                                    
         d0=d(1)                                                          
         e0=e(1)                                                          
         jmax=nmax-i                                                      
         do j=1,jmax                                                    
           d(j)=d(j+1)                                                      
           e(j)=e(j+1)                                                      
         enddo
                                                             
         call dnprod(b0,b,d0,d,p0,p,jmax)                                 
         call dnprod(c0,c,e0,e,q0,q,jmax)
	                                  
         d0=p0                                                            
         e0=q0
	                                                             
         do j=1,jmax                                                    
           d(j)=p(j)                                                        
           e(j)=q(j)                                                        
         enddo 
                                                            
         g(i)=d0                                                          
         h(i)=e0
	                                                           
       enddo
                                                                
       g(nmax)=b0*d(1)                                                  
       h(nmax)=c0*e(1)                                                  
       s= 1.0                                                           
       u=-g0                                                            
       r1=u                                                             
       t1=dabs(r1)                                                      
       i=0                                                              
    6  i=i+2                                                            
       if(i-nmax)7,7,9                                                  
    7  t2=dabs(g(i))                                                    
       if(t2-t1)8,8,9                                                   
    8  u=s*g(i)                                                         
       r1=r1+u                                                          
       s=-s                                                             
       t1=t2                                                            
       go to 6                                                          
    9  r1=r1-0.5*u                                                      
       s= 1.0                                                           
       u=-g(1)                                                          
       r2=u                                                             
       t1=dabs(r2)                                                      
       i=1                                                              
   11  i=i+2                                                            
       if(i-nmax)12,12,14                                               
   12  t2=dabs(g(i))                                                    
       if(t2-t1)13,13,14                                                
   13  u=s*g(i)                                                         
       r2=r2+u                                                          
       s=-s                                                             
       t1=t2                                                            
       go to 11                                                         
   14  r2=r2-0.5*u                                                      
       s=-1.0                                                           
       u=h0                                                             
       r3=u                                                             
       t1=dabs(r3)                                                      
       i=0                                                              
   16  i=i+2                                                            
       if(i-nmax)17,17,19                                               
   17  t2=dabs(h(i))                                                    
       if(t2-t1)18,18,19                                                
   18  u=s*h(i)                                                         
       r3=r3+u                                                          
       s=-s                                                             
       t1=t2                                                            
       go to 16                                                         
   19  r3=r3-0.5*u                                                      
       s=-1.0                                                           
       u=h(1)                                                           
       r4=u                                                             
       t1=dabs(r4)                                                      
       i=1                                                              
   21  i=i+2                                                            
       if(i-nmax)22,22,24                                               
   22  t2=dabs(h(i))                                                    
       if(t2-t1)23,23,24                                                
   23  u=s*h(i)                                                         
       r4=r4+u                                                          
       s=-s                                                             
       t1=t2                                                            
       go to 21                                                         
   24  r4=r4-0.5*u                                                      
       s1=dsin(phi1-phi2)                                               
       c1=dcos(phi1-phi2)                                               
       s2=dsin(phi1+phi2)                                               
       c2=dcos(phi1+phi2)                                               
       rem=0.5*(r1*s1+r2*c1+r3*s2+r4*c2)
                                       
       return 
                                                                 
      end                                                               
