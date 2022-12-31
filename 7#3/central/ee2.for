CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/maths/ee2.for,v 1.3 2007/04/23 13:40:45 allan Exp $ Date $Date: 2007/04/23 13:40:45 $
CX
       FUNCTION EE2(X)                                                  
C
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 FUNCTION: EE2 **************************
C
C PURPOSE: EVALUATES EXP(X)E2(X) WHERE E2 IS THE 2ND EXPONENTIAL
C          INTEGRAL
C
C CALLING PROGRAMS: GENERAL
C
C INPUT:  (R*8)  X       = INDEPENDENT VARIABLE
C
C OUTPUT: (R*8)  EE2     = EXP(X)E2(X) 
C
C ROUTINES:
C          ROUTINE    SOURCE   BRIEF DESCRIPTION
C          -------------------------------------------------------------
C          EEI        ADAS     EVALUATES 1ST EXPONENTIAL INTEGRAL
C
C UNIX-IDL PORT:
C
C VERSION: 1.1                          DATE: 11-07-95
C MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
C               - PUT UNDER S.C.C.S. CONTROL
C
C VERSION: 1.2                          DATE: 06-03-96
C MODIFIED: TIM HAMMOND 
C               - ADDED HEADERS
C
C VERSION: 1.3                          DATE: 17-4-07
C MODIFIED: HUGH SUMMERS
C           - COMPLETED COMMENT BLOCK DESCRIPTION
C
C-----------------------------------------------------------------------
C
       IMPLICIT REAL*8(A-H,O-Z)                                         
       IF(X-30.0D0)1,1,2                                                
    1  EE2=1.0D0-X*EEI(X)                                               
       GO TO 3                                                          
    2  X1=1.0D0/X                                                       
       EE2=X1*(1.0D0-X1*(2.0D0-X1*(6.0D0-X1*(24.0D0-X1*(120.0D0-X1*     
     1(720.0D0-X1*5040.0D0))))))                                        
    3  RETURN                                                           
      END                                                               
