CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/delta2.for,v 1.3 2004/07/06 13:30:22 whitefor Exp $ Date $Date: 2004/07/06 13:30:22 $
CX
       SUBROUTINE DELTA2(A,B,C,D,J)
       IMPLICIT REAL*8 (A-H,O-Z)
C-----------------------------------------------------------------------
C UNIX-IDL PORT:
C
C AUTHOR:  WILLIAM OSBORN (TESSELLA SUPPORT SERVICES PLC)
C
C DATE:    4TH JULY 1996
C
C VERSION: 1.1                          DATE: 04-07-96
C MODIFIED: WILLIAM OSBORN
C               - FIRST VERSION.
C
C VERSION: 1.2                          DATE: 19-12-01
C MODIFIED: Martin O'MULLANE
C               - Removed junk from > column 72.
C
C VERSION: 1.3                          DATE: 18-03-03
C MODIFIED: Richard Martin
C               - Increased dimensions to 500.
C
C-----------------------------------------------------------------------
       DIMENSION GAM(500),JGAM(500)
       COMMON /FAC/GAM,JGAM
       A1=A+B-C+0.001
       IF(A1)5,1,1
    1  A2=A-B+C+0.001
       IF(A2)5,2,2
    2  A3=B-A+C+0.001
       IF(A3)5,3,3
    3  J1=A1
       J2=A2
       J3=A3
       J4=J1+J2+J3+1
       J5=A1+A2+A3+1.501
       IF(J4-J5)5,4,5
    4  D=GAM(J1+1)*GAM(J2+1)*GAM(J3+1)/GAM(J4+1)
       J=JGAM(J1+1)+JGAM(J2+1)+JGAM(J3+1)-JGAM(J4+1)
       RETURN
    5  D=0.0
       J=0
       RETURN
      END
