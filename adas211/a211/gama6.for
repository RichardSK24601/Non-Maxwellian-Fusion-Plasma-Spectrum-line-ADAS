CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/gama6.for,v 1.2 2004/07/06 13:57:50 whitefor Exp $ Date $Date: 2004/07/06 13:57:50 $
CX
       REAL*8 FUNCTION GAMA6(X)
       IMPLICIT REAL*8(A-H,O-Z)
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
C VERSION  : 1.2                          
C DATE     : 19-12-2001
C MODIFIED : Martin O'Mullane
C               - Changed function defintion to a more standard form.
C               - Removed mainframe listing information beyond column 72.
C
C-----------------------------------------------------------------------
       Z=X
       T=1.0
    1  IF(Z-1.0)2,5,3
    2  T=T/Z
       Z=Z+1.0
       GO TO 1
    3  IF(Z-2.0)5,5,4
    4  Z=Z-1.0
       T=T*Z
       GO TO 3
    5  A=Z-1.0
       GAMA6=T*(1.0-A*(0.57710166-A*(0.98585399-A*(0.87642182-A*(0.83282
     212-A*(0.5684729-A*(0.25482049-A*0.0514993)))))))
       RETURN
      END
