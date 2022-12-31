CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/bbprs3.for,v 1.3 2007/05/17 17:22:23 allan Exp $ Date $Date: 2007/05/17 17:22:23 $
CX
       SUBROUTINE BBPRS3( STRING, IA, LCLSHL )
       IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: BBPRS3 *********************
C
C  PURPOSE:  TO ANALYSE A CONFIGURATION CHARACTER STRING INTO A  INTEGER
C            ARRAY OF OCCUPATION NUMBERS IN A STANDARD ORDER
C
C  CALLING PROGRAM: ADAS211
C
C  NOTES: THE STANDARD ORDER IS 1S,2S,2P,3S,3P,3D ......, 4F (15 VALUES)
C         CLOSED SHELLS WITHIN THE ACTIVE N-SHELLS ARE ASSUMED FULLY
C         OCCUPIED
C
C
C  SUBROUTINE:
C
C  INPUT : (C*(*))STRING   =  STRING TO BE PARSED
C          (L*4)  LCLSHL   =  SWITCH ON CLOSED SHELL ASSUMPTION
C
C  OUTPUT: (I*4)  IA()     =  SET OF OCCUPATION NUMBERS IN STANDARD
C                             ORDER
C
C ROUTINES: NONE
C
C AUTHOR:  HP SUMMERS
C          K1/1/57
C          JET EXT. 4941
C
C DATE:    29/06/92
C
C UPDATE:  W.J. DICKSON  7/10/92
C          ADDED PARAMETER LCLSHL TO SWITCH OFF CLOSED SHELL
C          APPROXIMATION
C UPDATE:  H. P. SUMMERS  1/10/96
C          PERMITTED LOWER AND UPPER CASE ORBITAL L-VALUES
C          IN CONFIGURATION STRINGS. DETECT RETURNED L<0
C
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
C VERSION: 1.2                          DATE: 14-10-96
C MODIFIED: WILLIAM OSBORN
C               - ADDED CHANGES DATED 1/10/96 ABOVE
C
C VERSION: 1.3                          DATE: 17-05-07
C MODIFIED: Allan Whiteford
C               - Removed non-standard control character from
C                 comments.
C
C-----------------------------------------------------------------------
       CHARACTER STRING*(*)
C-----------------------------------------------------------------------
       LOGICAL    LCLSHL
C-----------------------------------------------------------------------
       INTEGER*4  N      , L     , I     , NC    , IC    , IABT
       INTEGER*4  NMIN   , L1    , LUC   , LLC
       INTEGER*4  IA(15)
       INTEGER*4  I4FCTN , LEN   , IND   , INDEX , I4UNIT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
       IND(N,L)=(N*(N-1))/2+L+1
C-----------------------------------------------------------------------
       NMIN = 6
       DO 20 I=1,15
         IA(I)=0
   20  CONTINUE
C
       NC=LEN(STRING)
       IC=0
   30  IC=IC+1
       IF(IC.LE.NC-2)THEN
           IF(STRING(IC:IC).EQ.' ') GO TO 30
           N=I4FCTN(STRING(IC:IC),IABT)
           NMIN=MIN0(NMIN,N)
           IF(IABT.GT.0)THEN
            WRITE(I4UNIT(-1),*)'*** ERROR IN PRS2 AT I4FCTN, IABT=',IABT
               STOP
           ENDIF
           LUC=INDEX('SPDFGHIJK',STRING(IC+1:IC+1))-1
           LLC=INDEX('spdfghijk',STRING(IC+1:IC+1))-1
           L=MAX0(LUC,LLC)
           IF(N.GT.5.OR.L.GE.N.OR.L.LT.0)THEN
            WRITE(I4UNIT(-1),*)'*** ERROR IN PRS2 - N OR L OUR OF RANGE'
               STOP
           ENDIF
           IA(IND(N,L))=I4FCTN(STRING(IC+2:IC+2),IABT)
           IC=IC+2
           GO TO 30
       ENDIF
C
       NMIN=NMIN-1
       IF(NMIN.GT.0)THEN
           DO 50 N=1,NMIN
             DO 40 L1=1,N
               IF( LCLSHL ) THEN
                  IA(IND(N,L1-1))=4*L1-2
               ELSE
                  IA(IND(N,L1-1))=0
               ENDIF
   40        CONTINUE
   50      CONTINUE
       ENDIF
C
       RETURN
      END
