CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adas2xx/adas211/bbprs1.for,v 1.3 2010/04/09 13:54:33 mog Exp $ Date $Date: 2010/04/09 13:54:33 $
CX
       SUBROUTINE BBPRS1(NDMET,STRING,WNO,CPL,NPT,IPLA,ZPLA)
       IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: BBPRS1 *********************
C
C  PURPOSE:  TO ANALYSE THE TAIL CHARACTER STRING OF A LEVEL DATA LINE
C            OF A SPECIFIC ION FILE INTO WAVE-NUMBER AND SETS OF
C            (PARENT IDENTIFIER, EFFECTIVE ZETA FOR THE PARENT) PAIRS.
C
C  CALLING PROGRAM: ADAS208
C
C  NOTES: DETECT  -  LEVEL WAVE NUMBER WHICH PRECEEDS FIRST '{'
C                 -  SETS OF   PARENT INDEX CONTAINED IN '{.}'
C                              FOLLOWED BY EFFECTIVE ZETA
C         NB. 'X' AS FIRST PARENT ASSIGNMENT MEANS EXCLUDE IONISATION
C             FROM THIS LEVEL.
C             NO PARENT ASSIGNMENT MEANS TAKE LOWEST PARENT WITH
C             ZETA =1.
C             LOWEST PARENT BUT NO ZETA MEANS TAKE ZETA =1.
C             IF THERE IS MORE THAN ONE PARENT THEN ZETA'S MUST BE IN.
C
C
C  SUBROUTINE:
C
C  INPUT : (I*4)  NDMET    =  MAXIMUM NUMBER OF PARENTS
C  INPUT : (C*(*))STRING   =  STRING TO BE PARSED
C
C  OUTPUT: (R*8)  WNO      =  EXCITATION WAVE NUMBER OF LEVEL RELATIVE
C                             TO LOWEST PARENT
C  OUTPUT: (C*1)  CPL      =  LEAD PARENT FOR IONISATION  OR 'X'
C  OUTPUT: (I*4)  NPT      =  NUMBER OF BINDING WAVE NUMBERS DETECTED
C  OUTPUT: (I*4)  IPLA()   =  PARENT INDICES.
C  OUTPUT: (R*8)  ZPLA()   =  EFFECTIVE ZETA FOR PARENT IPLA()
C
C          (L*4)  LSET     =  .TRUE.  -  WAVE NUMBER PART SET
C                             .FALSE. -  WAVE NUMBER PART NOT SET
C          (L*4)  LWNO     =  .TRUE.  -  IN THE WAVE NUMBER PART
C                             .FALSE. -  NOT IN THE WAVE NUMBER PART
C          (L*4)  LPRNT    =  .TRUE.  -  IN A PARENT SPECIFIER
C                             .FALSE. -  NOT IN A PARENT SPECIFIER
C          (L*4)  LZETA    =  .TRUE.  -  IN A ZETA SPECIFIER
C                             .FALSE. -  NOT IN A ZETA SPECIFIER
C          (I*4)  IC       =  GENERAL USE
C          (I*4)  IABT     =  FAILURE NUMBER FROM R8FCTN
C          (I*4)  NCHAR    =  NUMBER OF CHARACTERS IN SUBSTRING
C          (C*15) SSTRNG   =  ISOLATED SUBSTRING
C
C ROUTINES: NONE
C
C AUTHOR:  HP SUMMERS
C          K1/1/57
C          JET EXT. 4941
C
C DATE:    22/06/92
C
C UPDATE:  24/06/96  HP SUMMERS - CHANGED SOUBRTOUNE NAME FROM B8PRS1 TO
C                                 BBPRS1
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
C VERSION: 1.2				    DATE: 17-03-03
C MODIFIED: RICHARD MARTIN
C		INITIALISED LWNO AS .FALSE. 
C
C VERSION : 1.3
C DATE    : 09-04-2010
C MODIFIED: Martin O'Mullane
C           - Change integer*4 to integer.
C
C-----------------------------------------------------------------------
       CHARACTER STRING*(*)   , SSTRNG*15  , CPL*1
C
       INTEGER   NDMET        , I4UNIT
       INTEGER   NPT          , IABT       , IC        , NCHAR
       INTEGER   IPLA(NDMET)
       INTEGER   I4FCTN
C
       LOGICAL   LSET         , LWNO         , LPRNT        , LZETA
C
       REAL*8    WNO
       REAL*8    ZPLA(NDMET)
       REAL*8    R8FCTN
C-----------------------------------------------------------------------
       LSET   = .FALSE.
       LPRNT  = .FALSE.
       LZETA  = .FALSE.
	 LWNO   = .FALSE.
       NPT    = 0
       NCHAR  = 0
C
       DO 10 IC=1,NDMET
        IPLA(IC)=0
        ZPLA(IC)=0.0D0
   10  CONTINUE
C
       DO 100 IC=1,LEN(STRING)
        IF(STRING(IC:IC).NE.' ')THEN
            IF(STRING(IC:IC).EQ.'{') THEN
                IF(LWNO) THEN
                    WNO=R8FCTN(SSTRNG(1:NCHAR),IABT)
                     IF(IABT.GT.0)THEN
                     WRITE(I4UNIT(-1),*)'*** ERROR IN BBPRS1 AT R8FCTN'
                         STOP
                     ENDIF
                    LSET=.TRUE.
                    LWNO=.FALSE.
                    LPRNT=.TRUE.
                    NCHAR=0
                ELSEIF(LZETA) THEN
                    ZPLA(NPT)=R8FCTN(SSTRNG(1:NCHAR),IABT)
                     IF(IABT.GT.0)THEN
                     WRITE(I4UNIT(-1),*)'*** ERROR IN BBPRS1 AT R8FCTN'
                         STOP
                     ENDIF
                    LZETA=.FALSE.
                    LPRNT=.TRUE.
                    NCHAR=0
                ELSEIF(LPRNT) THEN
           WRITE(I4UNIT(-1),*)'*** ERROR IN BBPRS1  -  LINE MISORDERED'
                    STOP
                ELSE
                    LPRNT=.TRUE.
                    NCHAR=0
                ENDIF
            ELSEIF(STRING(IC:IC).EQ.'}')THEN
                NPT=NPT+1
                IF(SSTRNG(1:1).EQ.'X') THEN
                     IPLA(NPT)=0
                     IF(NPT.EQ.1)CPL='X'
                ELSE
                     IPLA(NPT)=I4FCTN(SSTRNG(1:1),IABT)
                     IF(IABT.GT.0)THEN
                     WRITE(I4UNIT(-1),*)'*** ERROR IN BBPRS1 AT I4FCTN'
                         STOP
                     ENDIF
                     IF(NPT.EQ.1)CPL=SSTRNG(1:1)
                ENDIF
                LPRNT=.FALSE.
                NCHAR=0
            ELSEIF(LWNO.OR.LPRNT.OR.LZETA)THEN
                NCHAR=NCHAR+1
                SSTRNG(NCHAR:NCHAR)=STRING(IC:IC)
            ELSEIF(.NOT.LSET)THEN
                LWNO=.TRUE.
                SSTRNG(1:1)=STRING(IC:IC)
                NCHAR=1
            ELSEIF(LSET)THEN
                LZETA=.TRUE.
                SSTRNG(1:1)=STRING(IC:IC)
                NCHAR=1
            ENDIF
        ELSEIF(STRING(IC:IC).EQ.' '.AND.LWNO)THEN
            WNO=R8FCTN(SSTRNG(1:NCHAR),IABT)
            IF(IABT.GT.0)THEN
                WRITE(I4UNIT(-1),*)'*** ERROR IN BBPRS1 AT R8FCTN'
                STOP
            ENDIF
            LSET=.TRUE.
            LWNO=.FALSE.
            NCHAR=0
        ELSEIF(STRING(IC:IC).EQ.' '.AND.LZETA)THEN
            ZPLA(NPT)=R8FCTN(SSTRNG(1:NCHAR),IABT)
            IF(IABT.GT.0)THEN
                WRITE(I4UNIT(-1),*)'*** ERROR IN BBPRS1 AT R8FCTN'
                STOP
            ENDIF
            LZETA=.FALSE.
            NCHAR=0
        ENDIF
  100  CONTINUE
C
       IF(NPT.EQ.0) THEN
           NPT=1
           IPLA(1)=1
           ZPLA(1)=1.0D0
           CPL='1'
       ENDIF
       IF(NPT.EQ.1.AND.ZPLA(1).LE.0.0D0)THEN
           ZPLA(1)=1.0D0
       ENDIF
       IF(NPT.EQ.1.AND.IPLA(1).EQ.0)THEN
           IPLA(1)=1
           ZPLA(1)=0.0D0
       ENDIF
       RETURN
      END
