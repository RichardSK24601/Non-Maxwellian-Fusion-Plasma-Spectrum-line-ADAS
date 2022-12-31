CX  UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/system/xxname.for,v 1.4 2007/06/05 14:11:16 allan Exp $ Date $Date: 2007/06/05 14:11:16 $
CX

      SUBROUTINE XXNAME(REALNAME)
      IMPLICIT NONE

C-----------------------------------------------------------------------
C
C  ****************** FORTRAN77 SUBROUTINE: XXNAME *********************
C
C PURPOSE: TO DETERMINE THE REAL NAME OF THE USER BY EXAMINING THE
C          SYSTEM /etc/passwd FILE. A C PROGRAM READS THE FILE.
C   
C CALLING PROGRAM: GENERAL USE
C
C OUTPUT: (C*30)  REALNAME    = REAL NAME OF USER IF IT IS RECORDED
C                                OTHERWISE A DEFAULT STRING IS RETURNED
C ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          GETENV     SYSTEM    GETS USERNAME (8 LETTERS MAX)
C          XXSLEN     ADAS      FINDS NON BLANK POSITIONS IN STRING
C          PERS_F     ADAS      C ROUTINE TO INTEROGATE SYSTEM
C
C
C AUTHOR   : Martin O'Mullane
C DATE     : 11/08/98
C VERSION  : 1.1                          DATE: 11-08-98
C MODIFIED : Martin O'Mullane
C            FIRST VERSION
C
C VERSION:	1.1					DATE: 01-12-98
C MODIFIED: RICHARD MARTIN
C		- PUT UNDER SCCS CONTROL
C
C VERSION:	1.2					DATE: 15-12-98
C MODIFIED: Martin O'Mullane
C		- Replace USER with LOGNAME as location of username in the
C                 environment variables. Linux, at least at JET, does not
C                 assign USER. Other OS appear to set both.
C
C VERSION:	1.3					DATE: 19-01-06
C MODIFIED: Allan Whiteford
C		- Changed test on REALNAME to reflect changes in
C                 underlying C code. Also moved removal of last
C                 character to after 'Who produced this file' is
C                 possibly set.
C
C VERSION:	1.4					DATE: 06-12-06
C MODIFIED: Allan Whiteford
C		- Updated to allow for USERIDs > 8 characters (now set
C                 to 20).
C
C VERSION:	1.5					DATE: 05-07-07
C MODIFIED: Allan Whiteford
C		- Add on CHAR(0) to username as C style string
C                 terminator rather than '\0'
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      EXTERNAL      PERS_F
C-----------------------------------------------------------------------
      INTEGER       L1       , L2
C-----------------------------------------------------------------------
      CHARACTER*20  USERID
      CHARACTER*21  D
      CHARACTER*30  REALNAME , PERS_F
C-----------------------------------------------------------------------
      
      CALL GETENV("LOGNAME",USERID)
      CALL XXSLEN(USERID,L1,L2)

      IF (USERID.NE.' ') THEN
         D = USERID(1:L2)//CHAR(0) 
         REALNAME = PERS_F(D)
         IF (REALNAME(1:1).EQ.' ') REALNAME='Who produced this file?? '
         CALL XXSLEN(REALNAME,L1,L2)    ! remove C null terminator
         REALNAME = REALNAME(1:L2-1)
      ELSE
         REALNAME = 'Who produced this file?'
      ENDIF

      END
