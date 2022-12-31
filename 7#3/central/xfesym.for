CX UNIX PORT - SCCS info: Module @(#)$Header: /home/adascvs/fortran/adaslib/atomic/xfesym.for,v 1.2 2004/07/06 15:30:12 whitefor Exp $ Date $Date: 2004/07/06 15:30:12 $
CX
      FUNCTION XFESYM ( IZ0 )
      IMPLICIT NONE
C-----------------------------------------------------------------------
C
C ************** FORTRAN77 CHARACTER*2 FUNCTION: XFESYM ****************
C
C PURPOSE: TO RETURN THE SYMBOL FOR THE ELEMENT WITH NUCLEAR CHARGE IZ0
C          (CHARACTER*2 FUNCTION VERSION OF 'XXESYM')
C
C CALLING PROGRAM: GENERAL USE
C
C FUNCTION:
C
C          (C*2)  XFESYM  = FUNCTION NAME -
C                           SYMBOL OF ELEMENT WITH NUCLEAR CHARGE 'IZ0'
C          (I*4)  IZ0     = ELEMENT NUCLEAR CHARGE
C
C          (C*2)  SYMBOL()= SYMBOLS OF FIRST 92 ELEMENTS.
C                           ARRAY DIMENSION => NUCLEAR CHARGE
C
C NOTES:    IF NUCLEAR CHARGE IS OUT OF RANGE, I.E.NOT BETWEEN 1 & 92,
C           THEN THE CHARACTER STRING 'XFESYM' IS RETURNED BLANK.
C
C ROUTINES: NONE
C
C
C AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
C           K1/0/81
C           JET EXT. 4569
C
C DATE:     12/02/91
C
C UPDATES:  25/10/94 L. JALOTA (TESSELLA SUPPORT SERVICES PLC)
C		     CHANGED CASE OF SYMBOL TO LOWER CASE FOR UNIX
C
C VERSION:	1.2 
C UPDATES:  17/09/99 HUGH SUMMERS - INCREASED ELEMENT NUMBER TO 92 
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER      IZ0
C-----------------------------------------------------------------------
      CHARACTER*2  XFESYM , SYMBOL(92)
C-----------------------------------------------------------------------
      DATA SYMBOL/'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     &            'na','mg','al','si','p ','s ','cl','ar','k ','ca',
     &            'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     &            'ga','ge','as','se','br','kr','rb','sr','y ','zr',
     &            'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     &            'sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     &            'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     &            'lu','hf','ta','w ','re','os','ir','pt','au','hg',
     &            'tl','pb','bi','po','at','rn','fr','ra','ac','th',
     &            'pa','u '/
C-----------------------------------------------------------------------
         IF ( (IZ0.GT.92).OR.(IZ0.LT.0) ) THEN
            XFESYM = ' '
         ELSE
            XFESYM = SYMBOL(IZ0)
         ENDIF
C-----------------------------------------------------------------------
      RETURN
      END
