       IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ******************* FORTRAN77 PROGRAM: ADAS211 **********************
C
C  ORIGINAL NAME: ADASRRC (RADRCDW)
C
C  VERSION:  1.0
C
C  PURPOSE:  TO CALCULATE RADIATIVE RECOMBINATION RATES TO LS TERMS OF
C            AN ION IN A MANNER COMPATIBLE WITH THE ADAS STRUCTURE
C
C  NOTES:  USES DISTORTED WAVES IN A MODEL EXTENSION OVER THAT IN
C          JETSHP.RECOM.FORT(NRADRC).  IT CALLS BBPHOT AND SO GIIDW
C          INSTEAD OF PHOTO5 AND GIIBS.  EXTERNAL FUNCTION RESOLUTION
C          MUST BE FROM DLS.LOAD, ABLIB3.LOAD AND RECOM.LOAD IN THAT
C          ORDER TO FORCE ABLIB3 VERSIONS OF GAMAF AND WIG6J.
C
C  N.B.    OPTIONS 1,2 AND 3 ARE NOT WORKING AT THE MOMENT
C
C  DATA:   THE SOURCE DATA IS CONTAINED IN PARTITIONED DATASETS AS
C          FOLLOWS:-
C
C                  'JETSHP.RADREC.<SE>LIKE(<MEMBER>)'
C          WHERE, <SE> = ELEMENT SYMBOL FOR RECOMBINING ION ISO-
C                        ELECTRONIC SEQUENCE (ONE OR TWO CHARACTERS)
C
C  PROGRAM:
C          (I*4) NDMET     = PARAMETER = MAXIMUM NUMBER OF METASTABLES
C          (I*4) NDLEV     = PARAMETER = MAXIMUM NUMBER OF LEVELS
C          (I*4) NDTEM     = PARAMETER = MAXIMUM NUMBER OF TEMPERATURES
C          (I*4) IUNT05    = INPUT UNIT FROM BATCH FILE.
C          (I*4) IUNT17    = OUTPUT UNIT FOR RUN SUMMARY.
C          (I*4) IUNT09    = PARAMETER = UNIT NUMBER FOR I/O
C          (I*4) IUNT10    = PARAMETER = UNIT NUMBER FOR I/O
C          (I*4) IUNT11    = PARAMETER = UNIT NUMBER FOR I/O
C          (I*4) IUNT12    = PARAMETER = UNIT NUMBER FOR I/O
C          (I*4) IUNT37    = PARAMETER = UNIT NUMBER FOR I/O
C          (I*4) IUNT20    = TEMPORARY SCRATCH FILE
C          (I*4) NHCUT     = PARAMETER = HIGHEST N-SHELL TREATED NON-
C                                        HYDROGENICALLY.
C          (i*4) dist      = distribution type
C                              0 => Maxwellian
C                              1 => kappa
C                              2 => numeric
C                              3 => Druyvesteyn
C          (i*4) ndgnt     = max no of Gaunt values
C          (i*4) nemax     = max no of adf37 energies
C          (i*4) ntmax     = max no of adf37 temperatures
C          (i*4) icateg    = category of adf37 file
C                              1 => superposition
C                              2 => numerical
C          (i*4) nenerg    = icateg = 1 => number of distribution families
C                            icateg = 2 => number of energy points
C          (i*4) nblock    = icateg = 1 => number of members in output family
C                            icateg = 2 => number of effective temperatures
C          (i*4) nform1    = type of threshold behaviour
C                              1 => cutoff
C                              2 => energy^param1
C          (r*8) param1    = parameter of threshold form
C          (i*4) nform2    = type of high-energy behaviour
C                              1 => cutoff
C                              2 => energy^-param2(1)
C                              3 => exp(-param2(1)*energy)
C                              4 => exp(-param2(1)*energy^param2(2))
C          (r*8) param2()  = parameter of high-energy form
C          (i*4) maxe      = no of gaunt/vve pairs
C          (R*8) QDMIN     = PARAMETER = MINIMUM QUANTUM DEFECT AT WHICH
C                                        TO HYDROGENIC MATRX ELEMENTS OCCURS
C          (I*4) LEN3      = GENERAL CHARACTER STRING LENGTH INDEX
C          (I*4) LEN4      = GENERAL CHARACTER STRING LENGTH INDEX
C          (I*4) IZ0       = NUCLEAR CHARGE
C          (I*4) IZ1       = RECOMBINING ION CHARGE
C          (I*4) NPRNTI    = NUMBER OF INITIAL RECOMBINING PARENTS
C          (I*4) I         = GENERAL INDEX
C          (I*4) IP        = GENERAL PARENT INDEX
C          (I*4) ITRM      = GENERAL TERM INDEX
C          (I*4) NTRM      = NUMBER OF TERMS
C          (I*4) IT        = GENERAL TEMPERATURE INDEX
C          (I*4) MAXT      = NUMBER OF TEMPERATURES
C          (I*4) IRESOL    = NUMBER OF TEMPERATURES
C          (I*4) IRESOL    = 1 FOR  ((LP,SP)N L LT S,(LP,SP)L1 LT1 S)
C                          = 2 FOR  ((LP,SP)N L LT S,(LP,SP)L1 S)
C                                   ==> SUM OVER LT1 OF 1.
C                          = 3 FOR  ((LP,SP)N L S,(LP,SP)L1 S)
C                                   ==> SUM OVER LT  OF 2.
C                          = 4 FOR  ((LP,SP)N L,(LP,SP)L1)
C                                   ==> SUM OVER S  OF 3.
C                          = 5 FOR  NO L RESOLUTION USING GBF
C          (I*4) IBSOPT    = 1 => USE FITTED PEACH PHASE IN MATRIX
C                                 ELEMENT
C                          = 2 => USE BURGESS-SEATON PHASE IN MATRIX
C                                 ELEMENT
C                          = 3 => USE HYDROGENIC APPROXIMATION
C                          = 4 => USE DISTORTED WAVE MATRIX ELEMENT
C                          (IBSOPT ACTS THROUGH THE /BSPARS/
C                           COMMON BLOCK)
C          (I*4) NSHELL    = NUMBER OF ORBITAL SHELLS IN CONFIGURATION
C          (I*4) ISHL      = SHELL INDEX
C          (I*4) NL        = NL COMPOSITE FROM PARSING OF SHELL CODE
C          (I*4) LP        = PARENT TOTAL ORBITAL ANGULAR MOMENTUMDE
C          (I*4) ISP       = PARENT MULTIPLICITY (2*SP+1)
C          (I*4) LVCT      = L-CUT FOR HYDROGENIC APPROXIMATION
C          (I*4) ISCORE    = NUMBER OF IMPLICIT CLOSED CORE SHELLS
C          (I*4) N         = PRINCIPAL QUANTUM NUMBER
C          (I*4) L         = ORBITAL ANGULAR MOMENTUM QUANTUM NUMBER
C          (I*4) IPP       =
C          (I*4) IS        =
C          (I*4) LT        = TERM ORBITAL ANGULAR MOMENTUN QU. NO.
C          (I*4) NMAX      =
C          (I*4) LAM       = MULTIPOLE
C          (I*4) JH        = NUMBER OF RADIAL WAVE FUNCTION ORDINATES
C          (I*4) JEALFA    =
C          (I*4) JSN       =
C          (I*4) IFIRST    =
C          (I*4) IGONE     =
C          (I*4) IONCE     =
C          (I*4) IEXT      =
C          (I*4) IREPT     =
C          (I*4) LMIN      =
C          (I*4) LMAX      =
C          (I*4) L1        =
C          (I*4) LT1       =
C          (I*4) IWARN     =
C          (R*8) Z0        = NUCLEAR CHARGE
C          (R*8) Z1        = RECOMBINING ION CHARGE
C          (R*8) ACC       = ACCURACY OF SEARCH FOR ALPHA PARAMETER
C          (R*8) XMAX      = RANGE OF RADIAL WAVE FUNCTIONS
C          (R*8) H         = WAVE FUNCTION EVALUATION STEP LENGTH
C          (R*8) U1        = QUANTUM DEFECT EXPANSION COEFFICIENT
C          (R*8) U2        = QUANTUM DEFECT EXPANSION COEFFICIENT
C          (R*8) U3        = QUANTUM DEFECT EXPANSION COEFFICIENT
C          (R*8) ZETA      = QUANTUM DEFECT EXPANSION COEFFICIENT
C          (R*8) BWNR      = IONISATION POT. FOR CURRENT PARENT(CM-1)
C          (R*8) BWNO      = IONISATION POTENTIAL (CM-1)
C          (R*8) FRACPR    = FRACTIONAL PARENTAGE COEFFICIENT (APPROX?)
C          (R*8) EN        = PRINCIPAL QUANTUM NUMBER
C          (R*8) V         = EFFECTIVE PRINCIPAL QUANTUM NUMBER
C          (R*8) QD        = QUANTUM DEFECT
C          (R*8) TE        = ELECTRON TEMPERATURE (K)
C          (R*8) B         = 1.5789D5*Z**2/(V**2*TE)
C          (R*8) B1        = 1.5789D5*Z**2/(V**2*TR)
C                            (TE=ELECTRON TEMP. (K), TR=RAD. TEMP. (K),
C                             Z= BOUND STATE ION CHARGE +1)
C          (R*8) R1        = ESTIMATE OF ATOMIC RADIUS
C          (R*8) ANS       =
C          (R*8) PREC      = RADIATIVE RECOMBINATION INTEGRAL
C          (R*8) PION      = PHOTOIONISATION INTEGRAL
C          (R*8) PSTIM     = STIMULATED RECOMBINATION INTEGRAL
C          (r*8) dparam    = parameter of distribution
C                              dist = 1 => kappa
C                              dist = 3 => x
C          (r*8) rrcint    = radiative recombination integral
C          (r*8) rrcin()   = rrcint at adf37 temperatures
C          (r*8) rrcspl()  = rrcint at adf08 temperatures
C                ea(,)     = adf37 energy points of tabulation
C          (r*8) fa(,)     = distribution function tabulation
C          (r*8) teff()    = adf37 effective temperatures (eV)
C          (r*8) tin()     = adf37 effective temperatures (K)
C          (r*8) mode()    = most probable energy (eV)
C          (r*8) median()  = median energy (eV)
C          (r*8) ein()     = adf37 energies
C          (r*8) fin()     = distribution function at ein
C          (r*8) eout()    = Gaunt factor energies
C          (r*8) fout()    = distribution function at eout
C        (C*133) STRING    = GENERAL LONG STRING
C         (C*50) SSTRNG1   = GENERAL SHORT STRING
C        (C*133) SSTRNG    = GENERAL MEDIUM STRING
C          (C*8) BLANKS    = BLANK STRING
C         (C*20) BLNK20    = BLANK STRING
C          (C*2) SEQ       = ELEMENT CHEMICAL SYMBOL
C          (C*2) SEQUP     = ELEMENT CHEMICAL SYMBOL IN UPPERCASE
C        (C*100) SEQSS     = CONCATENATED ELEMENT SYMBOLS
C         (C*30) SHLSS     =
C          (C*2) SLSS      =
C         (C*80) DSNIN     = INPUT DATA FILE NAME
C         (C*80) DSNOUT    = OUTPUT DATA FILE NAME FOR STREAM 10
C         (C*80) DSNPAS    = OUTPUT DATA FILE NAME FOR STREAM 11
C         (C*80) DSNPAP    = OUTPUT DATA FILE NAME FOR STREAM 17
C         (C*80) DSN80     = GENERAL DATA SET NAME STRING
C         (C*10) DATE      =
C        (c*120) file37    = adf37 filename
C         (c*80) titl37    = header for adf37 file
C        (c*120) filout    = file name of output family
C        (c*120) filnam()  = file names of input families
C         (c*25) calgeb(,) = distribution function algebra
C         (c*25) ealgeb()  = energy parameter algebra
C          (L*4) OPEN09    = .TRUE.  => STREAM 09 DEVICE OPEN
C                          = .FALSE. => STREAM 09 DEVICE NOT OPEN
C          (L*4) OPEN10    = .TRUE.  => STREAM 10 DEVICE OPEN
C                          = .FALSE. => STREAM 10 DEVICE NOT OPEN
C          (L*4) OPEN11    = .TRUE.  => STREAM 11 DEVICE OPEN
C                          = .FALSE. => STREAM 11 DEVICE NOT OPEN
C          (l*4) open37    = .true.  => stream 37 device open
C                          = .false. => stream 37 device not open
C          (L*4) LCLSHL    = .TRUE.  => IMPLICIT CLOSED SHELLS PRESENT
C                          = .FALSE. => NO IMPLICIT CLOSED SHELLS
C          (l*4) lexist    = .true.  => adf37 file exists
C                          = .false. => adf37 file does not exist
C          (I*4) ISPA()    = PARENT TERM MULTIPLICITIES
C                            1ST DIM: PARENT INDEX
C          (I*4) LTPA()    = PARENT TERM ORBITAL QUANTUM NUMBER
C                            1ST DIM: PARENT INDEX
C          (I*4) ICNFPA(,) =
C                            1ST DIM: SHELL INDEX
C                            2ND DIM: PARENT INDEX
C          (I*4) ISA()     =  TERM MULTIPLICITIES
C                            1ST DIM: TERM INDEX
C          (I*4) LTA()     = TERM ORBITAL QUANTUM NUMBER
C                            1ST DIM: TERM INDEX
C          (I*4) ICNFA(,)  =
C                            1ST DIM: SHELL INDEX
C                            2ND DIM: TERM INDEX
C          (I*4) NPTA()    =
C                            1ST DIM: TERM INDEX
C          (I*4) IPLA(,)   = PARENT TERM ORBITAL QUANTUM NUMBER
C                            1ST DIM: PARENT INDEX
C                            2ND DIM: TERM INDEX
C          (I*4) NQLS()    = SHELL CODES FOR CONFIGUTATION
C                            1ST DIM: SHELL  INDEX
C          (I*4) NAA()     = PRINCIPAL QUANTUM NUMBER OF ACTIVE ELECTRON
C                            1ST DIM: INITIAL AND FINAL STATES
C          (I*4) LAA()     = ORBITAL QUANTUM NUMBER OF ACTIVE ELECTRON
C                            1ST DIM: INITIAL AND FINAL STATES
C          (I*4) J0AG()    =
C          (I*4) JCAG()    =
C          (I*4) J0AX()    =
C          (I*4) NA()      =
C          (R*8) WPA()     = PARENT STATISTICAL WEIGHTS
C                            1ST DIM: PARENT INDEX
C          (R*8) WNOPA()   = PARENT WAVE NUMBERS (CM-1)
C                            1ST DIM: PARENT INDEX
C          (R*8) WA()      = TERM STATISTICAL WEIGHTS
C                            1ST DIM: TERM INDEX
C          (R*8) WNOA()    = TERM WAVE NUMBERS (CM-1)
C                            1ST DIM: TERM INDEX
C          (R*8) FPLA(,)   =
C                            1ST DIM: PARENT INDEX
C                            2ND DIM: TERM INDEX
C          (R*8) EAA()     = ENERGY OF ACTIVE ELECTRON (RYD)
C                            1ST DIM: INITIAL AND FINAL STATES
C          (R*8) QDAA()    = QUANTUM DEFECT OF ACTIVE ELECTRON
C                            1ST DIM: INITIAL AND FINAL STATES
C          (R*8) ALFAA(,)  = SCREENING PARAMETERS FOR SHELLS
C                            1ST DIM: INITIAL AND FINAL STATES
C                            2ND DIM: SHELL INDEX
C          (R*8) REAX()    =
C          (R*8) RVAX()    =
C          (R*8) XA(,,)    =
C          (R*8) EAG()     =
C          (R*8) VAG()     =
C          (R*8) GA(,,)    =
C          (R*8) TEA()     = ELECTRON TEMPERATURES (K)
C                            1ST DIM: TEMPERATURE INDEX
C          (R*8) EWNA()    =
C          (R*8) ALF()     =
C          (R*8) SUM()     =
C          (r*8) gaunt()   = bound-free Gaunt factor
C          (r*8) vve()     = v**2*e 
C                              where e = (free electron energy)/z**2 (ryd)
c                                    v = effective principal quantum number 
C                                        of bound electron
C         (C*20) CNFPA()   =
C                            1ST DIM: PARENT INDEX
C         (C*20) CNFA()    =
C                            1ST DIM: TERM INDEX
C          (C*1) CPLA()    =
C                            1ST DIM: TERM INDEX
C          (C*3)  REP      = 'YES' => TERMINATE PROGRAM EXECUTION
C                            'NO'  => CONTINUE PROGRAM EXECUTION
C          (C*40) TITLE     = ISPF ENTERED GENERAL TITLE FOR RUN.
C          (C*120) CADAS    = ADAS HEADER: INCLUDES RELEASE, PROGRAM,
C                            TIME.
C          (C*8)  DATE     = CURRENT DATE (AS 'DD/MM/YY')
C          (L*4)  LCONT    = .TRUE.  => DISPLAY INPUT FILE DESCRIPTIVE
C                                       TEXT
C                            .FALSE. => - DO NOT DO THE ABOVE -
C          (L*4)  LPASS    = .TRUE.  => DISPLAY INPUT FILE DESCRIPTIVE
C                                       TEXT
C                            .FALSE. => - DO NOT DO THE ABOVE -
C          (L*4)  LPAPER   = .TRUE.  => OUTPUT TO PAPER.TXT FILE
C                            .FALSE. => - DO NOT DO THE ABOVE -
C          (L*4)  LBATCH    = FLAGS HOW PROGRAM IS BEING EXECUTED.
C                             .TRUE.  => PROGRAM RUNNING IN BATCH.
C                             .FALSE. => PROGRAM RUNNING IN FOREGROUND.
C          (L*4)  LPEND     = FLAGS IF END OF ANALYSIS REQUESTED.
C                             .TRUE.  => END ANALYSIS OF CURRENT DATA
C                                        SETS
C                             .FALSE. => CONTINUE PANALYSIS WITH CURRENT
C                                        DATA SETS
C          (L*4)  LBTSEL    = FLAGS IF PROGRAM TO BE RUN IN BATCH.
C                             .TRUE.  => RUN PROGRAM IN BATCH.
C                             .FALSE. => RUN PROGRAM IN FOREGROUND.
C          (L*4)  LEISS     = .TRUE.  => EISSNER CONFIG. FORM DETECTED
C                             .FALSE. => EISSNER CONFIG. FORM NOT DETECTED
C          (L*4)  LSTAN     = .TRUE.  => STANDARD CONFIG. FORM DETECTED
C                             .FALSE. => STANDARD CONFIG. FORM NOT DETECTED
C          (L*4)  LTERM     = .TRUE.  => TERM DATA SET TYPE DETECTED
C                             .FALSE. => TERM DATA SET TYPE NOT DETECTED
C          (L*4)  LLEVL     = .TRUE.  => LEVEL DATA SET TYPE DETECTED
C                             .FALSE. => LEVEL DATA SET TYPE NOT DETECTED
C          (L*4)  LPTERM    = .TRUE.  => TERM PARENT DATA SET TYPE DETECTED
C                             .FALSE. => TERM PARENT DATA SET TYPE NOT DETECTED
C          (L*4)  LPLEVL    = .TRUE.  => LEVEL PARENT DATA SET TYPE DETECTED
C                             .FALSE. => LEVEL PARENT DATA SET TYPE NOT DETECTED
C
C
C          (I*4)  LOGIC     = VALUE READ FROM PIPE FOR LOGICAL VARIABLES
C
C          (I*4)  PIPEIN    = STANDARD INPUT
C          (I*4)  PIPEOU    = STANDARD OUTPUT
C
C
C  ROUTINES:
C          ROUTINE    SOURCE    BRIEF DESCRIPTION
C          ------------------------------------------------------------
C          BBSPF0     ADAS      GATHERS INPUT FILE NAMES FROM IDL OR FILE
C          BBSPF1     ADAS      GATHERS OUTPUT FILE NAMES FROM IDL OR FILE
C          XXSLEN     ADAS      GET LENGTH OF STRING VARIABLE
C          BBPRS3     ADAS      PARSE SHELL CODES
C          BBPRS1     ADAS      PARSE LEVEL PARENTAGE STRING
C          BBLPOL     ADAS      EVALUATE MULTIPOLE MATRIX ELEMENT
C          BBPHOT     ADAS      EVALUATE PHOTO-INTEGRALS
C          XX0000     ADAS      SET MACHINE DEPENDANT ADAS CONFIGURATION
C          XXDTES     ADAS      DETECTS EISSNER CONFIGURATION FORM
C          XXCFTR     ADAS      CONVERTS BETWEEN EISSNER/STANDARD CONFIG. FORMS
C          I4UNIT     ADAS      FETCH UNIT NUMBER FOR OUTPUT OF MESSAGES
C          XXFLSH     ADAS      FLUSHES I/O BUFFER
C	   XXDATE     ADAS	OBTAINS DATE FROM IDL
C
C
C  AUTHOR:  HP SUMMERS
C           K1/1/57
C           JET EXT. 4941
C
C  UPDATE:  WJ DICKSON   8TH SEPTEMBER 1992
C           OPENED OUTPUT FILES WITH USER ID
C
C  UPDATE:  WJ DICKSON   7TH OCTOBER 1992
C           CREATED LOGICAL VARIABLE LCLSHL TO SET SWITCH ON
C           CLOSED SHELL APPROXIMATION
C           AJDUSTED AND RENAMED SUBROUTINE PRS2 TO PRS3
C
C  UPDATE:  WJ DICKSON   3RD NOVEMBER 1992
C           INCREASED LENGTH OF CHARACTER SSTRNG1 FROM 40 TO 50
C           TO ENSURE THAT INFO ON THREE PARENTS COULD BE READ
C
C  UPDATE:  WJ DICKSON  18TH JANUARY 1993
C           CORRECTED SUBROUTINE APHOTDW
C           TO CORRECT ERROR IN NUMERICAL QUADRATURE
C           (G-L WEIGHING FACTOR TYPED IN INCORRECTLY)
C
C  UPDATE:  HP SUMMERS  22/05/95  ADJUST FORMATTING FOR READING
C                                 TEMPERATURE LINE SO AS SAME AS OUTPUT
C
C  UPDATE:  HP SUMMERS  16/06/95  ALTER DEFINITION OF NLQS AS
C                                 1000*N+100*L+IQ TO AVOID PROBLEM WHEN\
C                                 NUMBER OF EQUIVALENT ELECTRONS IS 10.
C  UPDATE:  HP SUMMERS  19/06/96  RENAMED AS ADAS211.  ADJUSTED OUTPUT
C                                 FILE NAMES TO STANDARD
C  UPDATE:  HP SUMMERS  24/06/96  RENAMED B8PRS1 TO BBPRS1 AND BROUGHT
C                                 BACK DLS VERSION OF B8PRS1 AS BBPRS1
C                                 CORRECTED DSNOUT, DSNPAS TO *22
C  UPDATE:  HP SUMMERS  24/06/96  INTRODUCE NDTEM FOR MAX. NO. OF TEMPS.
C
C  UPDATE:  HP SUMMERS  02/07/96  EXTENDED TO ALLOW BACKGROUND OR FORE-.
C                                 RUNNING FROM ADAS
C  UPDATE:  HP SUMMERS  25/09/96  IMPLEMENT SWITCH TO HYDROGENIC GAUNT
C                                 FACTORS FOR SMALL AND NEGATIVE QUANTUM
C                                 DEFECTS.  DETECT EISSNER OR STANDARD
C                                 CONFIGURATION FORMS.  DETECT TERM OR
C                                 LEVEL ORGANISATION IN PARENTS AND
C                                 ENERGY LEVELS.
C  UPDATE:  HP SUMMERS  07/11/97  IMPLEMENT SWITCH TO HYDROGENIC GAUNT
C                                 FACTORS FOR N>NHCUT.
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
C VERSION: 1.2                          DATE: 19-08-96
C MODIFIED: WILLIAM OSBORN
C               - CHANGED UNIT 7 TO UNIT 17.
C		- COMMENTED-OUT DIAGNOSTIC PRINT STATEMENTS.
C
C VERSION: 1.3                          DATE: 19-08-96
C MODIFIED: WILLIAM OSBORN
C               - ADDED CHANGES DATED 25/09/96 ABOVE
C
C VERSION: 1.4                          DATE: 17-10-96
C MODIFIED: WILLIAM OSBORN
C               - REMOVED REFERENCES TO CSTRGO(X:20) SINCE IT IS OF
C                 TYPE CHARACTER*19
C
C VERSION: 1.5				DATE: 13-02-97
C MODIFIED: RICHARD MARTIN
C		- CHANGED NDTEM FROM 12 TO 14
C		- CHANGED SSTRNG AND FORMAT STATEMENTS 2012 & 2005 FROM CHAR*133
C		  TO CHAR*128 TO AVOID NULL CHARACTERS IN OUTPUT FILE
C		- ADDED CALL TO XXDATE
C		- CHANGED 'CALL BBSPF1(...,*200)' TO 'CALL BBSPF1(...,*201)'
C		  AND ADDED '201 CONTINUE'
C
C VERSION: 1.6				DATE: 01/05/97
C MODIFIED: KEITH NORMAN
C               - COMMENTED OUT ORIGINAL LINE 909
C                 IF(LEN(SSTRNG).GT.10*NDTEM) SSTRNG(10*NDTEM+1:LEN(SSTRNG))
C                      &           =' '
C                 SINCE PARAMETER (NDTEM=14) AND SSTRNG CHARACTER*128,
C                 HENCE 10*NDTEM+1 > LEN(SSTRNG)
C                 (GENERATES COMPILE TIME ERROR IF NOT COMMENTED OUT)
C
C VERSION: 1.6				DATE: 21/07/97
C MODIFIED: M O'MULLANE & RICHARD MARTIN
C              - CHANGE STRING & SSTRNG TO 151 AND FORMAT STATEMENTS
C		     2000, 2008 AND 2009 TO COPE WITH 14 TEMPERATURES.
C
C VERSION: 1.7				DATE: 09/11/97
C MODIFIED: HUGH SUMMERS & KEITH NORMAN
C              - ARRANGED TRAP FOR BARE NUCLEUS PARENT CONFIGURATION
C                NOT BEING IDENTIFIED AS EISSNER OR STANDARD TYPE
C                ALTERED XMAX AND JH TO CORRECT FAILURE FOR HELIKE C+5
C                ON THE N=5 SHELL. (HPS)
C
C              - ARRANGED SWITCH TO HYDROGENIC GAUNT FACTORS FOR N>NHCUT (HPS)
C
C              - COMMENTED OUT ORIGINAL LINE 909 (KN)
C                IF(LEN(SSTRNG).GT.10*NDTEM) SSTRNG(10*NDTEM+1:LEN(SSTRNG))
C              	    & 	    =' '
C                SINCE PARAMETER (NDTEM=14) AND SSTRNG CHARACTER*128,
C                HENCE 10*NDTEM+1 > LEN(SSTRNG)
C                (GENERATES COMPILE TIME ERROR IF NOT COMMENTED OUT)
C
C VERSION : 1.8				
C DATE    : 02/05/2001
C MODIFIED: Martin O'Mullane
C              - Added a trap for case where term energy is above
C                the ionisation potential. Causes certain binaries
C                to enter an infinite loop.		
C  		18/12/2001
C              - Remove junk at end of adf04 insert file.
C              - Add real name of producer via XXNAME.
C
C VERSION : 1.9				
C DATE    : 21/08/2002
C MODIFIED: Martin O'Mullane
C              - When calculating recombining ion charge from 
C                isoelectronic sequence identifer in adf08 file
C                convert it to uppercase before determining its 
C                atomic number.
C
C VERSION : 1.10				
C DATE    : 31/01/2005
C MODIFIED: Paul Bryans
C              - Non-Maxwellian parameters gathered from IDL input
C              - Calls XX_DATA37 for numerical distributions
C              - BBITRP interpolates distribution to Gaunt factor energies
C              - APHOTDW returns the Gaunt factors calculated therein
C              - Calculate radiative recombination coefficient for
C                non-Maxwellians (see BBRINT)
C              - Spline coefficient to temperatures in adf08 file
C
C VERSION : 1.11				
C DATE    : 02/02/2005
C MODIFIED: Allan Whiteford
C              - Changed call to bbspf0 so that it returns the
C                non-Maxwellian parameters.
C
C VERSION : 1.12
C DATE    : 04-09-2009
C MODIFIED: Martin O'Mullane
C              - Send signal to IDL to indicate termination.
C
C VERSION : 1.13
C DATE    : 09-04-2010
C MODIFIED: Martin O'Mullane
C           - Change integer*4 to integer.
C
C Version : 1.14
C Date    : 23-08-2010
C Modified: Martin O'Mullane
C           - Extra arguments for xxdtes.
C
C Version : 1.15
C Date    : 13-01-2011
C Modified: Martin O'Mullane
C           - Ensure the '---- ---' line is written after the TE line.
C
C Version : 1.16
C Date    : 21-08-2017
C Modified: Martin O'Mullane
C           - Increase number of levels to 200.
C
C-----------------------------------------------------------------------
       INTEGER   NDMET      , NDLEV       , NDTEM       , LOGIC
       INTEGER   IUNT05     , IUNT17      , IUNT09      , PIPEIN    ,
     &           IUNT10     , IUNT11      , IUNT12      , IUNT20    ,
     &           PIPEOU     , ONE         , iunt37
       INTEGER   NHCUT      , dist        , ndgnt       , nemax     ,
     &           ntmax      , icateg      , nenerg      , nblock    ,
     &           nform1     , nform2
C-----------------------------------------------------------------------
       REAL*8    QDMIN
C-----------------------------------------------------------------------
       PARAMETER ( NDMET  = 4  , NDLEV  = 200 , NDTEM  = 14 )
       PARAMETER ( NDGNT  = 5000 ) 
       PARAMETER ( IUNT05 = 5  , IUNT17 = 17  , IUNT09 = 9  )
       PARAMETER ( IUNT10 = 10 , IUNT11 = 11  , IUNT12 = 12 )
       PARAMETER ( IUNT20 = 20 , PIPEIN = 5   , PIPEOU = 6  , ONE = 1)
       PARAMETER ( NHCUT  = 5  , iunt37 = 37  )
       PARAMETER ( QDMIN  = 1.0D-3 )
       parameter ( nemax  = 200, ntmax  = 50  )
C-----------------------------------------------------------------------
       INTEGER   I4UNIT
       INTEGER   LEN3       , LEN4       ,
     &           IZ0        , IZ1         , NPRNTI     , I          ,
     &           IP         , ITRM        , NTRM       , IT         ,
     &           MAXT       , IRESOL      , IBSOPT     , NSHELL     ,
     &           ISHL       , NL          , LP         , ISP        ,
     &           LVCT       , ISCORE      , N          , L          ,
     &           IPP        , IS          , LT         , NMAX       ,
     &           LAM        , JH          , JEALFA     , JSN        ,
     &           IFIRST     , IGONE       , IONCE      , IEXT       ,
     &           IREPT      , LMIN        , LMAX       , L1         ,
     &           LT1        , IWARN       , NVLCE      , maxe       ,
     &           lvlce
C-----------------------------------------------------------------------
       REAL*8    Z0         , Z1          , ACC        , XMAX       ,
     &           H          , U1          , U2         , U3         ,
     &           ZETA       , BWNR        , BWNO       , FRACPR     ,
     &           EN         , V           , QD         , TE         ,
     &           B          , B1          , R1         , ANS        ,
     &           PREC       , PION        , PSTIM      , dparam     ,
     &           rrcint     , rrcin(ntmax), rrcspl(ndtem)           ,
     &           param1     , param2(2)   , ea(ntmax,nemax)         ,
     &           mode(ntmax), median(ntmax),fa(ntmax,nemax)         ,
     &           teff(ntmax), eout(ndgnt) , fout(ndgnt),
     &           ein(nemax) , fin(nemax)  , tin(ntmax)
C-----------------------------------------------------------------------
       CHARACTER STRING*151 , SSTRNG1*50  , SSTRNG*151 , BLANKS*8
       CHARACTER SEQ*2      , SEQSS*100   , SHLSS*30   , SLSS*5
       CHARACTER DSNIN*80   , DSNOUT*80   , DSNPAS*80  , DSN80*80
       CHARACTER DATE*8     , REP*3       , CADAS*120  , BLNK20*20  ,
     &           DSNPAP*80  , USER*30     , file37*120 , titl37*80
       CHARACTER TITLE*40   , filout*120  , calgeb(ntmax,4)*25
       CHARACTER CSTR21*21  , CSTR19*19   , CSTRGO*19
       CHARACTER SEQUP*2    , filnam(nemax)*120
       character ealgeb(ntmax)*25
       character cstr_top*19              , cstr_tail*99
C-----------------------------------------------------------------------
       LOGICAL   OPEN09     , OPEN10      , OPEN11     , open37
       LOGICAL   LCLSHL     , OPEN17
       LOGICAL   LCONT      , LPASS       , LPAPER
       LOGICAL   LBATCH     , LBTSEL      , LPEND
       LOGICAL   LTERM      , LLEVL       , LPTERM     , LPLEVL
       LOGICAL   LEISS      , LSTAN       , lexist     , lprnt      ,
     &           lbndl
C-----------------------------------------------------------------------
       INTEGER   ISPA(NDMET)     , LTPA(NDMET)     , ICNFPA(15,NDMET)
       INTEGER   ISA(NDLEV)      , LTA(NDLEV)      , ICNFA(15,NDLEV)
       INTEGER   NPTA(NDLEV)     , IPLA(NDMET,NDLEV)
       INTEGER   NLQS(10)        , NAA(2)          , LAA(2)
       INTEGER   J0AG(6)         , JCAG(6)
       INTEGER   J0AX(6)         , NA(3)
C-----------------------------------------------------------------------
       REAL*8    WPA(NDMET)      , WNOPA(NDMET)
       REAL*8    WA(NDLEV)       , WNOA(NDLEV)
       REAL*8    FPLA(NDMET,NDLEV)
       REAL*8    EAA(2)          , QDAA(2)         , ALFAA(2,10)
       REAL*8    REAX(17)        , RVAX(23)        , XA(17,23,6)
       REAL*8    EAG(16)         , VAG(22)         , GA(16,22,6)
       REAL*8    TEA(NDTEM)      , EWNA(3)         , ALF(NDTEM)     ,
     &           SUM(NDTEM)      , gaunt(ndgnt)    , vve(ndgnt)
C-----------------------------------------------------------------------
       CHARACTER CNFPA(NDMET)*20 , CNFA(NDLEV)*20 , CPLA(NDLEV)*1
C-----------------------------------------------------------------------
       COMMON /DWPARS/Z0,EAA,QDAA,ALFAA,ACC,XMAX,H,NLQS,NSHELL,NAA,LAA,
     &                JSN,JEALFA,IONCE
       COMMON /BSPARS/Z1,U1,U2,U3,ZETA,IBSOPT,IWARN
       COMMON /PCHGTB/GA,EAG,VAG,J0AG,JCAG,IGONE
       COMMON /PCHXTB/XA,REAX,RVAX,J0AX,IFIRST
C-----------------------------------------------------------------------
       DATA   SEQSS/'H HELIBEB C N O F NENAMGALSIP S CLARK CASCTIV CRMNF
     &ECONICUZNGAGEASSEBRKRRBSRY ZRNBMOTCRURHPDAGCDINSN'/
       DATA   SHLSS/'102021303132404142435051525354'/
       DATA   SLSS /'SPDFG'/
       DATA   BLANKS/'        '/
       DATA   BLNK20/'                    '/
       DATA   OPEN09/.FALSE./ , OPEN10/.FALSE./ , OPEN11/.FALSE./
       DATA   OPEN17/.FALSE./ , open37/.false./
       DATA   LCLSHL/.TRUE./  , LBTSEL/.FALSE./
       DATA   REP   / 'NO' /  , CADAS / ' ' /
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C *************************** MAIN PROGRAM ****************************
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C SET MACHINE DEPENDANT ADAS CONFIGURATION VALUES
C-----------------------------------------------------------------------
C
       CALL XX0000
       CALL XXNAME(USER)
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C GET PARAMETER LBATCH (IS THIS A BATCH RUN OR NOT?) FROM IDL OR FILE
C-----------------------------------------------------------------------
C
      READ( PIPEIN, *) LOGIC
      IF (LOGIC.EQ.0) THEN
         LBATCH	 = .FALSE.
      ELSE
         LBATCH = .TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C GET CURRENT DATE.
C-----------------------------------------------------------------------
C
      CALL XXDATE( DATE )
C
C-----------------------------------------------------------------------
C IF FILE IS ACTIVE ON UNIT 09 - CLOSE THE UNIT
C-----------------------------------------------------------------------
C
 50    IF (OPEN09) THEN
          CLOSE(9)
          OPEN09 = .FALSE.
       ENDIF
C
C-----------------------------------------------------------------------
C GET INPUT DATA SET NAME (FROM IDL OR FILE IF BATCH).
C-----------------------------------------------------------------------
C

       CALL BBSPF0( REP , DSNIN, DIST, DPARAM, FILE37)

C       IF(.NOT.LBATCH)CALL CHECK_PIPE()
C
C-----------------------------------------------------------------------
C IF PROGRAM RUN IS COMPLETED:  END RUN
C-----------------------------------------------------------------------
C
         IF (REP.EQ.'YES') THEN
            GOTO 9999
         ENDIF
C
C-----------------------------------------------------------------------
C OPEN INPUT DATA FILE - DSNIN
C-----------------------------------------------------------------------
C
      OPEN( UNIT=IUNT09 , FILE=DSNIN  , STATUS='UNKNOWN' )
      OPEN09=.TRUE.
C
C-----------------------------------------------------------------------
C PROCESS NON-MAXWELLIAN PARAMETERS
C-----------------------------------------------------------------------
C      
      if (dist.eq.2) then
        lexist = .false.
        if (open37) close(37)
        open37=.false.
        inquire( file=file37 , exist=lexist)
        if(lexist) then
	  open( unit=iunt37 , file=file37 , status='unknown')
          open37=.true.
	else
	  write(i4unit(-1),3040) file37
          stop
        endif
	call xxdata_37( iunt37 ,
     &                  nemax  , ntmax  ,
     &                  titl37 , icateg , nenerg , nblock ,
     &                  nform1 , param1 , nform2 , param2 ,
     &		        ea     , fa     , teff	 , mode   ,
     &		        median , filnam , filout , calgeb ,
     &                  ealgeb
     &                 )
      endif

C
C-----------------------------------------------------------------------
C GET USER-DEFINED VALUES (FROM IDL OR FILE IF BATCH).
C-----------------------------------------------------------------------
C
 55   LPEND = .FALSE.
C
CX      CALL BBISPF( IUNT20 , IPAN   ,
CX     &     LPEND  , LBTSEL , USERID , TITLE
CX     &     )
      CALL BBSPF1( OPEN10 , OPEN11  ,
     &     LCONT  , LPASS  , LPAPER , LBTSEL ,
     &     DSNOUT , DSNPAS , DSNPAP , LPEND, TITLE, *201
     &     )

C       IF(.NOT.LBATCH) CALL CHECK_PIPE()
C
C-----------------------------------------------------------------------
C IF PROGRAM RUN COMPLETE THEN RETURN TO INPUT DATA SET SELECTION PANEL.
C-----------------------------------------------------------------------
C
      IF (LPEND) GOTO 50

C-----------------------------------------------------------------------
C IF THE JOB HAS BEEN SELECTED TO RUN IN BATCH THEN GO BACK TO THE
C INPUT OPTIONS
C-----------------------------------------------------------------------
C
      IF (LBTSEL) GOTO 50

C-----------------------------------------------------------------------
C GET OUTPUT FILE NAMES AND OPEN THEM
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C OPEN PAPER.TXT FILE - IF REQUESTED
C-----------------------------------------------------------------------
C
      IF ((.NOT.OPEN17) .AND. LPAPER) THEN
         DSN80=' '
         DSN80=DSNPAP
         OPEN(UNIT=IUNT17,FILE=DSN80, STATUS='UNKNOWN')
         OPEN17=.TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C OPEN ADF08 FILE OUTPUT DATA SET - DSNOUT - IF REQUESTED
C-----------------------------------------------------------------------
C
      IF ((.NOT.OPEN10) .AND. LCONT) THEN
         DSN80=' '
         DSN80=DSNOUT
         OPEN(UNIT=IUNT10,FILE=DSN80, STATUS='UNKNOWN')
         OPEN10=.TRUE.
      ENDIF
C
C-----------------------------------------------------------------------
C OPEN ADF04 FILE OUTPUT DATA SET - DSNPAS - IF REQUESTED
C-----------------------------------------------------------------------
C
      IF ((.NOT.OPEN11) .AND. LPASS) THEN
         DSN80=' '
         DSN80=DSNPAS
         OPEN(UNIT=IUNT11,FILE=DSN80, STATUS='UNKNOWN')
         OPEN11=.TRUE.
      ENDIF

C
C-----------------------------------------------------------------------
C  FETCH INPUT DATA FROM IUNT09
C-----------------------------------------------------------------------
C
      READ(IUNT09,2000)STRING
      CALL XXSLEN(STRING,LEN3,LEN4)
      WRITE(IUNT10,'(a)')STRING(1:LEN4)
      IF(STRING(1:3).NE.'SEQ')THEN
         WRITE(I4UNIT(-1),*)'*** ERROR - FORMATTING IN LINE 1'
         STOP
      ELSE
         READ(STRING,2001)SEQ,IZ0
         SEQUP = SEQ
         CALL XXSTUC(SEQUP)
         IZ1=IZ0+1-(INDEX(SEQSS,SEQUP)+1)/2
         Z0=IZ0
         Z1=IZ1
C         WRITE(I4UNIT(-1),2201) SEQ, IZ0 ,IZ1
	 IF(OPEN17)THEN
            WRITE(IUNT17,2201) SEQ, IZ0 ,IZ1
         ENDIF
      ENDIF
C-------------------------------------------
C     SET SWITCH FOR CLOSED SHELL ASSUMPTION
C-------------------------------------------
      IF( IZ0 .EQ. IZ1 ) LCLSHL = .FALSE.
C-------------------------------------------
      READ(IUNT09,2000)STRING
      CALL XXSLEN(STRING,LEN3,LEN4)
      WRITE(IUNT10,'(a)')STRING(1:len4)
      READ(IUNT09,2000)STRING
      CALL XXSLEN(STRING,LEN3,LEN4)
      WRITE(IUNT10,'(a)')STRING(1:len4)
C-----------------------------------------------------------------------
C  DETECT IF PARENTS ARE TERMS OR LEVELS
C-----------------------------------------------------------------------
      LPTERM = .FALSE.
      LPLEVL = .FALSE.
      IF(STRING(4:15).EQ.'PARENT TERM ') THEN
          LPTERM = .TRUE.
          READ(STRING,2002) NPRNTI
      ELSEIF(STRING(4:15).EQ.'PARENT LEVEL') THEN
          LPLEVL = .TRUE.
          READ(STRING,2002) NPRNTI
      ELSE
          WRITE(I4UNIT(-1),*)'*** ERROR - PARENT COUPLING INDETERMINATE'
          STOP
      ENDIF
C-----------------------------------------------------------------------
C  READ PARENTS, CONVERT TO STANDARD IF EISSNER CONFIG. FORM
C-----------------------------------------------------------------------
      DO 90 I=1,3
         READ(IUNT09,2000)STRING
         CALL XXSLEN(STRING,LEN3,LEN4)
         WRITE(IUNT10,'(A)')STRING(1:LEN4)
 90   CONTINUE
      DO 100 IP=1,NPRNTI
         READ(IUNT09,2000)STRING
         CALL XXSLEN(STRING,LEN3,LEN4)
         WRITE(IUNT10,'(A)')STRING(1:LEN4)
         READ(STRING,2003)CSTR21,ISPA(IP),LTPA(IP),WPA(IP),WNOPA(IP)
C
         IF (IZ0.EQ.IZ1) THEN
             CNFPA(IP) = BLNK20
             DO 95 I=1,15
               ICNFPA(I,IP) = 0
   95        CONTINUE
         ELSE
             LEISS = .FALSE.
             LSTAN = .FALSE.
             CSTR19 = CSTR21(1:19)
             CALL XXDTES( CSTR19    , LEISS     , LSTAN   , 
     &                    lbndl     , lprnt     ,
     &                    cstr_top  , cstr_tail ,
     &                    nvlce     , lvlce     )
             IF (LEISS) THEN
                 CALL XXCFTR(3,CSTR19,CSTRGO)
                 CNFPA(IP) = CSTRGO(2:19)//CSTR21(21:21)
             ELSEIF (LSTAN) THEN
                 CNFPA(IP) = CSTR21(2:21)
             ELSE
                 WRITE(I4UNIT(-1),*)'*** ERROR - PARENT CONFIG. NOT ',
     &                              'RECOGNISED'
                 STOP
             ENDIF
C            WRITE(I4UNIT(-1),*)'IP=',IP,'  CNFPA(IP)=',CNFPA(IP),
C     &                  '  ISPA(IP)=',ISPA(IP),'  LTPA(IP)=',LTPA(IP),
C     &                  '  WPA(IP)=',WPA(IP),
C     &                  '  WNOPA(IP)=',WNOPA(IP)
             CALL BBPRS3(CNFPA(IP), ICNFPA(1,IP), LCLSHL)
         ENDIF
C        WRITE(I4UNIT(-1),2210) IP , (ICNFPA(I,IP),I=1,15)
         IF(OPEN17)THEN
            WRITE(IUNT17,2210) IP , (ICNFPA(I,IP),I=1,15)
         ENDIF
 100  CONTINUE
      READ(IUNT09,2000)STRING
      CALL XXSLEN(STRING,LEN3,LEN4)
      WRITE(IUNT10,'(A)')STRING(1:LEN4)
      READ(IUNT09,2000)STRING
      CALL XXSLEN(STRING,LEN3,LEN4)
      WRITE(IUNT10,'(A)')STRING(1:LEN4)
C-----------------------------------------------------------------------
C  DETECT IF ENERGY LEVELS ARE TERMS OR LEVELS
C-----------------------------------------------------------------------
      LTERM = .FALSE.
      LLEVL = .FALSE.
      IF(STRING(4:19).EQ.'LS RESOLVED TERM') THEN
          LTERM = .TRUE.
          READ(STRING,2004)BWNR,NTRM
C         WRITE(I4UNIT(-1),*)'BWNR=',BWNR,'  NTRM=',NTRM
      ELSEIF(STRING(4:19).EQ.'J RESOLVED LEVEL') THEN
          LLEVL = .TRUE.
          READ(STRING,2004)BWNR,NTRM
C         WRITE(I4UNIT(-1),*)'BWNR=',BWNR,'  NTRM=',NTRM
      ELSE
          WRITE(I4UNIT(-1),*)'*** ERROR - LEVEL COUPLING INDETERMINATE'
          STOP
      ENDIF
C-----------------------------------------------------------------------
C  READ LEVELS, CONVERT TO STANDARD IF EISSNER CONFIG. FORM
C-----------------------------------------------------------------------
      DO 110 I=1,3
         READ(IUNT09,2000)STRING
         CALL XXSLEN(STRING,LEN3,LEN4)
         WRITE(IUNT10,'(A)')STRING(1:LEN4)
 110  CONTINUE
      DO 120 ITRM=1,NTRM
         READ(IUNT09,2000)STRING
         CALL XXSLEN(STRING,LEN3,LEN4)
         WRITE(IUNT10,'(a)')STRING(1:len4)
         READ(STRING,2006) CSTR21,ISA(ITRM),LTA(ITRM),
     &        WA(ITRM),SSTRNG1
C
         LEISS = .FALSE.
         LSTAN = .FALSE.
         CSTR19 = CSTR21(1:19)
         CALL XXDTES( CSTR19    , LEISS     , LSTAN   , 
     &                lbndl     , lprnt     ,
     &                cstr_top  , cstr_tail ,
     &                nvlce     , lvlce     )
         IF (LEISS) THEN
             CALL XXCFTR(3,CSTR19,CSTRGO)
             CNFA(ITRM) = CSTRGO(2:19)//CSTR21(21:21)
         ELSEIF (LSTAN) THEN
             CNFA(ITRM) = CSTR21(2:21)
         ELSE
             WRITE(I4UNIT(-1),*)'*** ERROR - ENERGY  LEVEL CONFIG. ',
     &                          'NOT RECOGNISED'
             STOP
         ENDIF
C
C         WRITE(I4UNIT(-1),2206) ITRM, CNFA(ITRM), ISA(ITRM),
C     &        LTA(ITRM), WA(ITRM)
         IF(OPEN17)THEN
            WRITE(IUNT17,2206) ITRM, CNFA(ITRM), ISA(ITRM), LTA(ITRM),
     &           WA(ITRM)
         ENDIF
C
         CALL BBPRS1(NDMET,SSTRNG1,WNOA(ITRM),CPLA(ITRM),NPTA(ITRM),
     &        IPLA(1,ITRM),FPLA(1,ITRM))
C
C         WRITE(I4UNIT(-1),*)'WNOA(ITRM)=',WNOA(ITRM),
C     &        '  CPLA(ITRM)=',CPLA(ITRM),
C     &        '  IPLA(IP,ITRM)=',(IPLA(IP,ITRM),IP=1,NPRNTI),
C     &        '  FPLA(IP,ITRM)=',(FPLA(IP,ITRM),IP=1,NPRNTI)
         CALL BBPRS3( CNFA(ITRM), ICNFA(1,ITRM) , LCLSHL )
C         WRITE(I4UNIT(-1),2210) ITRM, (ICNFA(I,ITRM),I=1,15)
         IF(OPEN17)THEN
            WRITE(IUNT17,2210) ITRM, (ICNFA(I,ITRM),I=1,15)
         ENDIF
 120  CONTINUE
      READ(IUNT09,2000)STRING
      CALL XXSLEN(STRING,LEN3,LEN4)
      WRITE(IUNT10,'(A)')STRING(1:LEN4)
 125  READ(IUNT09,2000)STRING
      IF(STRING(9:10).NE.'TE')THEN
         GO TO 125
      ELSE
         READ(STRING,2005)(TEA(IT),IT=1,NDTEM)
         MAXT=NDTEM
         DO 130 IT=NDTEM,1,-1
            IF(TEA(IT).LE.0.0D0)THEN
               MAXT=MAXT-1
            ENDIF
 130     CONTINUE
C         WRITE(I4UNIT(-1),*)(TEA(IT),IT=1,MAXT)
         IF(MAXT.LE.0)THEN
            WRITE(I4UNIT(-1),*)'*** ERROR - NO TEMPERATURES SET'
            STOP
         ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C TELL THE IDL HOW MUCH COMPUTING THERE IS TO DO
C-----------------------------------------------------------------------
      IF(.NOT.LBATCH) THEN
         WRITE(PIPEOU,*)NPRNTI*NTRM
         CALL XXFLSH(PIPEOU)
      ENDIF

C
C-----------------------------------------------------------------------
C     SET UP PARAMETERS AND LOOPS FOR MAIN CALCULATION
C     NB. MUST SET IFIRST AND IGONE INITIALLY
C-----------------------------------------------------------------------
C
      IRESOL=2
C-----------------------------------------------------------------------
C
      DO 200 IP=1,NPRNTI
C
         WRITE(IUNT10,2007) IP,ISPA(IP),SLSS(LTPA(IP)+1:LTPA(IP)+1),
     &        ISPA(IP)
         WRITE(IUNT10,2008) (TEA(IT),IT=1,MAXT)
         WRITE(IUNT10,'(A)')'   ---- ---'
C
         LVCT=0
         DO 15 IT=1,MAXT
 15         SUM(IT)=0.0D0
C
         BWNO=BWNR+WNOPA(IP)
C         WRITE(I4UNIT(-1),2200) IP, BWNO
         NSHELL= 0
         DO 140 ISHL=1,15
            IF(ICNFPA(ISHL,IP).GT.0)THEN
               NSHELL=NSHELL+1
               READ(SHLSS(2*ISHL-1:2*ISHL),'(I2)')NL
               NLQS(NSHELL)=100*NL+ICNFPA(ISHL,IP)
               ALFAA(1,NSHELL)=1.0D0
               ALFAA(2,NSHELL)=1.0D0
            ENDIF
 140     CONTINUE
C         WRITE(I4UNIT(-1),*)'NLQS=',(NLQS(I),I=1,NSHELL)
C         WRITE(I4UNIT(-1),*)'ALFAA(1)=',(ALFAA(1,I),I=1,NSHELL)
C         WRITE(I4UNIT(-1),*)'ALFAA(2)=',(ALFAA(2,I),I=1,NSHELL)
C         WRITE(I4UNIT(-1),*)'IRESOL=',IRESOL,'  IBSOPT=',IBSOPT
C
         LP=LTPA(IP)
         ISP=ISPA(IP)
C-----------------------------------------------------------------------
         DO 190 ITRM=1,NTRM
C
            LVCT=LVCT+1
            DO 22 IT=1,MAXT
 22            ALF(IT)=0.0D0
C
            ISCORE=0
            DO 145 I=1,15
               IF(ICNFA(I,ITRM).NE.ICNFPA(I,IP))THEN
                  ISCORE=ISCORE+1
                  IF(ISCORE.GT.1.OR.(ICNFA(I,ITRM).NE.ICNFPA(I,IP)+1))
     &                 GO TO 189
                  READ(SHLSS(2*I-1:2*I-1),'(I1)')N
                  READ(SHLSS(2*I:2*I),'(I1)')L
               ENDIF
 145        CONTINUE
            IF(ISCORE.NE.1) GO TO 189
            DO 150 IPP=1,NPTA(ITRM)
               IF(IPLA(IPP,ITRM).EQ.IP)THEN
                  FRACPR=FPLA(IPP,ITRM)
                  GO TO 155
               ENDIF
 150        CONTINUE
            GO TO 189
 155        CONTINUE
C            WRITE(I4UNIT(-1),2202) IP, ITRM, N, L, FRACPR
            IS=ISA(ITRM)
            LT=LTA(ITRM)
            NMAX=1
            NA(1)=N
            EN=N
            EWNA(1)=WNOA(ITRM)
C            WRITE(I4UNIT(-1),*)'BWNO=',BWNO,'  EWNA(1)=',EWNA(1)

            IF (BWNO.LT.EWNA(1)) THEN
               WRITE(I4UNIT(-1),2240)IP, ITRM
               STOP
            ENDIF

            V=3.31266D2*Z1/DSQRT(BWNO-EWNA(1))
C
            QD=EN-V
C----------------------------------------------------------------------
C  SET HYDROGENIC CASE IF QUANTUM DEFECT IS BELOW THRESHOLD OR
C  N > NHCUT
C----------------------------------------------------------------------
            IF ((QD.LT.QDMIN). OR.(N.GT.NHCUT)) THEN
               IBSOPT = 3
               WRITE(I4UNIT(-1),2220) ITRM
               WRITE(I4UNIT(-1),2221)
            ELSE
               IBSOPT = 4
            ENDIF
C           WRITE(I4UNIT(-1),*)'Z1=',Z1,'  V=',V,' QD=',QD,
C    &                          ' IBSOPT=',IBSOPT
            EAA(1)=-Z1*Z1/(V*V)
            LAA(1)=L
            NAA(1)=N
            QDAA(1)=QD
            EAA(2)=0.2D0
            LAA(2)=LAA(1)
            NAA(2)=0
            QDAA(2)=QDAA(1)
            LAM=0
C     *******   ALTER MULTIPLIER ON NEXT LINE FOR  TESTS, ORIGINALLY 5.0
C     SET XMAX = 17 FOR HYDROGENIC IONS,  WJD 7/10/92
C-----------------------------------------------------------------------
C     XMAX=17.0*V*V/Z1
C-----------------------------------------------------------------------
C     *******   ALTER MULTIPLIER ON NEXT LINE FOR  CR AND MO RUNS
C     HPS 16/06/95
C-----------------------------------------------------------------------
c           XMAX=40.0*V*V/Z1  old values up to date 05/09/97
C           JH=512
C-----------------------------------------------------------------------
            XMAX=20.0D0*V*V/Z1
            JH=900
            H=JH
            H=XMAX/H
            JEALFA=1
            JSN=0
            ACC=0.01D0
            IFIRST=1
            IGONE=1
            IONCE=0
            IEXT=0
            IREPT=0
C
            IF (IBSOPT.EQ.4)THEN
                CALL ADWLPOL(Z0,NLQS,NSHELL,NAA,LAA,EAA,QDAA,ALFAA,
     &             JSN,JEALFA,ACC,XMAX,H,LAM,IREPT,IEXT,ANS,OPEN17)
C                CALL BBLPOL( Z0,NLQS,NSHELL,NAA,LAA,EAA,QDAA,ALFAA,
C     &            JSN,JEALFA,ACC,XMAX,H,LAM,IREPT,IEXT,ANS)
C
                IFIRST=0
                IGONE=0
                IONCE=1
                JEALFA=0
                DO 157 I=1,NSHELL
                   ALFAA(2,I)=ALFAA(1,I)
 157            CONTINUE
            ENDIF
C
            LMIN=L-1
            IF(LMIN.LT.0.OR.IRESOL.EQ.5)LMIN=LMIN+2
            LMAX=L+1
C
	    DO 60 L1=LMIN,LMAX,2
              IF(IRESOL.NE.5)THEN
                 NAA(2)=0
                 LAA(2)=L1
              ENDIF
C
	      maxe=0
	      DO 40 IT=1,MAXT
                 TE=TEA(IT)
                 B=1.5789D5*Z1*Z1/(TE*V*V)
                 B1=B
                 LT1=0
C
************************************************************************
C
C-----------------------------------------------------------------------
C  APHOTDW now returns each Gaunt factor it calculates as a function of
C  vve, and the number of Gaunt/vve pairs
C-----------------------------------------------------------------------
C
		 CALL APHOTDW(B,B1,V,N,L,L1,LP,ISP,LT,LT1,IS,PREC,PION,
     &                PSTIM,IRESOL,ndgnt,gaunt,vve,maxe)
C		  
		 R1=5.19678D-14*Z1*B**1.5D0*PREC*FRACPR
                 ALF(IT)=ALF(IT)+R1
                 SUM(IT)=SUM(IT)+R1
 40           CONTINUE
C
C-----------------------------------------------------------------------
C  sort Gaunt factors into ascending energy order
C-----------------------------------------------------------------------
C 	       
	      if (dist.ne.0) then
	        call xxsort(maxe,vve,gaunt)
	      endif
	      
	      if (dist.eq.1.or.dist.eq.3) then
	        do it=1,maxt
	          te=tea(it)
	          call bbrint( ndgnt  ,
     &		  	       gaunt  , vve    , z1     , v    ,
     &		               maxe   , te     , dparam , dist , 
     &			       fa     , rrcint
     &			     )
		  alf(it)=rrcint*fracpr
	        enddo
	      
	      elseif (dist.eq.2) then
	        do it=1,nblock
	          te=teff(it)*11605.4d0
		  tin(it) = te
C
C-----------------------------------------------------------------------
C  interpolate distribution function to Gaunt factor energy grid
C-----------------------------------------------------------------------
C		  
		  do i=1,maxe
		    eout(i) = (z1/v)**2d0*13.606d0*vve(i)
		  enddo
		  do i=1,nenerg
		    ein(i) = ea(it,i)
		    fin(i) = fa(it,i)
		  enddo
		  call bbitrp( nemax  , ndgnt  ,
     &                         nenerg , maxe   , te     ,
     &		               nform1 , param1 , nform2 , param2 ,
     &			       ein    , fin    , eout   , fout
     &			     )
C
C-----------------------------------------------------------------------
C  calculate radiative recombination coefficient at adf37 temperatures
C-----------------------------------------------------------------------
C     
		  call bbrint( ndgnt  ,
     &		  	       gaunt  , vve    , z1     , v    ,
     &		               maxe   , te     , dparam , dist , 
     &			       fout   , rrcint
     &			     )
		  rrcin(it) = rrcint
C
C-----------------------------------------------------------------------
C  spline to adf08 temperatures
C-----------------------------------------------------------------------
C
		enddo  
		call bbspln( ndtem   , ntmax   ,
     &                       nblock  , maxt    ,
     &                       tin     , tea     , 
     &                       rrcin   , rrcspl
     &                     ) 
	        do it = 1,maxt
		  alf(it)=rrcspl(it)*fracpr
		enddo
	      endif
	      
 60         CONTINUE
            WRITE(IUNT10,2009) ITRM , (ALF(IT),IT=1,MAXT)
            WRITE(STRING,2009) ITRM , (ALF(IT),IT=1,MAXT)

C Following appears to be necessary to avoid junk at end of string

            SSTRNG = ' '
            WRITE(SSTRNG(1:8),2011)'R',ITRM,'  +',IP
            SSTRNG(9:16)=BLANKS
            DO 185 IT=1,NDTEM
               IF(IT.LE.MAXT) THEN
                  SSTRNG(IT*8+9:IT*8+16)=STRING(IT*10+3:IT*10+7)//
     &                 STRING(IT*10+9:IT*10+11)
               ELSE
                  SSTRNG(IT*8+9:IT*8+16)=BLANKS
               ENDIF
 185        CONTINUE
C
CKN            IF(LEN(SSTRNG).GT.10*NDTEM) SSTRNG(10*NDTEM+1:LEN(SSTRNG))
CKN     &           =' '
C     WRITE(I4UNIT(-1),2012) SSTRNG

            IF(OPEN11) then
               CALL XXSLEN(SSTRNG,LEN3,LEN4)
               WRITE(IUNT11,'(A)')SSTRNG(1:LEN4)
	    ENDIF
	    		
C-----------------------------------------------------------------------
C TELL IDL THAT WE HAVE COMPLETED ONE STAGE FOR PROGRESS WINDOW
C-----------------------------------------------------------------------
 189        IF(.NOT.LBATCH) THEN
               WRITE(PIPEOU,*)ONE
               CALL XXFLSH(PIPEOU)
            ENDIF

 190     CONTINUE
 200  CONTINUE
 201  CONTINUE
      WRITE(IUNT10,2010)
C
      IF(OPEN11) THEN
         WRITE(IUNT11,3000)
         WRITE(IUNT11,3005)
         CALL XXSLEN(DSNIN,LEN3,LEN4)
         WRITE(IUNT11,3010) DSNIN(LEN3:LEN4)
         WRITE(IUNT11,3005)
         WRITE(IUNT11,3020)'ADAS211', USER, DATE
      ENDIF
C

C       IF(.NOT.LBATCH) CALL CHECK_PIPE()

C-----------------------------------------------------------------------
C GO BACK TO DATA-SET INPUT SCREEN
C-----------------------------------------------------------------------

       GO TO 50

C
 9999 CONTINUE
      CLOSE( IUNT17 )
      CLOSE( IUNT10 )
      CLOSE( IUNT11 )
      CLOSE( IUNT20 )
      write(pipeou,*)one
      call xxflsh(pipeou)
C
      STOP
C
C-----------------------------------------------------------------------
C
 1000 FORMAT(A5,3I5,F11.0,F10.5)
 1001 FORMAT(1H1,'RADIATIVE RECOMBINATION RATES'/1H0,
     &     A5,5X,'Z1=',I5,5X,'IBSOPT=',I5,5X,'IRESOL=',I5,5X,
     &     'BWNO=',F11.0,5X,'FRAC.PAR.=',F10.5)
 1002 FORMAT(I5)
 1003 FORMAT(1P,7E10.2)
 1004 FORMAT(3A4,6I3,3(F11.0,I3))
 1005 FORMAT(21X,I3,6X,3(F11.0,I3))
 1006 FORMAT(1H0,6X,'/TE(K)',1P,10D12.4/)
 1007 FORMAT(1H ,1A20,1P,10D11.3)
 1008 FORMAT(1H0,6X,'SUM = ',1P,10D12.4)
 1009 FORMAT(I5,1P,5D12.4)
 1010 FORMAT(1P,E12.4,9I5)
 1011 FORMAT(1H0,'TE=',1P,D12.4,3X,'Z1=',0P,F10.5,3X,'B=',1P,D12.4,3X,
     &     'B1=',1P,D12.4)
 1012 FORMAT(1H0,'SCREENING CONFIGURATION FOR DISTORTED WAVES'//1H ,
     &     'Z0 = ',F5.2//1H ,2X,'NLQ',5X,'ALFA(1)',5X,'ALFA(2)')
 1013 FORMAT(1H ,I5,2X,F10.5,2X,F10.5)
 2000 FORMAT(1A151)
 2001 FORMAT(5X,1A2,13X,I2)
 2002 FORMAT(62X,I2)
 2003 FORMAT(10X,1A21,1X,I1,1X,I1,1X,F4.1,1X,F15.1)
 2004 FORMAT(40X,F15.1,6X,I3)
 2005 FORMAT(11X,14E10.2)
 2006 FORMAT(10X,1A21,1X,I1,1X,I1,1X,F4.1,1X,1A50)
 2007 FORMAT(1H /1H /1H ,'---------------------------------'/
     &     1H ,'PRTI=',I2,2X,'TRMPRT= (',I1,1A1,')  SPNPRT=',I2/
     &     1H )
 2008 FORMAT(1H ,'  INDX TE=',1P,14(D10.2,:))
 2009 FORMAT(1H ,I5,5X,1P,14D10.2)
 2010 FORMAT(1H )
 2011 FORMAT(1A1,I3,1A3,I1)
 2012 FORMAT(1A128)
C
C-----------
C
 2200 FORMAT(1H ,'IP=',I3,3X,' BWNO=',F12.1)
 2201 FORMAT(1H ,'SEQ= ',1A3,3X,'IZ0=',I3,3X,'IZ1=',I3)
 2202 FORMAT(1H ,'IP=',I3,3X,'TRM=',I3,3X,'N=',I3,3X,'  L=',I3,
     &     '   FRACPR=',F6.3)
 2206 FORMAT(1H ,'ITRM=',I3,3X,'CNFA= ',1A20,3X,'ISA=',I3,3X,
     &     'LTA=',I3,3X,'WA=',F6.3)
 2210 FORMAT(I4,' : ',15I4)
 2220 FORMAT(1X,31('*'),' ADAS211 WARNING ',30('*')//
     &       1X,'SMALL OR NEGATIVE QUANTUM DEFECT DETECTED, INDEX : '
     &       ,I3)
 2221 FORMAT(/1X,28('*'),' SWITCH TO HYDROGENIC ',28('*'))
 2240 FORMAT(1X,31('*'),' ADAS211 FATAL ERROR ',30('*')//
     &       1X,'ENERGY OF TERM IS ABOVE IONISATION POTENTIAL - ',
     &          'PARENT : ', I3, ' TERM : ', I3)


 3000 FORMAT('C',79('-'))
 3005 FORMAT('C   ')
 3010 FORMAT('C  Generate radiative recombination from : ',/,'C',5x,A)
 3020   FORMAT('C',/,
     &         'C  CODE     : ',1A7/
     &         'C  PRODUCER : ',A30/
     &         'C  DATE     : ',1A8,/,'C',/,'C',79('-'))
 3030  format(5(1e15.8))
 3040  FORMAT(/ /,'FILE DOES NOT EXIST: ',1A80)
C----------------------------------------------------------------------------

      END
