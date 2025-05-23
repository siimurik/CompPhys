DOUBLE PRECISION FUNCTION D1MACH(I)
    INTEGER :: I

!  DOUBLE-PRECISION MACHINE CONSTANTS

!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.

!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.

!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.

!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.

!  D1MACH( 5) = LOG10(B)

!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.
!  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
!  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)

!  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
!  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.

!  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
!  TO SPECIFY THE CONSTANTS EXACTLY, WHICH HAS IN SOME CASES
!  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.

    EXTERNAL I1MACH
    INTEGER ::  I1MACH

    INTEGER :: SMALL(4)
    INTEGER :: LARGE(4)
    INTEGER :: RIGHT(4)
    INTEGER :: DIVER(4)
    INTEGER :: LOG10(4)

    DOUBLE PRECISION :: DMACH(5)

    EQUIVALENCE (DMACH(1),SMALL(1))
    EQUIVALENCE (DMACH(2),LARGE(1))
    EQUIVALENCE (DMACH(3),RIGHT(1))
    EQUIVALENCE (DMACH(4),DIVER(1))
    EQUIVALENCE (DMACH(5),LOG10(1))


!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.

    DATA SMALL(1),SMALL(2) /    1048576,          0 /
    DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
    DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
    DATA DIVER(1),DIVER(2) / 1018167296,          0 /
    DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /

!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
!     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
!     SIGNIFICANT BYTE IS STORED FIRST.

!       DATA SMALL(1),SMALL(2) /          0,    1048576 /
!       DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
!       DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
!       DATA DIVER(1),DIVER(2) /          0, 1018167296 /
!       DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /

!     MACHINE CONSTANTS FOR AMDAHL MACHINES.

!      DATA SMALL(1),SMALL(2) /    1048576,          0 /
!      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
!      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
!      DATA DIVER(1),DIVER(2) /  873463808,          0 /
!      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.

!      DATA SMALL(1) / ZC00800000 /
!      DATA SMALL(2) / Z000000000 /

!      DATA LARGE(1) / ZDFFFFFFFF /
!      DATA LARGE(2) / ZFFFFFFFFF /

!      DATA RIGHT(1) / ZCC5800000 /
!      DATA RIGHT(2) / Z000000000 /

!      DATA DIVER(1) / ZCC6800000 /
!      DATA DIVER(2) / Z000000000 /

!      DATA LOG10(1) / ZD00E730E7 /
!      DATA LOG10(2) / ZC77800DC0 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.

!      DATA SMALL(1) / O1771000000000000 /
!      DATA SMALL(2) / O0000000000000000 /

!      DATA LARGE(1) / O0777777777777777 /
!      DATA LARGE(2) / O0007777777777777 /

!      DATA RIGHT(1) / O1461000000000000 /
!      DATA RIGHT(2) / O0000000000000000 /

!      DATA DIVER(1) / O1451000000000000 /
!      DATA DIVER(2) / O0000000000000000 /

!      DATA LOG10(1) / O1157163034761674 /
!      DATA LOG10(2) / O0006677466732724 /

!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.

!      DATA SMALL(1) / O1771000000000000 /
!      DATA SMALL(2) / O7770000000000000 /

!      DATA LARGE(1) / O0777777777777777 /
!      DATA LARGE(2) / O7777777777777777 /

!      DATA RIGHT(1) / O1461000000000000 /
!      DATA RIGHT(2) / O0000000000000000 /

!      DATA DIVER(1) / O1451000000000000 /
!      DATA DIVER(2) / O0000000000000000 /

!      DATA LOG10(1) / O1157163034761674 /
!      DATA LOG10(2) / O0006677466732724 /

!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.

!      DATA SMALL(1) / 00604000000000000000B /
!      DATA SMALL(2) / 00000000000000000000B /

!      DATA LARGE(1) / 37767777777777777777B /
!      DATA LARGE(2) / 37167777777777777777B /

!      DATA RIGHT(1) / 15604000000000000000B /
!      DATA RIGHT(2) / 15000000000000000000B /

!      DATA DIVER(1) / 15614000000000000000B /
!      DATA DIVER(2) / 15010000000000000000B /

!      DATA LOG10(1) / 17164642023241175717B /
!      DATA LOG10(2) / 16367571421742254654B /


!     MACHINE CONSTANTS FOR CONVEX C-1

!       DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /
!       DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
!       DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /
!       DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /
!       DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /

!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.

!      DATA SMALL(1) / 201354000000000000000B /
!      DATA SMALL(2) / 000000000000000000000B /

!      DATA LARGE(1) / 577767777777777777777B /
!      DATA LARGE(2) / 000007777777777777776B /

!      DATA RIGHT(1) / 376434000000000000000B /
!      DATA RIGHT(2) / 000000000000000000000B /

!      DATA DIVER(1) / 376444000000000000000B /
!      DATA DIVER(2) / 000000000000000000000B /

!      DATA LOG10(1) / 377774642023241175717B /
!      DATA LOG10(2) / 000007571421742254654B /

!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200

!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
!     STATIC DMACH(5)

!      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
!      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
!      DATA LOG10/40423K,42023K,50237K,74776K/

!     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7

!      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
!      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
!      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
!      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
!      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /

!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.

!      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /

!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.

!      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
!      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
!      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
!      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /

!     MACHINE CONSTANTS FOR THE INTERDATA 8/32
!     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.

!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
!     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

!      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
!      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
!      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
!      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
!      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /

!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).

!      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
!      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
!      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
!      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
!      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /

!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).

!      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
!      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
!      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
!      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
!      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /

!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).

!      DATA SMALL(1),SMALL(2) /    8388608,           0 /
!      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
!      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
!      DATA DIVER(1),DIVER(2) /  620756992,           0 /
!      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /

!      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
!      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
!      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
!      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
!      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /

!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).

!      DATA SMALL(1),SMALL(2) /    128,      0 /
!      DATA SMALL(3),SMALL(4) /      0,      0 /

!      DATA LARGE(1),LARGE(2) /  32767,     -1 /
!      DATA LARGE(3),LARGE(4) /     -1,     -1 /

!      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
!      DATA RIGHT(3),RIGHT(4) /      0,      0 /

!      DATA DIVER(1),DIVER(2) /   9472,      0 /
!      DATA DIVER(3),DIVER(4) /      0,      0 /

!      DATA LOG10(1),LOG10(2) /  16282,   8346 /
!      DATA LOG10(3),LOG10(4) / -31493, -12296 /

!      DATA SMALL(1),SMALL(2) / O000200, O000000 /
!      DATA SMALL(3),SMALL(4) / O000000, O000000 /

!      DATA LARGE(1),LARGE(2) / O077777, O177777 /
!      DATA LARGE(3),LARGE(4) / O177777, O177777 /

!      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
!      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /

!      DATA DIVER(1),DIVER(2) / O022400, O000000 /
!      DATA DIVER(3),DIVER(4) / O000000, O000000 /

!      DATA LOG10(1),LOG10(2) / O037632, O020232 /
!      DATA LOG10(3),LOG10(4) / O102373, O147770 /

!     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
!     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
!     SUPPLIED BY IGOR BRAY.

!      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /
!      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /
!      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /
!      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /
!      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /

!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000

!      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
!      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
!      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
!      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
!      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /

!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.

!      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /

!     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER

!       DATA SMALL(1),SMALL(2) /        128,           0 /
!       DATA LARGE(1),LARGE(2) /     -32769,          -1 /
!       DATA RIGHT(1),RIGHT(2) /       9344,           0 /
!       DATA DIVER(1),DIVER(2) /       9472,           0 /
!       DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /

!     MACHINE CONSTANTS FOR THE VAX-11 WITH
!     FORTRAN IV-PLUS COMPILER

!      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /
!      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /
!      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /
!      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /

!     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2

!      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /
!      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
!      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /
!      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /
!      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /

!       THIS WILL CAUSE A COMPILER ERROR


    IF (I < 1  .OR.  I > 5) GOTO 999
    D1MACH = DMACH(I)
    RETURN
    999 WRITE(I1MACH(2),1999) I
    1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10)
    STOP
END FUNCTION

!	Caveat receptor.  (Jack) dongarra@anl-mcs, (Eric Grosse) research!ehg
!	Compliments of netlib   Wed Aug 19 10:08:28 CDT 1987
INTEGER FUNCTION I1MACH(I)
    INTEGER :: I

!  I/O UNIT NUMBERS.

!    I1MACH( 1) = THE STANDARD INPUT UNIT.

!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.

!    I1MACH( 3) = THE STANDARD PUNCH UNIT.

!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.

!  WORDS.

!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.

!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
!                 FOR  FORTRAN 77, THIS IS ALWAYS 1.  FOR FORTRAN 66,
!                 CHARACTER STORAGE UNIT = INTEGER STORAGE UNIT.

!  INTEGERS.

!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM

!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.

!    I1MACH( 7) = A, THE BASE.

!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.

!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.

!  FLOATING-POINT NUMBERS.

!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM

!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )

!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.

!    I1MACH(10) = B, THE BASE.

!  SINGLE-PRECISION

!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.

!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.

!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION

!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.

!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.

!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.

!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
!  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
!  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
!  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
!  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
!  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
!  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)

!  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
!  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
!  FOR IMACH(1) - IMACH(4).

    INTEGER :: IMACH(16),OUTPUT,SANITY

    EQUIVALENCE (IMACH(4),OUTPUT)



!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

    DATA IMACH( 1) /    5 /
    DATA IMACH( 2) /    6 /
    DATA IMACH( 3) /    7 /
    DATA IMACH( 4) /    0 /
    DATA IMACH( 5) /   32 /
    DATA IMACH( 6) /    4 /
    DATA IMACH( 7) /    2 /
    DATA IMACH( 8) /   31 /
    DATA IMACH( 9) / 2147483647 /
    DATA IMACH(10) /    2 /
    DATA IMACH(11) /   24 /
    DATA IMACH(12) / -125 /
    DATA IMACH(13) /  128 /
    DATA IMACH(14) /   53 /
    DATA IMACH(15) / -1021 /
    DATA IMACH(16) /  1024 /, SANITY/987/

!     MACHINE CONSTANTS FOR AMDAHL MACHINES.

!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  63 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  63 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.

!      DATA IMACH( 1) /    7 /
!      DATA IMACH( 2) /    2 /
!      DATA IMACH( 3) /    2 /
!      DATA IMACH( 4) /    2 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   33 /
!      DATA IMACH( 9) / Z1FFFFFFFF /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -256 /
!      DATA IMACH(13) /  255 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) / -256 /
!      DATA IMACH(16) /  255 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.

!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  48 /
!      DATA IMACH( 6) /   6 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  39 /
!      DATA IMACH( 9) / O0007777777777777 /
!      DATA IMACH(10) /   8 /
!      DATA IMACH(11) /  13 /
!      DATA IMACH(12) / -50 /
!      DATA IMACH(13) /  76 /
!      DATA IMACH(14) /  26 /
!      DATA IMACH(15) / -50 /
!      DATA IMACH(16) /  76 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.

!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  48 /
!      DATA IMACH( 6) /   6 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  39 /
!      DATA IMACH( 9) / O0007777777777777 /
!      DATA IMACH(10) /   8 /
!      DATA IMACH(11) /  13 /
!      DATA IMACH(12) / -50 /
!      DATA IMACH(13) /  76 /
!      DATA IMACH(14) /  26 /
!      DATA IMACH(15) / -32754 /
!      DATA IMACH(16) /  32780 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.

!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   60 /
!      DATA IMACH( 6) /   10 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   48 /
!      DATA IMACH( 9) / 00007777777777777777B /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   48 /
!      DATA IMACH(12) / -974 /
!      DATA IMACH(13) / 1070 /
!      DATA IMACH(14) /   96 /
!      DATA IMACH(15) / -927 /
!      DATA IMACH(16) / 1070 /, SANITY/987/


!     MACHINE CONSTANTS FOR CONVEX C-1.

!       DATA IMACH( 1) /    5 /
!       DATA IMACH( 2) /    6 /
!       DATA IMACH( 3) /    7 /
!       DATA IMACH( 4) /    6 /
!       DATA IMACH( 5) /   32 /
!       DATA IMACH( 6) /    4 /
!       DATA IMACH( 7) /    2 /
!       DATA IMACH( 8) /   31 /
!       DATA IMACH( 9) / 2147483647 /
!       DATA IMACH(10) /    2 /
!       DATA IMACH(11) /   24 /
!       DATA IMACH(12) / -128 /
!       DATA IMACH(13) /  127 /
!       DATA IMACH(14) /   53 /
!       DATA IMACH(15) /-1024 /
!       DATA IMACH(16) / 1023 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.

!      DATA IMACH( 1) /     5 /
!      DATA IMACH( 2) /     6 /
!      DATA IMACH( 3) /   102 /
!      DATA IMACH( 4) /     6 /
!      DATA IMACH( 5) /    64 /
!      DATA IMACH( 6) /     8 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    63 /
!      DATA IMACH( 9) /  777777777777777777777B /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    47 /
!      DATA IMACH(12) / -8189 /
!      DATA IMACH(13) /  8190 /
!      DATA IMACH(14) /    94 /
!      DATA IMACH(15) / -8099 /
!      DATA IMACH(16) /  8190 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.

!      DATA IMACH( 1) /   11 /
!      DATA IMACH( 2) /   12 /
!      DATA IMACH( 3) /    8 /
!      DATA IMACH( 4) /   10 /
!      DATA IMACH( 5) /   16 /
!      DATA IMACH( 6) /    2 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   15 /
!      DATA IMACH( 9) /32767 /
!      DATA IMACH(10) /   16 /
!      DATA IMACH(11) /    6 /
!      DATA IMACH(12) /  -64 /
!      DATA IMACH(13) /   63 /
!      DATA IMACH(14) /   14 /
!      DATA IMACH(15) /  -64 /
!      DATA IMACH(16) /   63 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.

!      DATA IMACH( 1) /       5 /
!      DATA IMACH( 2) /       6 /
!      DATA IMACH( 3) /       0 /
!      DATA IMACH( 4) /       6 /
!      DATA IMACH( 5) /      24 /
!      DATA IMACH( 6) /       3 /
!      DATA IMACH( 7) /       2 /
!      DATA IMACH( 8) /      23 /
!      DATA IMACH( 9) / 8388607 /
!      DATA IMACH(10) /       2 /
!      DATA IMACH(11) /      23 /
!      DATA IMACH(12) /    -127 /
!      DATA IMACH(13) /     127 /
!      DATA IMACH(14) /      38 /
!      DATA IMACH(15) /    -127 /
!      DATA IMACH(16) /     127 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.

!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /   43 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   63 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.

!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   7 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / Z7FFFFFFF /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  63 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  63 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE INTERDATA 8/32
!     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.

!     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
!     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

!      DATA IMACH( 1) /   5 /
!      DATA IMACH( 2) /   6 /
!      DATA IMACH( 3) /   6 /
!      DATA IMACH( 4) /   6 /
!      DATA IMACH( 5) /  32 /
!      DATA IMACH( 6) /   4 /
!      DATA IMACH( 7) /   2 /
!      DATA IMACH( 8) /  31 /
!      DATA IMACH( 9) / Z'7FFFFFFF' /
!      DATA IMACH(10) /  16 /
!      DATA IMACH(11) /   6 /
!      DATA IMACH(12) / -64 /
!      DATA IMACH(13) /  62 /
!      DATA IMACH(14) /  14 /
!      DATA IMACH(15) / -64 /
!      DATA IMACH(16) /  62 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).

!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    5 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / "377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   54 /
!      DATA IMACH(15) / -101 /
!      DATA IMACH(16) /  127 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).

!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    5 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / "377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   62 /
!      DATA IMACH(15) / -128 /
!      DATA IMACH(16) /  127 /, SANITY/987/

!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGER ARITHMETIC.

!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /, SANITY/987/

!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     16-BIT INTEGER ARITHMETIC.

!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   16 /
!      DATA IMACH( 6) /    2 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   15 /
!      DATA IMACH( 9) / 32767 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
!     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
!     SUPPLIED BY IGOR BRAY.

!      DATA IMACH( 1) /            1 /
!      DATA IMACH( 2) /            1 /
!      DATA IMACH( 3) /            2 /
!      DATA IMACH( 4) /            1 /
!      DATA IMACH( 5) /           32 /
!      DATA IMACH( 6) /            4 /
!      DATA IMACH( 7) /            2 /
!      DATA IMACH( 8) /           31 /
!      DATA IMACH( 9) / :17777777777 /
!      DATA IMACH(10) /            2 /
!      DATA IMACH(11) /           23 /
!      DATA IMACH(12) /         -127 /
!      DATA IMACH(13) /         +127 /
!      DATA IMACH(14) /           47 /
!      DATA IMACH(15) /       -32895 /
!      DATA IMACH(16) /       +32637 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.

!      DATA IMACH( 1) /     0 /
!      DATA IMACH( 2) /     0 /
!      DATA IMACH( 3) /     7 /
!      DATA IMACH( 4) /     0 /
!      DATA IMACH( 5) /    32 /
!      DATA IMACH( 6) /     1 /
!      DATA IMACH( 7) /     2 /
!      DATA IMACH( 8) /    31 /
!      DATA IMACH( 9) /  2147483647 /
!      DATA IMACH(10) /     2 /
!      DATA IMACH(11) /    24 /
!      DATA IMACH(12) /  -125 /
!      DATA IMACH(13) /   128 /
!      DATA IMACH(14) /    53 /
!      DATA IMACH(15) / -1021 /
!      DATA IMACH(16) /  1024 /, SANITY/987/

!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.

!     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
!     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
!     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.

!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    6 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /, SANITY/987/

!     MACHINE CONSTANTS FOR VAX.

!       DATA IMACH( 1) /    5 /
!       DATA IMACH( 2) /    6 /
!       DATA IMACH( 3) /    7 /
!       DATA IMACH( 4) /    0 /
!       DATA IMACH( 5) /   32 /
!       DATA IMACH( 6) /    4 /
!       DATA IMACH( 7) /    2 /
!       DATA IMACH( 8) /   31 /
!       DATA IMACH( 9) / 2147483647 /
!       DATA IMACH(10) /    2 /
!       DATA IMACH(11) /   24 /
!       DATA IMACH(12) / -127 /
!       DATA IMACH(13) /  127 /
!       DATA IMACH(14) /   56 /
!       DATA IMACH(15) / -127 /
!       DATA IMACH(16) /  127 /, SANITY/987/


!  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
    IF (SANITY /= 987) STOP 777
    IF (I < 1  .OR.  I > 16) GO TO 999
    I1MACH=IMACH(I)
! 6S
! 7S
    IF(I == 6) I1MACH=1
!/
    RETURN
    999 WRITE(OUTPUT,1999) I
    1999 FORMAT(' I1MACH - I OUT OF BOUNDS',I10)
    STOP
END FUNCTION I1MACH


subroutine dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr, &
    neval,ier,alist,blist,rlist,elist,iord,last)
!***begin prologue  dqage
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
!***keywords  automatic integrator, general-purpose,
!             integrand examinator, globally adaptive,
!             gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral   i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
!***description

!        computation of a definite integral
!        standard fortran subroutine
!        double precision version

!        parameters
!         on entry
!            f      - double precision
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.

!            a      - double precision
!                     lower limit of integration

!            b      - double precision
!                     upper limit of integration

!            epsabs - double precision
!                     absolute accuracy requested
!            epsrel - double precision
!                     relative accuracy requested
!                     if  epsabs.le.0
!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.

!            key    - integer
!                     key for choice of local integration rule
!                     a gauss-kronrod pair is used with
!                          7 - 15 points if key.lt.2,
!                         10 - 21 points if key = 2,
!                         15 - 31 points if key = 3,
!                         20 - 41 points if key = 4,
!                         25 - 51 points if key = 5,
!                         30 - 61 points if key.gt.5.

!            limit  - integer
!                     gives an upperbound on the number of subintervals
!                     in the partition of (a,b), limit.ge.1.

!         on return
!            result - double precision
!                     approximation to the integral

!            abserr - double precision
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)

!            neval  - integer
!                     number of integrand evaluations

!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier.gt.0 abnormal termination of the routine
!                             the estimates for result and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value
!                             of limit.
!                             however, if this yields no improvement it
!                             is rather advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined(e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. if possible, an appropriate
!                             special-purpose integrator should be used
!                             which is designed for handling the type of
!                             difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             (epsabs.le.0 and
!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                             result, abserr, neval, last, rlist(1) ,
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.

!            alist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)

!            blist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)

!            rlist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the
!                      integral approximations on the subintervals

!            elist   - double precision
!                      vector of dimension at least limit, the first
!                       last  elements of which are the moduli of the
!                      absolute error estimates on the subintervals

!            iord    - integer
!                      vector of dimension at least limit, the first k
!                      elements of which are pointers to the
!                      error estimates over the subintervals,
!                      such that elist(iord(1)), ...,
!                      elist(iord(k)) form a decreasing sequence,
!                      with k = last if last.le.(limit/2+2), and
!                      k = limit+1-last otherwise

!            last    - integer
!                      number of subintervals actually produced in the
!                      subdivision process

!***references  (none)
!***routines called  d1mach,dqk15,dqk21,dqk31,
!                    dqk41,dqk51,dqk61,dqpsrt
!***end prologue  dqage

    double precision :: a,abserr,alist,area,area1,area12,area2,a1,a2,b, &
    blist,b1,b2,dabs,defabs,defab1,defab2,dmax1,d1mach,elist,epmach, &
    epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f, &
    resabs,result,rlist,uflow
    integer :: ier,iord,iroff1,iroff2,k,key,keyf,last,limit,maxerr,neval, &
    nrmax

    dimension alist(limit),blist(limit),elist(limit),iord(limit), &
    rlist(limit)

    external f

!            list of major variables
!            -----------------------

!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                      (alist(i),blist(i))
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest
!                       error estimate
!           errmax    - elist(maxerr)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision


!           machine dependent constants
!           ---------------------------

!           epmach  is the largest relative spacing.
!           uflow  is the smallest positive magnitude.

!***first executable statement  dqage
    epmach = d1mach(4)
    uflow = d1mach(1)

!           test on validity of parameters
!           ------------------------------

    ier = 0
    neval = 0
    last = 0
    result = 0.0d+00
    abserr = 0.0d+00
    alist(1) = a
    blist(1) = b
    rlist(1) = 0.0d+00
    elist(1) = 0.0d+00
    iord(1) = 0
    if(epsabs <= 0.0d+00 .AND. &
    epsrel < dmax1(0.5d+02*epmach,0.5d-28)) ier = 6
    if(ier == 6) go to 999

!           first approximation to the integral
!           -----------------------------------

    keyf = key
    if(key <= 0) keyf = 1
    if(key >= 7) keyf = 6
    neval = 0
    if(keyf == 1) call dqk15(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 2) call dqk21(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 3) call dqk31(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 4) call dqk41(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 5) call dqk51(f,a,b,result,abserr,defabs,resabs)
    if(keyf == 6) call dqk61(f,a,b,result,abserr,defabs,resabs)
    last = 1
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1

!           test on accuracy.

    errbnd = dmax1(epsabs,epsrel*dabs(result))
    if(abserr <= 0.5d+02*epmach*defabs .AND. abserr > errbnd) ier = 2
    if(limit == 1) ier = 1
    if(ier /= 0 .OR. (abserr <= errbnd .AND. abserr /= resabs) &
     .OR. abserr == 0.0d+00) go to 60

!           initialization
!           --------------


    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    nrmax = 1
    iroff1 = 0
    iroff2 = 0

!           main do-loop
!           ------------

    do 30 last = 2,limit
    
    !           bisect the subinterval with the largest error estimate.
    
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf == 1) call dqk15(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 2) call dqk21(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 3) call dqk31(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 4) call dqk41(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 5) call dqk51(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 6) call dqk61(f,a1,b1,area1,error1,resabs,defab1)
        if(keyf == 1) call dqk15(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 2) call dqk21(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 3) call dqk31(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 4) call dqk41(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 5) call dqk51(f,a2,b2,area2,error2,resabs,defab2)
        if(keyf == 6) call dqk61(f,a2,b2,area2,error2,resabs,defab2)
    
    !           improve previous approximations to integral
    !           and error and test for accuracy.
    
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1 == error1 .OR. defab2 == error2) go to 5
        if(dabs(rlist(maxerr)-area12) <= 0.1d-04*dabs(area12) &
         .AND. erro12 >= 0.99d+00*errmax) iroff1 = iroff1+1
        if(last > 10 .AND. erro12 > errmax) iroff2 = iroff2+1
        5 rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = dmax1(epsabs,epsrel*dabs(area))
        if(errsum <= errbnd) go to 8
    
    !           test for roundoff error and eventually set error flag.
    
        if(iroff1 >= 6 .OR. iroff2 >= 20) ier = 2
    
    !           set error flag in the case that the number of subintervals
    !           equals limit.
    
        if(last == limit) ier = 1
    
    !           set error flag in the case of bad integrand behaviour
    !           at a point of the integration range.
    
        if(dmax1(dabs(a1),dabs(b2)) <= (0.1d+01+0.1d+03* &
        epmach)*(dabs(a2)+0.1d+04*uflow)) ier = 3
    
    !           append the newly-created intervals to the list.
    
        8 if(error2 > error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
        10 alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
    
    !           call subroutine dqpsrt to maintain the descending ordering
    !           in the list of error estimates and select the subinterval
    !           with the largest error estimate (to be bisected next).
    
        20 call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    ! ***jump out of do-loop
        if(ier /= 0 .OR. errsum <= errbnd) go to 40
    30 end do

!           compute final result.
!           ---------------------

    40 result = 0.0d+00
    do 50 k=1,last
        result = result+rlist(k)
    50 end do
    abserr = errsum
    60 if(keyf /= 1) neval = (10*keyf+1)*(2*neval+1)
    if(keyf == 1) neval = 30*neval+15
    999 return
end subroutine dqage


subroutine dqk15(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk15
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description

!           integration rules
!           standard fortran subroutine
!           double precision version

!           parameters
!            on entry
!              f      - double precision
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.

!              a      - double precision
!                       lower limit of integration

!              b      - double precision
!                       upper limit of integration

!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the7-point gauss rule(resg).

!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk15

    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    !external f

    dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule

!           wgk    - weights of the 15-point kronrod rule

!           wg     - weights of the 7-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.129484966168869693270611432679082d0 /
    data wg  (  2) / 0.279705391489276667901467771423780d0 /
    data wg  (  3) / 0.381830050505118944950369775488975d0 /
    data wg  (  4) / 0.417959183673469387755102040816327d0 /

    data xgk (  1) / 0.991455371120812639206854697526329d0 /
    data xgk (  2) / 0.949107912342758524526189684047851d0 /
    data xgk (  3) / 0.864864423359769072789712788640926d0 /
    data xgk (  4) / 0.741531185599394439863864773280788d0 /
    data xgk (  5) / 0.586087235467691130294144838258730d0 /
    data xgk (  6) / 0.405845151377397166906606412076961d0 /
    data xgk (  7) / 0.207784955007898467600689403773245d0 /
    data xgk (  8) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.022935322010529224963732008058970d0 /
    data wgk (  2) / 0.063092092629978553290700663189204d0 /
    data wgk (  3) / 0.104790010322250183839876322541518d0 /
    data wgk (  4) / 0.140653259715525918745189590510238d0 /
    data wgk (  5) / 0.169004726639267902826583426598550d0 /
    data wgk (  6) / 0.190350578064785409913256402421014d0 /
    data wgk (  7) / 0.204432940075298892414161999234649d0 /
    data wgk (  8) / 0.209482141084727828012999174891714d0 /


!           list of major variables
!           -----------------------

!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)

!           machine dependent constants
!           ---------------------------

!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

!***first executable statement  dqk15
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.

    fc = f(centr)
    resg = fc*wg(4)
    resk = fc*wgk(8)
    resabs = dabs(resk)
    do j=1,3
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,4
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5d+00
    resasc = wgk(8)*dabs(fc-reskh)
    do j=1,7
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.0d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
end subroutine dqk15


subroutine dqk21(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk21
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  21-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description

!           integration rules
!           standard fortran subroutine
!           double precision version

!           parameters
!            on entry
!              f      - double precision
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the driver program.

!              a      - double precision
!                       lower limit of integration

!              b      - double precision
!                       upper limit of integration

!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 21-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 10-point gauss rule (resg).

!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk21

    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    !external f

    dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 21-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point gauss rule

!           wgk    - weights of the 21-point kronrod rule

!           wg     - weights of the 10-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.066671344308688137593568809893332d0 /
    data wg  (  2) / 0.149451349150580593145776339657697d0 /
    data wg  (  3) / 0.219086362515982043995534934228163d0 /
    data wg  (  4) / 0.269266719309996355091226921569469d0 /
    data wg  (  5) / 0.295524224714752870173892994651338d0 /

    data xgk (  1) / 0.995657163025808080735527280689003d0 /
    data xgk (  2) / 0.973906528517171720077964012084452d0 /
    data xgk (  3) / 0.930157491355708226001207180059508d0 /
    data xgk (  4) / 0.865063366688984510732096688423493d0 /
    data xgk (  5) / 0.780817726586416897063717578345042d0 /
    data xgk (  6) / 0.679409568299024406234327365114874d0 /
    data xgk (  7) / 0.562757134668604683339000099272694d0 /
    data xgk (  8) / 0.433395394129247190799265943165784d0 /
    data xgk (  9) / 0.294392862701460198131126603103866d0 /
    data xgk ( 10) / 0.148874338981631210884826001129720d0 /
    data xgk ( 11) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.011694638867371874278064396062192d0 /
    data wgk (  2) / 0.032558162307964727478818972459390d0 /
    data wgk (  3) / 0.054755896574351996031381300244580d0 /
    data wgk (  4) / 0.075039674810919952767043140916190d0 /
    data wgk (  5) / 0.093125454583697605535065465083366d0 /
    data wgk (  6) / 0.109387158802297641899210590325805d0 /
    data wgk (  7) / 0.123491976262065851077958109831074d0 /
    data wgk (  8) / 0.134709217311473325928054001771707d0 /
    data wgk (  9) / 0.142775938577060080797094273138717d0 /
    data wgk ( 10) / 0.147739104901338491374841515972068d0 /
    data wgk ( 11) / 0.149445554002916905664936468389821d0 /


!           list of major variables
!           -----------------------

!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point gauss formula
!           resk   - result of the 21-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)


!           machine dependent constants
!           ---------------------------

!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

!***first executable statement  dqk21
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 21-point kronrod approximation to
!           the integral, and estimate the absolute error.

    resg = 0.0d+00
    fc = f(centr)
    resk = wgk(11)*fc
    resabs = dabs(resk)
    do j=1,5
        jtw = 2*j
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,5
        jtwm1 = 2*j-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5d+00
    resasc = wgk(11)*dabs(fc-reskh)
    do j=1,10
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.0d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
end subroutine dqk21


subroutine dqk31(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk31
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  31-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description

!           integration rules
!           standard fortran subroutine
!           double precision version

!           parameters
!            on entry
!              f      - double precision
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.

!              a      - double precision
!                       lower limit of integration

!              b      - double precision
!                       upper limit of integration

!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 31-point
!                       gauss-kronrod rule (resk), obtained by optimal
!                       addition of abscissae to the 15-point gauss
!                       rule (resg).

!              abserr - double precison
!                       estimate of the modulus of the modulus,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk31
    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    !external f

    dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 31-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 15-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 15-point gauss rule

!           wgk    - weights of the 31-point kronrod rule

!           wg     - weights of the 15-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.030753241996117268354628393577204d0 /
    data wg  (  2) / 0.070366047488108124709267416450667d0 /
    data wg  (  3) / 0.107159220467171935011869546685869d0 /
    data wg  (  4) / 0.139570677926154314447804794511028d0 /
    data wg  (  5) / 0.166269205816993933553200860481209d0 /
    data wg  (  6) / 0.186161000015562211026800561866423d0 /
    data wg  (  7) / 0.198431485327111576456118326443839d0 /
    data wg  (  8) / 0.202578241925561272880620199967519d0 /

    data xgk (  1) / 0.998002298693397060285172840152271d0 /
    data xgk (  2) / 0.987992518020485428489565718586613d0 /
    data xgk (  3) / 0.967739075679139134257347978784337d0 /
    data xgk (  4) / 0.937273392400705904307758947710209d0 /
    data xgk (  5) / 0.897264532344081900882509656454496d0 /
    data xgk (  6) / 0.848206583410427216200648320774217d0 /
    data xgk (  7) / 0.790418501442465932967649294817947d0 /
    data xgk (  8) / 0.724417731360170047416186054613938d0 /
    data xgk (  9) / 0.650996741297416970533735895313275d0 /
    data xgk ( 10) / 0.570972172608538847537226737253911d0 /
    data xgk ( 11) / 0.485081863640239680693655740232351d0 /
    data xgk ( 12) / 0.394151347077563369897207370981045d0 /
    data xgk ( 13) / 0.299180007153168812166780024266389d0 /
    data xgk ( 14) / 0.201194093997434522300628303394596d0 /
    data xgk ( 15) / 0.101142066918717499027074231447392d0 /
    data xgk ( 16) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.005377479872923348987792051430128d0 /
    data wgk (  2) / 0.015007947329316122538374763075807d0 /
    data wgk (  3) / 0.025460847326715320186874001019653d0 /
    data wgk (  4) / 0.035346360791375846222037948478360d0 /
    data wgk (  5) / 0.044589751324764876608227299373280d0 /
    data wgk (  6) / 0.053481524690928087265343147239430d0 /
    data wgk (  7) / 0.062009567800670640285139230960803d0 /
    data wgk (  8) / 0.069854121318728258709520077099147d0 /
    data wgk (  9) / 0.076849680757720378894432777482659d0 /
    data wgk ( 10) / 0.083080502823133021038289247286104d0 /
    data wgk ( 11) / 0.088564443056211770647275443693774d0 /
    data wgk ( 12) / 0.093126598170825321225486872747346d0 /
    data wgk ( 13) / 0.096642726983623678505179907627589d0 /
    data wgk ( 14) / 0.099173598721791959332393173484603d0 /
    data wgk ( 15) / 0.100769845523875595044946662617570d0 /
    data wgk ( 16) / 0.101330007014791549017374792767493d0 /


!           list of major variables
!           -----------------------
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 15-point gauss formula
!           resk   - result of the 31-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)

!           machine dependent constants
!           ---------------------------
!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.
!***first executable statement  dqk31
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 31-point kronrod approximation to
!           the integral, and estimate the absolute error.

    fc = f(centr)
    resg = wg(8)*fc
    resk = wgk(16)*fc
    resabs = dabs(resk)
    do j=1,7
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,8
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5d+00
    resasc = wgk(16)*dabs(fc-reskh)
    do j=1,15
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.0d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
end subroutine dqk31


subroutine dqk41(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk41
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  41-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description

!           integration rules
!           standard fortran subroutine
!           double precision version

!           parameters
!            on entry
!              f      - double precision
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.

!              a      - double precision
!                       lower limit of integration

!              b      - double precision
!                       upper limit of integration

!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 41-point
!                       gauss-kronrod rule (resk) obtained by optimal
!                       addition of abscissae to the 20-point gauss
!                       rule (resg).

!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integal of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk41

    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    !external f

    dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 41-point gauss-kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 20-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 20-point gauss rule

!           wgk    - weights of the 41-point gauss-kronrod rule

!           wg     - weights of the 20-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.017614007139152118311861962351853d0 /
    data wg  (  2) / 0.040601429800386941331039952274932d0 /
    data wg  (  3) / 0.062672048334109063569506535187042d0 /
    data wg  (  4) / 0.083276741576704748724758143222046d0 /
    data wg  (  5) / 0.101930119817240435036750135480350d0 /
    data wg  (  6) / 0.118194531961518417312377377711382d0 /
    data wg  (  7) / 0.131688638449176626898494499748163d0 /
    data wg  (  8) / 0.142096109318382051329298325067165d0 /
    data wg  (  9) / 0.149172986472603746787828737001969d0 /
    data wg  ( 10) / 0.152753387130725850698084331955098d0 /

    data xgk (  1) / 0.998859031588277663838315576545863d0 /
    data xgk (  2) / 0.993128599185094924786122388471320d0 /
    data xgk (  3) / 0.981507877450250259193342994720217d0 /
    data xgk (  4) / 0.963971927277913791267666131197277d0 /
    data xgk (  5) / 0.940822633831754753519982722212443d0 /
    data xgk (  6) / 0.912234428251325905867752441203298d0 /
    data xgk (  7) / 0.878276811252281976077442995113078d0 /
    data xgk (  8) / 0.839116971822218823394529061701521d0 /
    data xgk (  9) / 0.795041428837551198350638833272788d0 /
    data xgk ( 10) / 0.746331906460150792614305070355642d0 /
    data xgk ( 11) / 0.693237656334751384805490711845932d0 /
    data xgk ( 12) / 0.636053680726515025452836696226286d0 /
    data xgk ( 13) / 0.575140446819710315342946036586425d0 /
    data xgk ( 14) / 0.510867001950827098004364050955251d0 /
    data xgk ( 15) / 0.443593175238725103199992213492640d0 /
    data xgk ( 16) / 0.373706088715419560672548177024927d0 /
    data xgk ( 17) / 0.301627868114913004320555356858592d0 /
    data xgk ( 18) / 0.227785851141645078080496195368575d0 /
    data xgk ( 19) / 0.152605465240922675505220241022678d0 /
    data xgk ( 20) / 0.076526521133497333754640409398838d0 /
    data xgk ( 21) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.003073583718520531501218293246031d0 /
    data wgk (  2) / 0.008600269855642942198661787950102d0 /
    data wgk (  3) / 0.014626169256971252983787960308868d0 /
    data wgk (  4) / 0.020388373461266523598010231432755d0 /
    data wgk (  5) / 0.025882133604951158834505067096153d0 /
    data wgk (  6) / 0.031287306777032798958543119323801d0 /
    data wgk (  7) / 0.036600169758200798030557240707211d0 /
    data wgk (  8) / 0.041668873327973686263788305936895d0 /
    data wgk (  9) / 0.046434821867497674720231880926108d0 /
    data wgk ( 10) / 0.050944573923728691932707670050345d0 /
    data wgk ( 11) / 0.055195105348285994744832372419777d0 /
    data wgk ( 12) / 0.059111400880639572374967220648594d0 /
    data wgk ( 13) / 0.062653237554781168025870122174255d0 /
    data wgk ( 14) / 0.065834597133618422111563556969398d0 /
    data wgk ( 15) / 0.068648672928521619345623411885368d0 /
    data wgk ( 16) / 0.071054423553444068305790361723210d0 /
    data wgk ( 17) / 0.073030690332786667495189417658913d0 /
    data wgk ( 18) / 0.074582875400499188986581418362488d0 /
    data wgk ( 19) / 0.075704497684556674659542775376617d0 /
    data wgk ( 20) / 0.076377867672080736705502835038061d0 /
    data wgk ( 21) / 0.076600711917999656445049901530102d0 /


!           list of major variables
!           -----------------------

!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 20-point gauss formula
!           resk   - result of the 41-point kronrod formula
!           reskh  - approximation to mean value of f over (a,b), i.e.
!                    to i/(b-a)

!           machine dependent constants
!           ---------------------------

!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

!***first executable statement  dqk41
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 41-point gauss-kronrod approximation to
!           the integral, and estimate the absolute error.

    resg = 0.0d+00
    fc = f(centr)
    resk = wgk(21)*fc
    resabs = dabs(resk)
    do j=1,10
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,10
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5d+00
    resasc = wgk(21)*dabs(fc-reskh)
    do j=1,20
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
end subroutine dqk41
    

subroutine dqk51(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk51
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  51-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description

!           integration rules
!           standard fortran subroutine
!           double precision version

!           parameters
!            on entry
!              f      - double precision
!                       function subroutine defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.

!              a      - double precision
!                       lower limit of integration

!              b      - double precision
!                       upper limit of integration

!            on return
!              result - double precision
!                       approximation to the integral i
!                       result is computed by applying the 51-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 25-point gauss rule (resg).

!              abserr - double precision
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)

!              resabs - double precision
!                       approximation to the integral j

!              resasc - double precision
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)

!***references  (none)
!***routines called  d1mach
!***end prologue  dqk51

    double precision :: a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    !external f

    dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)

!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.

!           xgk    - abscissae of the 51-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 25-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 25-point gauss rule

!           wgk    - weights of the 51-point kronrod rule

!           wg     - weights of the 25-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.011393798501026287947902964113235d0 /
    data wg  (  2) / 0.026354986615032137261901815295299d0 /
    data wg  (  3) / 0.040939156701306312655623487711646d0 /
    data wg  (  4) / 0.054904695975835191925936891540473d0 /
    data wg  (  5) / 0.068038333812356917207187185656708d0 /
    data wg  (  6) / 0.080140700335001018013234959669111d0 /
    data wg  (  7) / 0.091028261982963649811497220702892d0 /
    data wg  (  8) / 0.100535949067050644202206890392686d0 /
    data wg  (  9) / 0.108519624474263653116093957050117d0 /
    data wg  ( 10) / 0.114858259145711648339325545869556d0 /
    data wg  ( 11) / 0.119455763535784772228178126512901d0 /
    data wg  ( 12) / 0.122242442990310041688959518945852d0 /
    data wg  ( 13) / 0.123176053726715451203902873079050d0 /

    data xgk (  1) / 0.999262104992609834193457486540341d0 /
    data xgk (  2) / 0.995556969790498097908784946893902d0 /
    data xgk (  3) / 0.988035794534077247637331014577406d0 /
    data xgk (  4) / 0.976663921459517511498315386479594d0 /
    data xgk (  5) / 0.961614986425842512418130033660167d0 /
    data xgk (  6) / 0.942974571228974339414011169658471d0 /
    data xgk (  7) / 0.920747115281701561746346084546331d0 /
    data xgk (  8) / 0.894991997878275368851042006782805d0 /
    data xgk (  9) / 0.865847065293275595448996969588340d0 /
    data xgk ( 10) / 0.833442628760834001421021108693570d0 /
    data xgk ( 11) / 0.797873797998500059410410904994307d0 /
    data xgk ( 12) / 0.759259263037357630577282865204361d0 /
    data xgk ( 13) / 0.717766406813084388186654079773298d0 /
    data xgk ( 14) / 0.673566368473468364485120633247622d0 /
    data xgk ( 15) / 0.626810099010317412788122681624518d0 /
    data xgk ( 16) / 0.577662930241222967723689841612654d0 /
    data xgk ( 17) / 0.526325284334719182599623778158010d0 /
    data xgk ( 18) / 0.473002731445714960522182115009192d0 /
    data xgk ( 19) / 0.417885382193037748851814394594572d0 /
    data xgk ( 20) / 0.361172305809387837735821730127641d0 /
    data xgk ( 21) / 0.303089538931107830167478909980339d0 /
    data xgk ( 22) / 0.243866883720988432045190362797452d0 /
    data xgk ( 23) / 0.183718939421048892015969888759528d0 /
    data xgk ( 24) / 0.122864692610710396387359818808037d0 /
    data xgk ( 25) / 0.061544483005685078886546392366797d0 /
    data xgk ( 26) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.001987383892330315926507851882843d0 /
    data wgk (  2) / 0.005561932135356713758040236901066d0 /
    data wgk (  3) / 0.009473973386174151607207710523655d0 /
    data wgk (  4) / 0.013236229195571674813656405846976d0 /
    data wgk (  5) / 0.016847817709128298231516667536336d0 /
    data wgk (  6) / 0.020435371145882835456568292235939d0 /
    data wgk (  7) / 0.024009945606953216220092489164881d0 /
    data wgk (  8) / 0.027475317587851737802948455517811d0 /
    data wgk (  9) / 0.030792300167387488891109020215229d0 /
    data wgk ( 10) / 0.034002130274329337836748795229551d0 /
    data wgk ( 11) / 0.037116271483415543560330625367620d0 /
    data wgk ( 12) / 0.040083825504032382074839284467076d0 /
    data wgk ( 13) / 0.042872845020170049476895792439495d0 /
    data wgk ( 14) / 0.045502913049921788909870584752660d0 /
    data wgk ( 15) / 0.047982537138836713906392255756915d0 /
    data wgk ( 16) / 0.050277679080715671963325259433440d0 /
    data wgk ( 17) / 0.052362885806407475864366712137873d0 /
    data wgk ( 18) / 0.054251129888545490144543370459876d0 /
    data wgk ( 19) / 0.055950811220412317308240686382747d0 /
    data wgk ( 20) / 0.057437116361567832853582693939506d0 /
    data wgk ( 21) / 0.058689680022394207961974175856788d0 /
    data wgk ( 22) / 0.059720340324174059979099291932562d0 /
    data wgk ( 23) / 0.060539455376045862945360267517565d0 /
    data wgk ( 24) / 0.061128509717053048305859030416293d0 /
    data wgk ( 25) / 0.061471189871425316661544131965264d0 /
    data wgk ( 26) / 0.061580818067832935078759824240066d0 /
!       note: wgk (26) was calculated from the values of wgk(1..25)


!           list of major variables
!           -----------------------

!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 25-point gauss formula
!           resk   - result of the 51-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)

!           machine dependent constants
!           ---------------------------

!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

!***first executable statement  dqk51
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 51-point kronrod approximation to
!           the integral, and estimate the absolute error.

    fc = f(centr)
    resg = wg(13)*fc
    resk = wgk(26)*fc
    resabs = dabs(resk)
    do j=1,12
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j = 1,13
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc)
        fval2 = f(centr+absc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5d+00
    resasc = wgk(26)*dabs(fc-reskh)
    do j=1,25
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.0d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
end subroutine dqk51

subroutine dqk61(f,a,b,result,abserr,resabs,resasc)
!***begin prologue  dqk61
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  61-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                           estimate
!                       j = integral of dabs(f) over (a,b)
!***description

!        integration rule
!        standard fortran subroutine
!        double precision version


!        parameters
!         on entry
!           f      - double precision
!                    function subprogram defining the integrand
!                    function f(x). the actual name for f needs to be
!                    declared e x t e r n a l in the calling program.

!           a      - double precision
!                    lower limit of integration

!           b      - double precision
!                    upper limit of integration

!         on return
!           result - double precision
!                    approximation to the integral i
!                    result is computed by applying the 61-point
!                    kronrod rule (resk) obtained by optimal addition of
!                    abscissae to the 30-point gauss rule (resg).

!           abserr - double precision
!                    estimate of the modulus of the absolute error,
!                    which should equal or exceed dabs(i-result)

!           resabs - double precision
!                    approximation to the integral j

!           resasc - double precision
!                    approximation to the integral of dabs(f-i/(b-a))


!***references  (none)
!***routines called  d1mach
!***end prologue  dqk61

    double precision :: a,dabsc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1, &
    d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
    integer :: j,jtw,jtwm1
    !external f

    dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)

!           the abscissae and weights are given for the
!           interval (-1,1). because of symmetry only the positive
!           abscissae and their corresponding weights are given.

!           xgk   - abscissae of the 61-point kronrod rule
!                   xgk(2), xgk(4)  ... abscissae of the 30-point
!                   gauss rule
!                   xgk(1), xgk(3)  ... optimally added abscissae
!                   to the 30-point gauss rule

!           wgk   - weights of the 61-point kronrod rule

!           wg    - weigths of the 30-point gauss rule


! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

    data wg  (  1) / 0.007968192496166605615465883474674d0 /
    data wg  (  2) / 0.018466468311090959142302131912047d0 /
    data wg  (  3) / 0.028784707883323369349719179611292d0 /
    data wg  (  4) / 0.038799192569627049596801936446348d0 /
    data wg  (  5) / 0.048402672830594052902938140422808d0 /
    data wg  (  6) / 0.057493156217619066481721689402056d0 /
    data wg  (  7) / 0.065974229882180495128128515115962d0 /
    data wg  (  8) / 0.073755974737705206268243850022191d0 /
    data wg  (  9) / 0.080755895229420215354694938460530d0 /
    data wg  ( 10) / 0.086899787201082979802387530715126d0 /
    data wg  ( 11) / 0.092122522237786128717632707087619d0 /
    data wg  ( 12) / 0.096368737174644259639468626351810d0 /
    data wg  ( 13) / 0.099593420586795267062780282103569d0 /
    data wg  ( 14) / 0.101762389748405504596428952168554d0 /
    data wg  ( 15) / 0.102852652893558840341285636705415d0 /

    data xgk (  1) / 0.999484410050490637571325895705811d0 /
    data xgk (  2) / 0.996893484074649540271630050918695d0 /
    data xgk (  3) / 0.991630996870404594858628366109486d0 /
    data xgk (  4) / 0.983668123279747209970032581605663d0 /
    data xgk (  5) / 0.973116322501126268374693868423707d0 /
    data xgk (  6) / 0.960021864968307512216871025581798d0 /
    data xgk (  7) / 0.944374444748559979415831324037439d0 /
    data xgk (  8) / 0.926200047429274325879324277080474d0 /
    data xgk (  9) / 0.905573307699907798546522558925958d0 /
    data xgk ( 10) / 0.882560535792052681543116462530226d0 /
    data xgk ( 11) / 0.857205233546061098958658510658944d0 /
    data xgk ( 12) / 0.829565762382768397442898119732502d0 /
    data xgk ( 13) / 0.799727835821839083013668942322683d0 /
    data xgk ( 14) / 0.767777432104826194917977340974503d0 /
    data xgk ( 15) / 0.733790062453226804726171131369528d0 /
    data xgk ( 16) / 0.697850494793315796932292388026640d0 /
    data xgk ( 17) / 0.660061064126626961370053668149271d0 /
    data xgk ( 18) / 0.620526182989242861140477556431189d0 /
    data xgk ( 19) / 0.579345235826361691756024932172540d0 /
    data xgk ( 20) / 0.536624148142019899264169793311073d0 /
    data xgk ( 21) / 0.492480467861778574993693061207709d0 /
    data xgk ( 22) / 0.447033769538089176780609900322854d0 /
    data xgk ( 23) / 0.400401254830394392535476211542661d0 /
    data xgk ( 24) / 0.352704725530878113471037207089374d0 /
    data xgk ( 25) / 0.304073202273625077372677107199257d0 /
    data xgk ( 26) / 0.254636926167889846439805129817805d0 /
    data xgk ( 27) / 0.204525116682309891438957671002025d0 /
    data xgk ( 28) / 0.153869913608583546963794672743256d0 /
    data xgk ( 29) / 0.102806937966737030147096751318001d0 /
    data xgk ( 30) / 0.051471842555317695833025213166723d0 /
    data xgk ( 31) / 0.000000000000000000000000000000000d0 /

    data wgk (  1) / 0.001389013698677007624551591226760d0 /
    data wgk (  2) / 0.003890461127099884051267201844516d0 /
    data wgk (  3) / 0.006630703915931292173319826369750d0 /
    data wgk (  4) / 0.009273279659517763428441146892024d0 /
    data wgk (  5) / 0.011823015253496341742232898853251d0 /
    data wgk (  6) / 0.014369729507045804812451432443580d0 /
    data wgk (  7) / 0.016920889189053272627572289420322d0 /
    data wgk (  8) / 0.019414141193942381173408951050128d0 /
    data wgk (  9) / 0.021828035821609192297167485738339d0 /
    data wgk ( 10) / 0.024191162078080601365686370725232d0 /
    data wgk ( 11) / 0.026509954882333101610601709335075d0 /
    data wgk ( 12) / 0.028754048765041292843978785354334d0 /
    data wgk ( 13) / 0.030907257562387762472884252943092d0 /
    data wgk ( 14) / 0.032981447057483726031814191016854d0 /
    data wgk ( 15) / 0.034979338028060024137499670731468d0 /
    data wgk ( 16) / 0.036882364651821229223911065617136d0 /
    data wgk ( 17) / 0.038678945624727592950348651532281d0 /
    data wgk ( 18) / 0.040374538951535959111995279752468d0 /
    data wgk ( 19) / 0.041969810215164246147147541285970d0 /
    data wgk ( 20) / 0.043452539701356069316831728117073d0 /
    data wgk ( 21) / 0.044814800133162663192355551616723d0 /
    data wgk ( 22) / 0.046059238271006988116271735559374d0 /
    data wgk ( 23) / 0.047185546569299153945261478181099d0 /
    data wgk ( 24) / 0.048185861757087129140779492298305d0 /
    data wgk ( 25) / 0.049055434555029778887528165367238d0 /
    data wgk ( 26) / 0.049795683427074206357811569379942d0 /
    data wgk ( 27) / 0.050405921402782346840893085653585d0 /
    data wgk ( 28) / 0.050881795898749606492297473049805d0 /
    data wgk ( 29) / 0.051221547849258772170656282604944d0 /
    data wgk ( 30) / 0.051426128537459025933862879215781d0 /
    data wgk ( 31) / 0.051494729429451567558340433647099d0 /

!           list of major variables
!           -----------------------

!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           dabsc  - abscissa
!           fval*  - function value
!           resg   - result of the 30-point gauss rule
!           resk   - result of the 61-point kronrod rule
!           reskh  - approximation to the mean value of f
!                    over (a,b), i.e. to i/(b-a)

!           machine dependent constants
!           ---------------------------

!           epmach is the largest relative spacing.
!           uflow is the smallest positive magnitude.

    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(b+a)
    hlgth = 0.5d+00*(b-a)
    dhlgth = dabs(hlgth)

!           compute the 61-point kronrod approximation to the
!           integral, and estimate the absolute error.

!***first executable statement  dqk61
    resg = 0.0d+00
    fc = f(centr)
    resk = wgk(31)*fc
    resabs = dabs(resk)
    do j=1,15
        jtw = j*2
        dabsc = hlgth*xgk(jtw)
        fval1 = f(centr-dabsc)
        fval2 = f(centr+dabsc)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(dabs(fval1)+dabs(fval2))
    end do
    do j=1,15
        jtwm1 = j*2-1
        dabsc = hlgth*xgk(jtwm1)
        fval1 = f(centr-dabsc)
        fval2 = f(centr+dabsc)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(dabs(fval1)+dabs(fval2))
    end do
    reskh = resk*0.5d+00
    resasc = wgk(31)*dabs(fc-reskh)
    do j=1,30
        resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
    end do
    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = dabs((resk-resg)*hlgth)
    if(resasc /= 0.0d+00 .AND. abserr /= 0.0d+00) &
    abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
    if(resabs > uflow/(0.5d+02*epmach)) abserr = dmax1 &
    ((epmach*0.5d+02)*resabs,abserr)
    return
end subroutine dqk61

subroutine dqpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
!***begin prologue  dqpsrt
!***refer to  dqage,dqagie,dqagpe,dqawse
!***routines called  (none)
!***revision date  810101   (yymmdd)
!***keywords  sequential sorting
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  this routine maintains the descending ordering in the
!            list of the local error estimated resulting from the
!            interval subdivision process. at each call two error
!            estimates are inserted using the sequential search
!            method, top-down for the largest error estimate and
!            bottom-up for the smallest error estimate.
!***description

!           ordering routine
!           standard fortran subroutine
!           double precision version

!           parameters (meaning at output)
!              limit  - integer
!                       maximum number of error estimates the list
!                       can contain

!              last   - integer
!                       number of error estimates currently in the list

!              maxerr - integer
!                       maxerr points to the nrmax-th largest error
!                       estimate currently in the list

!              ermax  - double precision
!                       nrmax-th largest error estimate
!                       ermax = elist(maxerr)

!              elist  - double precision
!                       vector of dimension last containing
!                       the error estimates

!              iord   - integer
!                       vector of dimension last, the first k elements
!                       of which contain pointers to the error
!                       estimates, such that
!                       elist(iord(1)),...,  elist(iord(k))
!                       form a decreasing sequence, with
!                       k = last if last.le.(limit/2+2), and
!                       k = limit+1-last otherwise

!              nrmax  - integer
!                       maxerr = iord(nrmax)

!***end prologue  dqpsrt

    double precision :: elist,ermax,errmax,errmin
    integer :: i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr, &
    nrmax
    dimension elist(last),iord(last)

!           check whether the list contains more than
!           two error estimates.

!***first executable statement  dqpsrt
    if(last > 2) go to 10
    iord(1) = 1
    iord(2) = 2
    go to 90

!           this part of the routine is only executed if, due to a
!           difficult integrand, subdivision increased the error
!           estimate. in the normal case the insert procedure should
!           start after the nrmax-th largest error estimate.

    10 errmax = elist(maxerr)
    if(nrmax == 1) go to 30
    ido = nrmax-1
    do i = 1,ido
        isucc = iord(nrmax-1)
    ! ***jump out of do-loop
        if(errmax <= elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
    end do

!           compute the number of elements in the list to be maintained
!           in descending order. this number depends on the number of
!           subdivisions still allowed.

    30 jupbn = last
    if(last > (limit/2+2)) jupbn = limit+3-last
    errmin = elist(last)

!           insert errmax by traversing the list top-down,
!           starting comparison from the element elist(iord(nrmax+1)).

    jbnd = jupbn-1
    ibeg = nrmax+1
    if(ibeg > jbnd) go to 50
    do i=ibeg,jbnd
        isucc = iord(i)
    ! ***jump out of do-loop
        if(errmax >= elist(isucc)) go to 60
        iord(i-1) = isucc
    end do
    50 iord(jbnd) = maxerr
    iord(jupbn) = last
    go to 90

!           insert errmin by traversing the list bottom-up.

    60 iord(i-1) = maxerr
    k = jbnd
    do j=i,jbnd
        isucc = iord(k)
    ! ***jump out of do-loop
        if(errmin < elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
    end do
    iord(i) = last
    go to 90
    80 iord(k+1) = last

!           set maxerr and ermax.

    90 maxerr = iord(nrmax)
    ermax = elist(maxerr)
    return
end subroutine dqpsrt
    