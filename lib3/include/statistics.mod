GFORTRAN module version '10' created from /home/jaemoo/seadas-7.4/ocssw/build/lib3/src/healpix/Healpix_3.11/src/f90/mod/statistics.f90
MD5:9228ef672875395b44e17d7de061caaa -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

(('compute_statistics' 'statistics' 2 3) ('median' 'statistics' 4 5) (
'tstats' 'statistics' 6))

()

()

()

(6 'Tstats' 'statistics' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 0 ((7 'ntot'
(INTEGER 8 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (8 'nvalid' (INTEGER 8 0 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (9 'mind' (REAL 8 0 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (10
'maxd' (REAL 8 0 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (11 'average' (REAL 8 0 0 0 REAL
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (12 'absdev' (REAL 8 0 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (13
'rms' (REAL 8 0 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (14 'var' (REAL 8 0 0 0 REAL ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (15 'skew' (REAL 8 0 0 0 REAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (16
'kurt' (REAL 8 0 0 0 REAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
54237575)
17 'print_statistics' 'statistics' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ()) 18
0 (19) () 0 () () () 0 0)
20 'tstats' 'statistics' '' 1 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC
DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () ()
0 () () () 0 0)
2 'comp_stats_d' 'statistics' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 21 0 (22 23 24) () 0 () () () 0 0)
3 'comp_stats_s' 'statistics' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 25 0 (26 27 28) () 0 () () () 0 0)
4 'median_d' 'statistics' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ()) 29 0 (
30 31 32) () 33 () () () 0 0)
5 'median_s' 'statistics' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (REAL 4 0 0 0 REAL ()) 34 0 (
35 36 37) () 38 () () () 0 0)
19 'stats' '' '' 18 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 6 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
22 'data' '' '' 21 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0') ()) 0 () () () 0 0)
23 'stats' '' '' 21 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 6 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
24 'badval' '' '' 21 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
26 'data' '' '' 25 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0') ()) 0 () () () 0 0)
27 'stats' '' '' 25 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 6 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
28 'badval' '' '' 25 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
30 'data' '' '' 29 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION TARGET DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
31 'badval' '' '' 29 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
32 'even' '' '' 29 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
33 'med' '' '' 29 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 RESULT ALWAYS_EXPLICIT) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () () ()
0 0)
35 'data' '' '' 34 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION TARGET DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
36 'badval' '' '' 34 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () () 0 0)
37 'even' '' '' 34 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
38 'med' '' '' 34 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 RESULT ALWAYS_EXPLICIT) (REAL 4 0 0 0 REAL ()) 0 0 () () 0 () () ()
0 0)
)

('Tstats' 0 6 'print_statistics' 0 17 'tstats' 0 20)
