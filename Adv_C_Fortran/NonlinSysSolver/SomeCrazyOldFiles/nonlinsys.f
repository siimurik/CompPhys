DIMENSION X(2),F(2),AJINV(2,2),W(100)
X(1)    = 15.
X(2)    = -2.
N       = 2
STEP    = 0.01
ACC     = 0.000001
MAXFUN  = 100
DMAX    = 10.
CALL NSOIA (2,X,F,AJINV,STEPtDMAX,ACCtMAXFUN,l,W)
CALL EXIT
END