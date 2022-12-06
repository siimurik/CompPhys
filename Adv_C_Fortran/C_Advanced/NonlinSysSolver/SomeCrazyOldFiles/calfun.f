SUBROUTINE CALFUN (N,X,F)
DIMENSION X(1), F(1)
F(l) = -13.+X(l) + ((-X(2) + 5.)*X(2)-2. )*X(2)
F(2) = -29.+X(1) + (( X(2) + l.)*X(2)-14.)*X(2)
RETURN
END