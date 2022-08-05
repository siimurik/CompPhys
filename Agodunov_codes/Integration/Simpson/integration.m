clear
format short g
disp("Numerical integration based on Gaussian quadrature.");
tic
[q, ier, nfun, err] = quad("f", 0, pi)
time = toc
disp("\nNumerical integration using an adaptive vectorized Simpson’s rule.");
tic
[q1, nfun1] = quadv("f", 0, pi)
time = toc 
disp('\n');
diff = 3.15593897269467898E-04 - 0.31559390E-03
