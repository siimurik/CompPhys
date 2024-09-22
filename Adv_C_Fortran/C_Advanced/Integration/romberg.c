// Compile and execute with
// 	$ gcc polintd.c romberg.c -lnrutil -lm -o rom
// 	$ ./rom
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#define EPS 1.0e-14
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5
#define FUNC(x) ((*func)(x))
#define CLOCK_MONOTONIC 1

double func(double x);
double trapzd(double (*func)(double), double a, double b, int n);
double qromb(double (*func)(double), double a, double b);

int main()
{
	double a, b, s;
    struct timespec start, stop;

	a = -2.0;
	b =  2.0;

	// Get starting time
    clock_gettime(CLOCK_MONOTONIC, &start);

	s = qromb(func, a, b);

    // Get end time
    clock_gettime(CLOCK_MONOTONIC, &stop);

	printf("Result: %.16f\n", s);

	// Calculate the elapsed time in seconds
    double time_taken = (stop.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (stop.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Elapsed time: %.4e.\n", time_taken);

	return 0;
}

double func(double x)
{
	return exp(-x*x);
}

double trapzd(double (*func)(double), double a, double b, int n)
// This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
// as a pointer to the function to be integrated between limits a and b, also input. When called with
// n=1, the routine returns the crudest estimate of ∫ b
// a f (x)dx. Subsequent calls with n=2,3,...
// (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
{
	double x, tnm, sum, del;
	static double s;
	int it, j;
	if (n == 1) {
		return (s = 0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1, j=1; j<n-1; j++) it <<= 1;
		tnm = it;
		del = (b-a)/tnm; // This is the spacing of the points to be added.
		x = a+0.5*del;
		for (sum = 0.0, j=1; j<=it; j++, x+=del) sum += FUNC(x);
		s = 0.5*(s+(b-a)*sum/tnm); //This replaces s by its refined value.
		return s;
	}
}
	
// Here EPS is the fractional accuracy desired, as determined by the extrapolation error estimate;
// JMAX limits the total number of steps; K is the number of points used in the extrapolation.
double qromb(double (*func)(double), double a, double b)
// Returns the integral of the function func from a to b. Integration is performed by Romberg’s
// method of order 2K, where, e.g., K=2 is Simpson’s rule.
{
	void polintd(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	double ss, dss;
	double s[JMAXP],h[JMAXP+1]; // These store the successive trapezoidal approxi-
	int j; //mations and their relative stepsizes.
	
	h[1] = 1.0;
	for (j=1; j<=JMAX; j++) {
		s[j] = trapzd(func, a, b, j);
		if (j >= K) {
			polintd(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1] = 0.25*h[j];
		// This is a key step: The factor is 0.25 even though the stepsize is decreased by only
		// 0.5. This makes the extrapolation a polynomial in h2 as allowed by equation (4.2.1),
		// not just a polynomial in h.
	}
	nrerror("Too many steps in routine qromb");
	return 0.0; // Never get here.
}
