//gcc nrgaujac.c -o gj -lnrutil -lm
#include <stdio.h>
#include <math.h>
#define EPS 3.0e-14 // Increase EPS if you don’t have this preci-
#define MAXIT 10    // sion.

double gammln(double xx)
// Returns the value ln[Γ(xx)] for xx > 0.
{
    // Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
    // accuracy is good enough.
    double x, y, tmp, ser;
    static double cof[6]={  76.18009172947146,      -86.50532032941677,
                            24.01409824083091,      -1.231739572450155,
                            0.1208650973866179e-2,  -0.5395239384953e-5};
    int j;

    y = x = xx;
    tmp  = x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser  =1.000000000190015;
    for (j=0; j<=5; j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

void gaujac(double x[], double w[], int n, double alf, double bet)
// Given alf and bet, the parameters α and β of the Jacobi polynomials, this routine returns
// arrays x[1..n] and w[1..n] containing the abscissas and weights of the n-point Gauss-Jacobi
// quadrature formula. The largest abscissa is returned in x[1], the smallest in x[n].
{
    double gammln(double xx);
    void nrerror(char error_text[]);
    int i, its, j;
    double alfbet, an, bn, r1, r2, r3;
    double a, b, c, p1, p2, p3, pp, temp, z, z1; // High precision is a good idea for this rou-
                                // tine.
    for (i=1; i<=n; i++) {      // Loop over the desired roots.
        if (i == 1) {           // Initial guess for the largest root.
            an = alf/n;
            bn = bet/n;
            r1 = (1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
            r2 = 1.0 + 1.48*an + 0.96*bn + 0.452*an*an + 0.83*an*bn;
            z  = 1.0 - r1/r2;
        } else if (i == 2) {    // Initial guess for the second largest root.
            r1 = (4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
            r2 = 1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
            r3 = 1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
            z -= (1.0-z)*r1*r2*r3;
        } else if (i == 3) {    // Initial guess for the third largest root.
            r1 = (1.67+0.28*alf)/(1.0+0.37*alf);
            r2 = 1.0+0.22*(n-8.0)/n;
            r3 = 1.0+8.0*bet/((6.28+bet)*n*n);
            z -= (x[1]-z)*r1*r2*r3;
        } else if (i == n-1) {  // Initial guess for the second smallest root.
            r1 = (1.0+0.235*bet)/(0.766+0.119*bet);
            r2 = 1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
            r3 = 1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
            z += (z-x[n-3])*r1*r2*r3;
        } else if (i == n) {    // Initial guess for the smallest root.
            r1 = (1.0+0.37*bet)/(1.67+0.28*bet);
            r2 = 1.0/(1.0+0.22*(n-8.0)/n);
            r3 = 1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
            z += (z-x[n-2])*r1*r2*r3;
        } else {                // Initial guess for the other roots.
            z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
        }
        alfbet=alf+bet;
        for (its=1; its<=MAXIT; its++) {    // Refinement by Newton’s method.
            temp = 2.0+alfbet;              // Start the recurrence with P0 and P1 to avoid
                                            // a division by zero when α + β = 0 or −1.
            p1 = (alf-bet+temp*z)/2.0;
            p2 = 1.0;
            for (j=2; j<=n; j++) {  // Loop up the recurrence relation to get the
                p3 = p2;            // Jacobi polynomial evaluated at z.
                p2 = p1;
                temp = 2*j+alfbet;
                a  = 2*j*(j+alfbet)*(temp-2.0);
                b  = (temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
                c  = 2.0*(j-1+alf)*(j-1+bet)*temp;
                p1 = (b*p2-c*p3)/a;
            }
            pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
            // p1 is now the desired Jacobi polynomial. We next compute pp, its derivative, by
            // a standard relation involving also p2, the polynomial of one lower order.
            z1 = z;
            z  = z1-p1/pp; // Newton’s formula.
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) nrerror("too many iterations in gaujac");
        x[i] = z;   // Store the root and the weight.
        w[i] = exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
                gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);
        // For debugging
        //printf("x[%d] = %f \t w[%d] = %f\n", i, x[i], i, w[i]);
    }
}

double f(double x)
{
    return exp(-x*x);
}

double compute_integral(double (*func)(double), double a, double b, int n)
{
    int i;
    double alf =  0.0, bet = 0.0; // Initialize alf and bet
    //double a   = -2.0, b   = 2.0; // Define the lower and upper bounds
    double x[n], w[n];
    double m = (b - a) / 2.0;
    double c = (b + a) / 2.0;
    double integral = 0.0;

    gaujac(x, w, n, alf, bet);
    for (i = 0; i < n; i++)
    {
        integral += w[i] * func(m * x[i] + c);
    }
    integral *= m;

    return integral;
}

int main()
{   
    int i, n = 128;
    double a, b, integral; // Define the lower and upper bounds
    int nint;

    // Integration limits and accuraccy.
    a   = -2.0;
    b   =  2.0;

    integral = compute_integral(f, a, b, n);
    printf("Approximate integral: %f\n", integral);

    return 0;
}
