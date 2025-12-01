#include "nrutil.h"

void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
void (*derivs)(float, float [], float []))
// Given values for the variables y[1..n] and their derivatives dydx[1..n] known at x, use the
// fourth-order Runge-Kutta method to advance the solution over an interval h and return the
// incremented variables as yout[1..n], which need not be a distinct array from y. The user
// supplies the routine derivs(x,y,dydx), which returns derivatives dydx at x.
{
    int i;
    float xh, hh, h6, *dym, *dyt,* yt;

    dym=vector(1,n);
    dyt=vector(1,n);
    yt=vector(1,n);

    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;

    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];   // First step.
    (*derivs)(xh,yt,dyt);                       // Second step.
    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
    (*derivs)(xh,yt,dym);                       // Third step.

    for (i=1;i<=n;i++) {
        yt[i]=y[i]+h*dym[i];
        dym[i] += dyt[i];
    }

    (*derivs)(x+h,yt,dyt);                      // Fourth step.
    for (i=1;i<=n;i++)                              // Accumulate increments with proper
        yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);    // weights.

    free_vector(yt,1,n);
    free_vector(dyt,1,n);
    free_vector(dym,1,n);
}

#include "nrutil.h"
float **y,*xx; // For communication back to main.
void rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
void (*derivs)(float, float [], float []))
// Starting from initial values vstart[1..nvar] known at x1 use fourth-order Runge-Kutta
// to advance nstep equal increments to x2. The user-supplied routine derivs(x,v,dvdx)
// evaluates derivatives. Results are stored in the global variables y[1..nvar][1..nstep+1]
// and xx[1..nstep+1].
{
    void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
    void (*derivs)(float, float [], float []));
    int i,k;
    float x,h;
    float *v,*vout,*dv;

    v=vector(1,nvar);
    vout=vector(1,nvar);
    dv=vector(1,nvar);

    for (i=1;i<=nvar;i++) { // Load starting values.
        v[i]=vstart[i];
        y[i][1]=v[i];
    }

    xx[1]=x1;
    x=x1;
    h=(x2-x1)/nstep;

    for (k=1;k<=nstep;k++) { // Take nstep steps.
        (*derivs)(x,v,dv);
        rk4(v,dv,nvar,x,h,vout,derivs);
        if ((float)(x+h) == x) nrerror("Step size too small in routine rkdumb");
        x += h;
        xx[k+1]=x; // Store intermediate steps.
        for (i=1;i<=nvar;i++) {
            v[i]=vout[i];
            y[i][k+1]=v[i];
        }
    }

    free_vector(dv,1,nvar);
    free_vector(vout,1,nvar);
    free_vector(v,1,nvar);
}