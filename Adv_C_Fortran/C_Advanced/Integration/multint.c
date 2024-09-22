// Compile and execute with:
//  $ gcc multint.c -o mi -lm
//  $ ./mi

#include <stdio.h>
#include <math.h>

// Declare the common block variables as global variables
static float xsav, ysav;
static float (*nrfunc)(float, float, float);

// Function prototypes
float quad3d(float (*func)(float, float, float), float x1, float x2);
float qgaus(float (*func)(float), float a, float b);
float f1(float x);
float f2(float y);
float f3(float z);
float func(float x, float y, float z);
float yy1(float x);
float yy2(float x);
float z1(float x, float y);
float z2(float x, float y);

int main() {
    float x1 = -1.0, x2 = 1.0;
    float result = quad3d(func, x1, x2);
    printf("Result: %f\n", result);
    return 0;
}

// Function definitions
float quad3d(float (*func)(float, float, float), float x1, float x2) {
    nrfunc = func;
    return qgaus(f1, x1, x2);
}

float f1(float x) {
    xsav = x;
    return qgaus(f2, yy1(x), yy2(x));
}

float f2(float y) {
    ysav = y;
    return qgaus(f3, z1(xsav, y), z2(xsav, y));
}

float f3(float z) {
    return (*nrfunc)(xsav, ysav, z);
}

// 64th degree Legendre-Gauss quadrature integral approximation
// There are only half the points because negative values repeat
// and can be included with a minus sign.
float qgaus(float (*func)(float), float a, float b) {
    int j;
    float xr, xm, dx, s;

    static float x[32] = { 
        0.0243502926634244, 0.0729931217877990, 0.1214628192961206, 0.1696444204239928,
        0.2174236437400071, 0.2646871622087674, 0.3113228719902110, 0.3572201583376681,
        0.4022701579639916, 0.4463660172534641, 0.4894031457070530, 0.5312794640198946,
        0.5718956462026340, 0.6111553551723933, 0.6489654712546573, 0.6852363130542333,
        0.7198818501716109, 0.7528199072605319, 0.7839723589433414, 0.8132653151227975,
        0.8406292962525803, 0.8659993981540928, 0.8893154459951141, 0.9105221370785028,
        0.9295691721319396, 0.9464113748584028, 0.9610087996520538, 0.9733268277899110,
        0.9833362538846260, 0.9910133714767443, 0.9963401167719553, 0.9993050417357722 
    };

    static float w[32] = { 
        0.0486909570091397, 0.0485754674415034, 0.0483447622348030, 0.0479993885964583,
        0.0475401657148303, 0.0469681828162100, 0.0462847965813144, 0.0454916279274181,
        0.0445905581637566, 0.0435837245293235, 0.0424735151236536, 0.0412625632426235,
        0.0399537411327203, 0.0385501531786156, 0.0370551285402400, 0.0354722132568824,
        0.0338051618371416, 0.0320579283548516, 0.0302346570724025, 0.0283396726142595,
        0.0263774697150547, 0.0243527025687109, 0.0222701738083833, 0.0201348231535302,
        0.0179517157756973, 0.0157260304760247, 0.0134630478967186, 0.0111681394601311,
        0.0088467598263639, 0.0065044579689784, 0.0041470332605625, 0.0017832807216964 
    };

    xm = 0.5 * (b + a);
    xr = 0.5 * (b - a);
    s = 0;

    for (j = 0; j < 32; j++) {
        dx = xr * x[j];
        s += w[j] * ((*func)(xm + dx) + (*func)(xm - dx));
    }
    return s * xr;
}

// User defined functions start here
float func(float x, float y, float z)
{
    return x*x + y*y - z*z;
}

// yy1 is defined in the m,ath.h libarary as the Bessel function
// hence why we have to use yy1 instead
float yy1(float x)
{   
    // Lower limit of integration along the y-axis
    return -sqrt(1.0 - x*x);
}

float yy2(float x)
{
    // Upper limit of integration along the y-axis
    return sqrt(1.0 - x*x);
}

float z1(float x, float y)
{
    // Lower limit of integration along the z-axis
    return -sqrt(1.0 - x*x - y*y);
}

float z2(float x, float y)
{
    // Upper limit of integration along the z-axis
    return sqrt(1.0 - x*x - y*y);
}


