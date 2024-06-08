#include <math.h>
#include <string.h>
#include <stdio.h> // Perform input and output operations
#include <stdlib.h>

#define PI 4.0*atan(1.0)
#define g 9.81

void write_to_file(FILE *file, double t, double x, double y, double vx, double vy) {
    // Set precision for each variable using %.17g for high precision
    fprintf(file, "%.17g,\t%.17g,\t%.17g,\t%.17g,\t%.17g\n", t, x, y, vx, vy);
}

int main(){
    //----------------------------------------
    // Declaration of variables
    double x0, y0, R, x, y, vx, vy, t, t0, tf, dt;
    double theta, v0x, v0y, v0;
    char str[4] = "buf";
    //----------------------------------------
    // Ask user for input:
    printf("# Enter initial position (x0 & y0) values:\n");
    scanf("%lf, %lf", &x0, &y0);
    printf("# Enter v0, theta (in degrees):\n");
    scanf("%lf, %lf", &v0, &theta);
    printf("# Enter t0, tf, dt:\n");
    scanf("%lf, %lf, %lf", &t0, &tf, &dt);
    printf("# v0  = %.3f", v0);
    printf("\ttheta = %.3f degrees\n", theta);
    printf("# t0  = %6.3f", 0.0);
    printf("\tdt    = %.3f", dt);
    printf("\ttf = %.3f\n", tf);
    //----------------------------------------
    // Initialize
    if (v0 <= 0.0) {
        fprintf(stderr, "Illegal value of v0 <= 0\n");
        exit(1);
    }
    if (theta <= 0.0) {
        fprintf(stderr, "Illegal value of theta <= 0\n");
        exit(1);
    }
    if (theta >= 90.0) {
        fprintf(stderr, "Illegal value of theta >= 90.0\n");
        exit(1);
    }
    theta = (PI/180.0)*theta; //convert to radians
    v0x   = v0*cos(theta);
    v0y   = v0*sin(theta);
    printf("# v0x = %9.6lf\tv0y = %lf\n", v0x, v0y);

    // Set initial position
    // This is necessary bc if not defined before the 
    // 'y' variable in the while loop can hold any 
    // arbitrary value, which can mess up the condition 
    // for the while to start properly.
    x = x0;
    y = y0;

    // Open file for writing
    FILE *file = fopen("output.csv", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        exit(1);
    }

    //Compute
    t = t0;
    while (y >= 0) { //'>=' meaning while the projectile is above ground
        x = v0x * t;
        y = v0y*t - 0.5*g*t*t;
        vx = v0x;
        vy = v0y - g*t;
        printf("%lf\t%lf\t%lf\t%lf\t%lf\n", t, x, y, vx, vy);
        //printf("%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n", t, x, y, vx, vy); // Set precision for console output
        write_to_file(file, t, x, y, vx, vy); // Write to file with precision
        t = t + dt;
    }

    fclose(file);
    return 0;
}