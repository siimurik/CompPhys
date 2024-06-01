#include <math.h>
#include <string.h>
#include <stdio.h> // Perform input and output operations
#include <stdlib.h>

#define PI 4.0*atan(1.0)
#define g 9.81

//void matrix_csv(char *filename,double a[][N],int n,int m);

int main(){
    //----------------------------------------
    // Declaration of variables
    double x0, y0, R, x, y, vx, vy, t, tf, dt;
    double theta, v0x, v0y, v0;
    char str[4] = "buf";
    //----------------------------------------
    // Ask user for input:
    printf("# Enter v0, theta (in degrees):\n");
    scanf("%lf, %lf", &v0, &theta);
    printf("# Enter tf, dt:\n");
    scanf("%lf, %lf", &tf, &dt);
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

    //Compute
    t = 0.0;
    while (y <= -0.001) {
        x = v0x * t;
        y = v0y*t - 0.5*g*t*t;
        vx = v0x;
        vy = v0y - g*t;
        printf("%lf\t%lf\t%lf\t%lf\n", x, y, vx, vy);
        t = t + dt;
    }
    return 0;
}

//================================================================================
// Function that creates and stores the temperature matrix values
//================================================================================
/* Example case:
    // Save matrix output to CSV format 
    char str[10] = "colors";
    matrix_csv(str,c1,N,N);
*/
/*
void matrix_csv(char *filename, double a[][N],int n,int m){
    printf("\n Creating %s.csv file for matrix.",filename);
    FILE *fp;
    int i,j;
    filename=strcat(filename,".csv");
    fp=fopen(filename,"w+");
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            fprintf(fp,"%20.14f, ",a[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("\n %s file created.\n",filename);
}
*/