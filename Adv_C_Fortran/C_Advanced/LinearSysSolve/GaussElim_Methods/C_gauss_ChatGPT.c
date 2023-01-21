/*
 Compile and execute with:
	$ gcc C_gauss_ChatGPT.c -o gauss
	$ ./gauss
*/
#include <stdio.h>

#define N 3 // number of equations
#define M 4 // number of unknowns

// function prototypes
void print_matrix(double a[N][M]);
void gauss_elimination(double a[N][M]);
void back_substitution(double a[N][M]);

int main()
{
    static double a[N][M] = {
        {3.0,    2.0,    4.0,   4.0},
        {2.0,   -3.0,    1.0,   2.0},
        {1.0,    1.0,    2.0,   3.0},
    };

    // print the original matrix
    printf("Original matrix:\n");
    print_matrix(a);

    // apply Gaussian elimination
    gauss_elimination(a);

    // print the resulting matrix
    printf("\nAfter Gaussian elimination:\n");
    print_matrix(a);

    // print the resulting vector
    back_substitution(a);
    return 0;
}

// function to print a matrix
void print_matrix(double a[N][M])
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            printf(" %5.4f ", a[i][j]);
        }
        printf("\n");
    }
}

// function to apply Gaussian elimination
void gauss_elimination(double a[N][M]){
    int i, j, k;
    double p;
    
    // loop over rows
    for(k = 0; k <= N; k++){
        // loop over columns
	    for(i = k+1; i < M; i++){
            // calculate the multiplier for each row
	        p = a[i][k]/a[k][k];
            // loop over columns again
            for(j = k; j <= M; j++){
                // subtract the appropriate multiple of the current row
                // from each of the other rows
                a[i][j] = a[i][j] - (p * a[k][j]);
		    //printf("\n a[%d][%d] = %f",i,j,a[i][j]);
	        }
	    }
    }
}

void back_substitution(double a[N][M]){
    int i, j, k;
    static double x[N];
    double s;

    x[N-1] = a[N-1][N] / a[N-1][N-1];
    for(i = N-2; i >= 0; i--){
        s=0;
        for(j = i+1; j < N; j++){
            s += (a[i][j]*x[j]);
            x[i] = (a[i][N]-s)/a[i][i];
        }
    }
    // Printing
    printf("\nThe result is:");
    for(i = 0; i<N ; i++){
        printf("\nx[%d] = %5.4f ", i+1, x[i]);
    }
    printf("\n\n");
}

// function to apply Gaussian elimination
/*void gauss_elimination(double a[N][M])
{
    int i, j, k;
    double m;

    // loop over rows
    for (k = 0; k < N-1; k++)
    {
        // loop over columns
        for (i = k+1; i < M; i++)
        {
            // calculate the multiplier for each row
            m = a[i][k] / a[k][k];

            // loop over columns again
            for (j = k; j <= M; j++)
            {
                // subtract the appropriate multiple of the current row
                // from each of the other rows
                a[i][j] = a[i][j] - (m * a[k][j]);
            }
        }
    }
}*/

