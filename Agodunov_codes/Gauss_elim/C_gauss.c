//============================================================
// Compile and execute with:
//  > gcc C_gauss.c -o c_gauss
//  > ./c_gauss
//============================================================
// A system of linear equations A*x = b
// Method: calls the basic elimination
//------------------------------------------------------------
// Matrix format for A*x = b:
// (Example for 3x3 matrices)
// [a[1][1] a[1][2] a[1][3]]   [x[1]]   [b[1]]
// [a[2][1] a[2][2] a[2][3]] * [x[2]] = [b[2]]
// [a[3][1] a[3][2] a[3][3]]   [x[3]]   [b[3]]
//------------------------------------------------------------
// Instructions:
// INPUT:
// * The co-efficients of matrix A and vector b
//------------------------------------------------------------
// OUTPUT:
// * Answers in vector x
//============================================================
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define   SIZE   10
int main(void) {
    int i,j,k,n;
    float a[SIZE][SIZE],x[SIZE];
    double s,p;
//============================================================
    printf("======================================================");
    printf("\n        Gauss elimination for the format Ax=b");
    printf("\n======================================================");
    printf("\n\nEnter the number of equations : ");
    scanf("%d",&n);
    printf("\nEnter the co-efficients of the matrix A:\n\n");
// Reading in the matrix A coefficients
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%f", &a[i][j]);
        }
    }
// Reading in the vector b coefficients
    printf("\nEnter the co-efficients of the vector b:\n\n");
    for (int i = 0; i < n; i++) {
        scanf("%f",&a[i][n]) ;
    }
// Gaussian elimination
    clock_t begin = clock();
    for(k=0; k<=n-1; k++){
	    for(i=k+1; i<n; i++){
	        p = a[i][k]/a[k][k];
            for(j=k; j<=n; j++){
		        a[i][j] = a[i][j] - (p * a[k][j]);
		    //printf("\n a[%d][%d] = %f",i,j,a[i][j]);
	        }
	    }
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
// Printing out results after elimination
    printf("\nMatrix A and vector b after elimination:\n");	
	for(i = 0; i < n; i++){
		printf("[");
		for(j = 0; j < n+1; j++){
			printf(" %f ", a[i][j]);
		}
		printf("]\n");
	}
	printf("\n");
// Final results
    x[n-1] = a[n-1][n] / a[n-1][n-1];
    for(i=n-2; i>=0; i--){
        s=0;
        for(j=i+1; j<n; j++){
            s += (a[i][j]*x[j]);
            x[i] = (a[i][n]-s)/a[i][i];
        }
    }
// Printing
    printf("\nThe result is:");
    for(i=0; i<n ; i++){
        printf("\nx[%d] = %.4f",i+1,x[i]);
    }
    printf ("\n\nCalculation time = %.3e sec.\n", time_spent);
    printf("\n");
    return 0;
}