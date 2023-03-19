/*========================================
 Compile and excute with:
    $ gcc sum_check.c -o sum_check -lm
    $ ./sum_check
========================================*/
#include <stdio.h>
#include <math.h>
#include <time.h>

int main()
{
    double sum = 0.0, tol = 1e-10;
    int k = 1;
    clock_t start_time = clock(); // start the timer
    double elapsed_time = 0.0;    // initialize the elapsed time to zero
    double term;
    do {
        term = sin(2.0*k)/(k + pow(3.0, k));
        //term = 3.0*k/((k*k-2.0)*(log(k)));
        sum += term;
        k++;
        elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC; // update elapsed time
        if (elapsed_time >= 5) { // stop if elapsed time is >= 5 seconds
            printf("\rElapsed time: %.2f seconds", elapsed_time); // print elapsed time without newline
            printf("\nElapsed time exceeds %.2f seconds.\n", elapsed_time);
            printf("The series is most likely divergent.\n");
            fflush(stdout); // flush stdout to force printing
            break;
        }
        printf("\rElapsed time: %.2f seconds", elapsed_time); // print elapsed time without newline
        fflush(stdout); // flush stdout to force printing
    } while (fabs(term) > tol);
    printf("\nsum = %.16f\n", sum);
    printf("k = %i\n", k);
    return 0;
}
