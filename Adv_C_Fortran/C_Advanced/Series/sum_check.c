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
    double sum = 0.0, term = 1.0, tol = 1e-10;
    int k = 4;
    clock_t start_time = clock(); // start the timer
    double elapsed_time = 0.0;    // initialize the elapsed time to zero
    while (fabs(term) > tol) {
        sum += term;
        k++;
        //term = 1.0 / log(k); // compute the next term
        term = pow(2, k) / tgamma(k+1.0); // converges quite nicely
        // term = 1.0/( k*log(k)*pow(log(log(k)), 0.5) );

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
    }
    printf("\nsum = %.16f\n", sum);
    printf("k = %i\n", k);
    return 0;
}
