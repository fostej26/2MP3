#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    //Handle the inputs^^^
    const char *filename = argv[1];

    CSRMatrix A;
    ReadMMtoCSR(filename, &A);

    //Initialize all the vector b (in Ax= b)
    double *b = (double *)malloc(A.num_rows * sizeof(double));
    //Set all elements of b to 1
    for (int i = 0; i < A.num_rows; ++i) 
    {
        b[i] = 1.0;
    }

    double *x = (double *)malloc(A.num_rows * sizeof(double)); 
    int max_iterations = 1000;
    double tolerance = 1e-12;

    //Start solver
    clock_t start_time = clock();

    // Call your solver function
    solver(&A, b, x, max_iterations, tolerance);

    //Finish solver time
    clock_t end_time = clock();
    double cpu_time= ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    //Print matrix data
    printf("Matrix name: %s\n", filename);
    printf("Dimension of the matrix: %d by %d\n", A.num_rows, A.num_cols);
    printf("Number of non-zero entries: %d\n", A.num_non_zeros);

    //Print solver time
    printf("CPU time for solver: %f\n", cpu_time);

    //Compute residual and its norm
    double *residual = (double *)malloc(A.num_rows * sizeof(double));
    compute_residual(&A, b, x, residual);
    double residual_norm = compute_norm(residual, A.num_rows);

    printf("Residual Norm: %.16f\n", residual_norm);

    // Free allocated memory
    free(A.col_ind);
    free(A.csr_data);
    free(A.row_ptr);
    free(b);
    free(x);
    free(residual);
   
    return EXIT_SUCCESS;
}