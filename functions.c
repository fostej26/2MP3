#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix)
{
    FILE *file = fopen(filename, "r"); // open file in read mode
    char line[1024];
    
    //Check for errors opening the file
    if (file != NULL){

    //moves scanf pointer to first relevant line of code
    do{
        fgets(line, sizeof(line), file);         
    } while (line[0] == '%' || line[0] == '\n'); 

    //scans first line containing dimensions and nnz
    sscanf(line, "%d %d %d", &(matrix->num_rows), &(matrix->num_cols), &(matrix->num_non_zeros));

    //assume symmetric matrix and fill in rest of vals with symmetry
    matrix->num_non_zeros += (matrix->num_non_zeros - matrix->num_cols);
    
    //allocate memory for matrix properties & temporary row array
    int *temp_row = (int *)malloc((matrix->num_non_zeros) * sizeof(int)); 
    matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
    matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));

    //ensure memory is properly allocated, exit if failure
    if (temp_row == NULL || matrix->col_ind == NULL || matrix->csr_data == NULL){
        printf("Unable to allocate memory.\n");
        exit(1); 
    }

    //initialize temp vars
    int row, col;
    double val;


    for (int i = 0; i < matrix->num_non_zeros; ++i){

        //scan a line in the file for each entry
        fscanf(file, "%d %d %lf\n", &row, &col, &val);

        //ensure that col_ind and row ind are in index notation not cartesian notation
        matrix->csr_data[i] = val;
        matrix->col_ind[i] = col - 1;
        temp_row[i] = row - 1;


        //if an entry is not along main diagonal, fill in missing entries with symmetry: AI
        if ((matrix->col_ind[i] != temp_row[i])){
            
            temp_row[i + 1] = matrix->col_ind[i];
            matrix->col_ind[i + 1] = temp_row[i];
            matrix->csr_data[i + 1] = matrix->csr_data[i];
            i++;
        }
    }


    for (int i = 0; i < matrix->num_non_zeros - 1; i++) {
        for (int j = 0; j < matrix->num_non_zeros - i - 1; j++) {
        if (temp_row[j] > temp_row[j + 1]) {
            //swap elements in temp_row
            int temp = temp_row[j];
            temp_row[j] = temp_row[j + 1];
            temp_row[j + 1] = temp;

            //swap corresponding elements in matrix->col_ind and csr_vals
            temp = matrix->col_ind[j];
            matrix->col_ind[j] = matrix->col_ind[j + 1];
            matrix->col_ind[j + 1] = temp;

            double temp_double = matrix->csr_data[j];
            matrix->csr_data[j] = matrix->csr_data[j + 1];
            matrix->csr_data[j + 1] = temp_double;
        }
        }
    }

    // this for loop converts the rows array to the proper array
    int ind_count = 1;
    int val_count = 1;
    for (int i = 1; i < matrix->num_non_zeros; i++)
    {
        if (temp_row[i] != temp_row[i - 1])
        {
            temp_row[ind_count] = i;
            ind_count++;
            val_count++;
        }
    }
    val_count++;
    //adds one more value for final value of row ptr to display row_ptrs and nnz
    temp_row[ind_count] = matrix->num_non_zeros;

    //allocates row_ptr to n+1 size (val_count = n+1)
    matrix->row_ptr = (int *)malloc((val_count) * sizeof(int));
    for (int i = 0; i < val_count; i++){
        matrix->row_ptr[i] = temp_row[i];
    }
    
    free(temp_row);
    fclose(file);
    // printf("\nCSR data:\n");

    // for(int i =0; i<matrix->num_non_zeros; i++){
    //     printf("%.3f \n", matrix->csr_data[i]);
    // }
    // printf("col_ind:\n");
    // for(int i =0; i<matrix->num_non_zeros; i++){
    //     printf("%d ", matrix->col_ind[i]);
    // }
    // printf("\n\nrow_ptr:\n");

    // for(int i =0; i<val_count; i++){
    //     printf("%d ", matrix->row_ptr[i]);
    // }
    // printf("\n\n");
    
    //return the original number of non zeros to be printed with mtx details
    matrix->num_non_zeros -= (matrix->num_non_zeros/2 - matrix->num_cols/2);
    

    }else{
        printf("Error opening file");
        exit(1);
    }
}

void spmv_csr(const CSRMatrix *A, const double *x, double *y) {

    //initialize resultant vector to 0
    for (int i = 0; i < A->num_rows; i++){
         y[i]=0.0; 
     }
    for (int i = 0; i < A->num_rows; i++){ 
        //initialize start and end indicies
        int start_row = A->row_ptr[i]; 
        int end_row = A->row_ptr[i+1];
        //multiply and add to cumulative sum
        for (int j = start_row; j < end_row; j++){
            y[i] += A->csr_data[j] *x[A->col_ind[j]]; //sum products of corresponding rows/columns

        }
    }

}

void solver(const CSRMatrix *A, const double *b, double *x, int max_iterations, double tolerance) {

    int n = A->num_rows;
    int iterations = 0;
    double max_change;

    do {
        max_change = 0;

        for (int i = 0; i < n; i++) {
            double sum = 0;
            double diagonal = 0;
            int start_row = A->row_ptr[i];
            int end_row = A->row_ptr[i + 1];

            //for every entry in a row, check not on main diag.
            for (int j = start_row; j < end_row; j++) {
                if (A->col_ind[j] == i) {
                    diagonal = A->csr_data[j];
                //otherwise add to cumulative sum
                } else {
                    sum += A->csr_data[j] * x[A->col_ind[j]];
                }
                double updated_value = (b[i] - sum) / diagonal; //since x1 = 2x2 + 3x3.../a1

            // Max change between old and new solution
            max_change = fmax(max_change, fabs(updated_value - x[i]));

            // Use the updated value immediately in the same iteration
            x[i] = updated_value;

            }

            
        }
        iterations++;
        //if any change in xn < tolerance, stop the program
    } while (max_change > tolerance && iterations < max_iterations);

    // Optionally, print or use the final solution
    for (int i = 0; i < n; i++) {
        // printf("%f \n", x[i]);
    }
}

void compute_residual(const CSRMatrix *A, const double *b, const double *x, double *residual) {
    // Compute the residual: residual = b - Ax
    for (int i = 0; i < A->num_rows; i++) {
        residual[i] = b[i];
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
            int col = A->col_ind[j];
            residual[i] -= A->csr_data[j] * x[col];
        }
    }
}

double compute_norm(const double *vector, int size) {

    //euclidian norm: sqrt(a[0] + a[1]... a[n])
    double norm = 0.0;
    for (int i = 0; i < size; ++i) {
        norm += vector[i] * vector[i];
    }
    return sqrt(norm);
}

