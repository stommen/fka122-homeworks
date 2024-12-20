#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>

#include "tools.h"

void
elementwise_addition(
                     double *res,
                     double *v1,
                     double *v2,
                     unsigned int len
                    )
{
    for (int i = 0; i < len; i++) {
        res[i] = v1[i] + v2[i];
    }
}

void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          )
{
    for (int i = 0; i < len; i++) {
        res[i] = v1[i] * v2[i];
    }
}

int
int_sum(int *v,int len)
{
    int sum = 0;
    for(int i = 0;i < len; i++)
    {
        sum += v[i];
    }
    return sum;
}

void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len)
{
    for (int i = 0; i < len; i++) {
        res[i] = v[i] + constant;
    }
}

void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len)
{
    for (int i = 0; i < len; i++) {
        res[i] = v[i] * constant;
    }
}

double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           )
{
    double res = 0.;
    for (int i = 0; i < len; i++) {
        res += v1[i] * v2[i];
    }
    return res;
}

double **
create_2D_array(
                unsigned int row_size,
                unsigned int column_size
               )
{
    // Allocate memory for row pointers
    double **array = malloc(row_size * sizeof(double *));
    if (array == NULL) {
        fprintf(stderr, "Memory allocation failed for row pointers.\n");
        return NULL;
    }

    // Allocate a single contiguous block for all elements, i.e. all matrix elements
    // get stored in a single block of memory, array[0] is the start of the block
    // i.e. here we have array[0][index] as the way to access the elements
    array[0] = malloc(row_size * column_size * sizeof(double));
    if (array[0] == NULL) {
        fprintf(stderr, "Memory allocation failed for data block.\n");
        free(array);
        return NULL;
    }

    // Set the row pointers to the appropriate positions in the block, in order to
    // allow for the array[i][j] syntax. Now array[i] points to the start of the
    // i-th row.
    for (int i = 1; i < row_size; i++) {
        array[i] = array[0] + i * column_size;
    }

    return array;
}

void
destroy_2D_array(
                 double **array
                )
{
    if (array != NULL) {
        free(array[0]);  // Free the contiguous block
        free(array);      // Free the row pointers
    }
}

void
matrix_vector_multiplication(
                             double *result,
                             double **A,
                             double *b,
                             unsigned int n,
                             unsigned int m
                            )
{
    for (int i = 0; i < n; i++) {
        result[i] = 0.;
        for (int j = 0; j < m; j++) {
            result[i] += A[i][j] * b[j];
        }
    }
}

void
matrix_matrix_multiplication(
                             double **result,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            )
{
    for (int i = 0; i < n; i++) {
        for (int kappa = 0; kappa < k; kappa++) {
            result[i][kappa] = 0.;
            for (int j = 0; j < m; j++) {
                result[i][kappa] += A[i][j] * B[j][kappa];
            }
        }
    }
}

double
vector_norm(
            double *v1,
            unsigned int len
           )
{   
    double res = 0.;
    for (int i = 0; i < len; i++) {
        res += v1[i] * v1[i];
    }
    return sqrt(res);
}


void
normalize_vector(
                 double *v,
                 unsigned int len
                )
{
    double norm = vector_norm(v, len);
    multiplication_with_constant(v, v, 1. / norm, len);
}

double
average(
        double *v,
        unsigned int len
       )
{
    double res = 0.;
    for (int i = 0; i < len; i++) {
        res += v[i];
    }
    return res / len;
}


double
standard_deviation(
                   double *v,
                   unsigned int len
                  )
{   
    double avg = average(v, len);
    double res = 0.;
    for (int i = 0; i < len; i++) {
        res += (v[i] - avg) * (v[i] - avg);
    }
    return sqrt(res / len);
}

double
variance(
        double *v,
        unsigned int len
        )
{
    double var;
    var = pow(standard_deviation(v, len), 2);

    return var;
}

double
autocorrelation(
			   double *data,
			   int data_len,
			   int time_lag_ind
               )
{   
    double corr = 0.;
    double avg = average(data, data_len);
    double var = variance(data, data_len);

    for (int i = 0; i < data_len - time_lag_ind; i++)
    {
        corr += (data[i]*data[i+time_lag_ind] - avg*avg) / var; 
    }

    return corr / (data_len - time_lag_ind);
}

double 
block_average(
              double *data,
              int data_len,
              int block_size
             )
{
    int num_blocks = data_len / block_size;
    double *block_averages = (double *)calloc(num_blocks, sizeof(double));
    for (int i = 0; i < num_blocks; i++) {
        double sum = 0.0;
        for (int j = 0; j < block_size; j++) {
            sum += data[i * block_size + j];
        }
        block_averages[i] = sum / block_size;
    }
    double block_avg = block_size * variance(block_averages, num_blocks) / variance(data, data_len);
    free(block_averages);

    return block_avg;
}

double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        )
{
    double res = 0.;
    // With previous defined functions
    double *diff = malloc(len * sizeof(double));
    multiplication_with_constant(v1, v1, -1., len);
    elementwise_addition(diff, v1, v2, len);
    res = vector_norm(diff, len);
    free(diff);
    return res;
}

void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int v_len
                      )
{
    double sum = 0.;
    res[0] = 0.;
    for (int i = 1; i < v_len; i++) {
        sum = 0.5 * (v[i - 1] + v[i]) * dx;
        res[i] = res[i - 1] + sum;
    }
}

void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms
         )
{
    fprintf(fp, "%i\nLattice=\"%f 0.0 0.0 0.0 %f 0.0 0.0 0.0 %f\" ", natoms, alat, alat, alat);
    fprintf(fp, "Properties=species:S:1:pos:R:3:vel:R:3 pbc=\"T T T\"\n");
    for(int i = 0; i < natoms; ++i){
        fprintf(fp, "%s %f %f %f %f %f %f\n",
                symbol, positions[i][0], positions[i][1], positions[i][2],
                velocities[i][0], velocities[i][1], velocities[i][2]);
    }
}

void fft_freq(
              double *res,
              int n,
              double timestep
             )
{
    for (int i = 0; i < n; i++) {
        if (i < n / 2) {
            res[i] = 2 * M_PI * i / (n * timestep);
        } 
        else {
            res[i] = 2 * M_PI * (i - n) / (n * timestep);
        }
    }
}

/* Freely given functions */
void
skip_line(
          FILE *fp
         )
{
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
}

void
read_xyz(
         FILE *fp,
         char *symbol,
         double **positions,
         double **velocities,
         double *alat
        )
{
    int natoms;
    if(fscanf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" ", &natoms, alat, alat, alat) == 0){
        perror("Error");
    }
    skip_line(fp);
    for(int i = 0; i < natoms; ++i){
        fscanf(fp, "%s %lf %lf %lf ",
                symbol, &positions[i][0], &positions[i][1], &positions[i][2]);
        fscanf(fp, "%lf %lf %lf\n",
                &velocities[i][0], &velocities[i][1], &velocities[i][2]);
    }
}

void powerspectrum(
                   double *res,
                   double *signal,
                   int n,
                   double timestep
                  )
{
    /* Declaration of variables */
    double *complex_coefficient = malloc(sizeof(double) * 2*n); // array for the complex fft data
    double *data_cp = malloc(sizeof(double) * n);

    /*make copy of data to avoid messing with data in the transform*/
    for (int i = 0; i < n; i++) {
    data_cp[i] = signal[i];
    }

    /* Declare wavetable and workspace for fft */
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;

    /* Allocate space for wavetable and workspace for fft */
    work = gsl_fft_real_workspace_alloc(n);
    real = gsl_fft_real_wavetable_alloc(n);

    /* Do the fft*/
    gsl_fft_real_transform(data_cp, 1, n, real, work);

    /* Unpack the output into array with alternating real and imaginary part */
    gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

    /*fill the output powspec_data with the powerspectrum */
    for (int i = 0; i < n; i++) {
    res[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1]);
    res[i] *= timestep / n;
    }

    /* Free memory of wavetable and workspace */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    free(complex_coefficient);
    free(data_cp);
}
