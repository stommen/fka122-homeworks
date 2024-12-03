#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "lattice.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int
run(
    int argc,
    char *argv[]
   )
{
    // ------------------------------- Constants --------------------------------- //
    // double E_cucu = -0.436; // eV
    // double E_znzn = -0.113; // eV
    // double E_cuzn = -0.294; // eV
    // double k_B = 8.617333262145e-5; // eV/K
    // double N = 2500000;
    double a = 2.949; // [Å]

    // FILE *file = fopen("P_T.csv", "w");
    // double delta_E = E_cucu + E_znzn - 2 * E_cuzn;
    // double *T = (double *)malloc(N * sizeof(double));
    // double *P = (double *)malloc(N * sizeof(double));
    // P[0] = -1.;
    // double step_size = 0.00000001;
    // for (int i = 1; i < N-1; i++)
    // {   
    //     P[i] = P[i-1] + step_size;
    //     if (fabs(P[i]) > 0.99)
    //     {   
    //         step_size = 0.00000001;
    //     }
    //     else if (fabs(P[i]) <= 0.99 )
    //     {
    //         step_size = 0.001;
    //     }
    //     else if (fabs(P[i]) > 1)
    //     {
    //         break;
    //     }
    //     T[i] = 4 * P[i] * delta_E / k_B / log((1 + P[i]) / (1 - P[i]));
    //     fprintf(file, "%f, %f, %f\n", P[i], T[i], step_size);
    // }

    // fclose(file);
    // free(T);
    // free(P);

    double **sub_A = create_2D_array(1000, 3);
    // double **sub_B = create_2D_array(1000, 3);
    int **neighbors = (int **)create_2D_array(1000,8);

    char filename[100];
    sprintf(filename, "neighborsB.csv");
    FILE *fn_B = fopen(filename, "w");
    init_sc(sub_A, 10, a, (double[3]){0, 0, 0}, neighbors,8, fn_B);
    fclose(fn_B);
    //init_sc(sub_B, 10, a, (double[1][3]){{0.5, 0.5, 0.5}});

    return 0;
}

gsl_rng *
init_gsl_rng(int seed){
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default; // default random number generator
    r = gsl_rng_alloc(T); // allocate memory for the random number generator

    if (!r) {
        fprintf(stderr, "Error: Could not allocate memory for RNG.\n");
        exit(EXIT_FAILURE); // Exit if allocation fails
    }

    // Set the seed
    gsl_rng_set(r, seed);

    return r;
}



