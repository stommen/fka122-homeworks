#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "lattice.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double PT_fun(double P, double T, double delta_E, double k_B);
void newton_raphson(double *P, double T, double delta_E, double k_B, FILE *file);

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

    double **pos_1 = create_2D_array(1000, 3);
    init_sc(pos_1, 10, 2.949, (double[1][3]){{0, 0, 0}});
    FILE *fp = fopen("pos.csv", "a");
    for (int i = 0; i < 1000; i++)
    {
        fprintf(fp, "%f, %f, %f\n", pos_1[i][0], pos_1[i][1], pos_1[i][2]);
    }

    double **pos_2 = create_2D_array(729, 3);
    init_sc(pos_2, 9, 2.949, (double[1][3]){{0.5, 0.5, 0.5}});
    for (int i = 0; i < 729; i++)
    {
        fprintf(fp, "%f, %f, %f\n", pos_2[i][0], pos_2[i][1], pos_2[i][2]);
    }

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
