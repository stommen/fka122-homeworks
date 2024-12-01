#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "run.h"
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
    double E_cucu = -0.436; // eV
    double E_znzn = -0.113; // eV
    double E_cuzn = -0.294; // eV
    double k_B = 8.617333262145e-5; // eV/K

    double *T = (double *)malloc(500 * sizeof(double));
    for (int i = 0; i < 500; i++)
    {
        T[i] = i * 2;
    }
    FILE *file = fopen("P_T.csv", "w");
    double delta_E = E_cucu + E_znzn - 2 * E_cuzn;
    for (int i = 0; i < 500; i++)
    {   
        double P = 0.5;
        newton_raphson(&P, T[i], delta_E, k_B, file);
    }

    fclose(file);
    free(T);

    return 0;
}

double PT_fun(double P, double T, double delta_E, double k_B)
{
    return -4. * P * delta_E + k_B * T *(log(1 + P) - log(1 - P));
}

void newton_raphson(double *P, double T, double delta_E, double k_B, FILE *file)
{
    double tol = 1e-6;
    double f = PT_fun(*P, T, delta_E, k_B);
    double df = -4. * delta_E + k_B * T * (1 / (1 + *P) + 1 / (1 - *P));
    double delta_P = f / df;
    while (fabs(delta_P) > tol)
    {
        *P -= delta_P;
        f = PT_fun(*P, T, delta_E, k_B);
        df = -4. * delta_E + k_B * T * (1 / (1 + *P) + 1 / (1 - *P));
        delta_P = f / df;
    }
    fprintf(file, "%f, %f\n", T, *P);
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