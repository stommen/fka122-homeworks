#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"

double DMC(
           double* walkmen_pos,
           int N_its,
           int N_0,
           double dtau,
           double tau,
           double E_T_start,
           gsl_rng *U,
           double gamma,
           FILE *fp_ET,
           FILE *fp_w
          );

double Mooooooorse(
                   double x
                  );

void wavef_ground_state_anal(
                             double *psi, 
                             double *x, 
                             int n
                            );

double weight(
              double dtau, 
              double ET, 
              double x
             ); 

gsl_rng * init_gsl_rng(
                       int seed
                      ); 

void init_walkmen(
                  double *walkers, 
                  int N_walkers
                 );

int
run(
    int argc,
    char *argv[]
   )
{
    // --------------------------------- Task 1a --------------------------------- //
    int N_0 = 200;
    double dtau = 0.02;
    double tau = 5000;
    int N_its = tau / dtau;
    double E_T_start = 0.5;
    gsl_rng *U = init_gsl_rng(19);
    double gamma = 0.5;
    int its_eq = 20 / dtau;

    double* walkmen_pos = (double*)calloc(N_0 * 100, sizeof(double));
    init_walkmen(walkmen_pos, N_0);

    E_T_start = DMC(walkmen_pos, its_eq, N_0, dtau, tau, E_T_start, U, gamma, NULL, NULL);

    FILE* fp_ET = fopen("data/task_1/ET_Nwalk_non_eq.csv","w");
    FILE* fp_w = fopen("data/task_1/I_was_walkin_in_morse.csv","w");

    DMC(walkmen_pos, N_its, N_0, dtau, tau, E_T_start, U, gamma, fp_ET, fp_w);

    free(walkmen_pos);
    fclose(fp_ET);
    fclose(fp_w);

    // --------------------------------- Task 1b --------------------------------- //
    int N_0 = 1000;
    double dtau = 0.01;
    double tau = 5000;
    int N_its = tau / dtau;
    double E_T_start = 0.5;
    gsl_rng *U = init_gsl_rng(19);
    double gamma = 0.5;
    int its_eq = 20 / dtau;

    // --------------------------------- Task 2a --------------------------------- //

    // --------------------------------- Task 2b --------------------------------- //
    
    return 0;
}

double 
DMC(
    double* walkmen_pos,
    int N_its,
    int N_0,
    double dtau,
    double tau,
    double E_T_start,
    gsl_rng *U,
    double gamma,
    FILE *fp_ET,
    FILE *fp_w
   )
{
    double E_T = E_T_start;
    double E_T_avg = E_T_start;
    int N_sprinters = N_0;

    for(int i = 0; i < N_its ; i++)
    {   
        int *num_walkers = (int *)calloc(N_sprinters, sizeof(int));
        for(int j = 0; j < N_sprinters; j++)
        {   
            int m = (int)(weight(dtau, E_T, walkmen_pos[j]) + gsl_rng_uniform(U) * 1.);
            num_walkers[j] = m;
        }
        // Number of surviving walkers
        int N_survived = int_sum(num_walkers, N_sprinters);

        // Giving birth to new walkers 
        double *walkmen_pos_new = (double *)malloc(N_survived * sizeof(double));

        // Index for the new array
        int M = 0;
        // Populate the new array with the walkers
        for (int k = 0; k < N_sprinters; k++) 
        {
            // If the walker survives (num_walkers[j] > 0)
            for (int m = 0; m < num_walkers[k]; m++) {
                walkmen_pos_new[M] = walkmen_pos[k]; // Copy walker's position
                M++; // Increment the new array index
            }
        }

        // Generate new positions, x' = x + sqrt(dtau)*G
        for(int l = 0; l < N_survived; l++)
        {
            walkmen_pos[l] = walkmen_pos_new[l] + gsl_ran_gaussian(U, 1.) * sqrt(dtau); 
            if (i > 9 * N_its / 10 && fp_w != NULL)
            {
                fprintf(fp_w,"%lf,\n", walkmen_pos[l]);
            }
        }  

        // Updating ET
        double sprinter_ratio = (double)N_survived / (double)N_0;
        E_T = E_T_avg - gamma * log(sprinter_ratio); // E_T at i+1
        E_T_avg = E_T / (i+1) + E_T_avg * i / (i+1); // E_T_avg at i+1

        if (fp_ET != NULL)
        {
            fprintf(fp_ET, "%lf, %i, %f, %f\n", E_T_avg, N_survived, (100. * (double)N_survived) / (double)N_sprinters, E_T);
        }
        
        N_sprinters = N_survived;
        free(walkmen_pos_new);
        free(num_walkers);
    }

    return E_T;
}

void
init_walkmen(
             double* walkers, 
             int N_walkers
            )
{
    int j = 0;
    for(double i = -5; i < 5; i += 10./N_walkers)
    {
        walkers[j] = i;
        j+=1;
    }
}

double
Mooooooorse(
            double x
           )
{
    double V = 0.5 * pow((1  - exp(-x)),2);
    
    return V;
}

double
weight(
       double dtau, 
       double ET, 
       double x
      )
{   
    double W = exp((ET - Mooooooorse(x)) * dtau);
    return W;
}

void
wavef_ground_state_anal(
                        double *psi, 
                        double *x, 
                        int n
                       )
{
    for(int i = 0;i < n; i++)
    {
        psi[i] = sqrt(2.) * exp(-exp(-x[i]) - x[i]/2.);
    }
}

gsl_rng *
init_gsl_rng(
             int seed
            )
{
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
