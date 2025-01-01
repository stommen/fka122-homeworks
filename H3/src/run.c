#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"

double DMC(
           double **walkmen_pos,
           int N_its,
           int N_0,
           double dtau,
           double tau,
           double E_T_start,
           gsl_rng *U,
           double gamma,
           FILE *fp_ET,
           FILE *fp_w,
           int ndim
          );

double mooooooorse(
                   double x
                  );

void wavef_ground_state_anal(
                             double *psi, 
                             double *x, 
                             int n
                            );

double weight(
              double dtau, 
              double E_T, 
              double V
             ); 

gsl_rng * init_gsl_rng(
                       int seed
                      ); 

void init_walkmen_1D(
                     double **walkers, 
                     int N_walkers
                    );

void init_walkmen_6D(
                     double **walkers,
                     int N_walkers,
                     gsl_rng *U
                    );

int
run(
    int argc,
    char *argv[]
   )
{
    // --------------------------------- Task 1a --------------------------------- //
    // int N_0 = 200;
    // double dtau = 0.02;
    // double tau = 5000;
    // int N_its = tau / dtau;
    // double E_T_start = 0.5;
    // gsl_rng *U = init_gsl_rng(19);
    // double gamma = 0.5;
    // int its_eq = 20 / dtau;

    // double **walkmen_pos = create_2D_array(N_0 * 100, 1);
    // init_walkmen_1D(walkmen_pos, N_0);

    // E_T_start = DMC(walkmen_pos, its_eq, N_0, dtau, tau, E_T_start, U, gamma, NULL, NULL, 1);

    // FILE* fp_ET = fopen("data/task_1/1D/ET_Nwalk_non_eq.csv", "w");
    // FILE* fp_w = fopen("data/task_1/1D/I_was_walkin_in_morse.csv", "w");

    // DMC(walkmen_pos, N_its, N_0, dtau, tau, E_T_start, U, gamma, fp_ET, fp_w, 1);

    // free(walkmen_pos);
    // fclose(fp_ET);
    // fclose(fp_w);

    // --------------------------------- Task 1b --------------------------------- //
    int N_0 = 1000;
    double dtau = 0.01;
    double tau = 1000;
    int N_its = tau / dtau;
    double E_T_start = -3;
    gsl_rng *U = init_gsl_rng(19);
    double gamma = 0.5;
    int its_eq = 100 / dtau;

    double **walkmen_pos = create_2D_array(N_0 * 10, 6);
    init_walkmen_6D(walkmen_pos, N_0, U);

    E_T_start = DMC(walkmen_pos, its_eq, N_0, dtau, tau, E_T_start, U, gamma, NULL, NULL, 6);

    FILE* fp_ET = fopen("data/task_1/6D/ET_Nwalk_non_eq.csv", "w");
    // FILE* fp_w = fopen("data/task_1/6D/I_was_walkin_in_morse.csv", "w");

    DMC(walkmen_pos, N_its, N_0, dtau, tau, E_T_start, U, gamma, fp_ET, NULL, 6);

    free(walkmen_pos);
    fclose(fp_ET);
    // fclose(fp_w);

    // --------------------------------- Task 2a --------------------------------- //

    // --------------------------------- Task 2b --------------------------------- //
    
    return 0;
}

double 
DMC(
    double **walkmen_pos,
    int N_its,
    int N_0,
    double dtau,
    double tau,
    double E_T_start,
    gsl_rng *U,
    double gamma,
    FILE *fp_ET,
    FILE *fp_w,
    int ndim
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
            if (ndim == 1)
            {
                int m = (int)(weight(dtau, E_T, mooooooorse(walkmen_pos[j][0])) + gsl_rng_uniform(U) * 1.);
                num_walkers[j] = m;
            }
            else if (ndim == 6)
            {
                double r_1 = walkmen_pos[j][0];
                double r_2 = walkmen_pos[j][3];
                double r_12 = sqrt(pow(r_1, 2) + pow(r_2, 2) - 2 * r_1 * r_2 * cos(walkmen_pos[j][2] - walkmen_pos[j][5]));
                double V_tot = - 2 / r_1 - 2 / r_2 + 1 / r_12;
                int m = (int)(weight(dtau, E_T, V_tot) + gsl_rng_uniform(U) * 1.);
                num_walkers[j] = m;
            }
        }
        // Number of surviving walkers
        int N_survived = int_sum(num_walkers, N_sprinters);

        // Giving birth to new walkers 
        double **walkmen_pos_new = create_2D_array(N_survived, ndim);

        // Index for the new array
        int M = 0;
        // Populate the new array with the walkers
        for (int k = 0; k < N_sprinters; k++) 
        {
            // If the walker survives (num_walkers[j] > 0)
            for (int m = 0; m < num_walkers[k]; m++)
            {   
                for(int n = 0; n < ndim; n++)
                {
                    walkmen_pos_new[M][n] = walkmen_pos[k][n]; // Copy walker's position
                }
                M++; // Increment the new array index
            }
        }

        // Generate new positions, x' = x + sqrt(dtau)*G
        for(int l = 0; l < N_survived; l++)
        {
            for (int n = 0; n < ndim; n++)
            {   
                walkmen_pos[l][n] = walkmen_pos_new[l][n] + gsl_ran_gaussian(U, 1.) * sqrt(dtau);
                if (i > 9 * N_its / 10 && fp_w != NULL)
                {   
                    if (ndim == 1)
                    {
                        fprintf(fp_w, "%lf\n", walkmen_pos[l][0]);
                    }
                    else
                    {   
                        if (n < ndim - 1)
                        {
                            fprintf(fp_w, "%lf, ", walkmen_pos[l][n]);
                        }
                        else
                        {
                            fprintf(fp_w, "%lf\n", walkmen_pos[l][n]);
                        }
                    }
                }
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
        destroy_2D_array(walkmen_pos_new);
        free(num_walkers);
    }

    return E_T;
}

void
init_walkmen_1D(
             double **walkers, 
             int N_walkers
            )
{
    int j = 0;
    for (double i = -5; i < 5; i+=10./N_walkers)
    {
        walkers[j][0] = i;
        j+=1;
    }
}

void
init_walkmen_6D(
                double **walkers,
                int N_walkers,
                gsl_rng *U
               )
{
    for(int i = 0; i < N_walkers; i++)
    {
        walkers[i][0] = 0.7 + gsl_rng_uniform(U);
        walkers[i][1] = acos(2 * gsl_rng_uniform(U) - 1);
        walkers[i][2] = 2 * M_PI * gsl_rng_uniform(U);
        walkers[i][3] = 0.7 + gsl_rng_uniform(U);
        walkers[i][4] = acos(2 * gsl_rng_uniform(U) - 1);
        walkers[i][5] = 2 * M_PI * gsl_rng_uniform(U);
    }
}

double
mooooooorse(
            double x
           )
{
    double V = 0.5 * pow((1  - exp(-x)),2);
    
    return V;
}

double
weight(
       double dtau, 
       double E_T, 
       double V
      )
{   
    double W = exp((E_T - V) * dtau);

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
