#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"


double Mooooooorse(double x); void wavef_ground_state_anal(double *Psi,double *x, int n);
double weight(double dtau,double ET, double x); 
gsl_rng * init_gsl_rng(int seed); void init_walkmen(double* walkers, int N_walkers);

int
run(
    int argc,
    char *argv[]
   )
{
    int N0 = 200;
    int N_sprinters = N0;
    // double E0 = 3./8.;
    double dtau = 0.02; double tau = 5000.;
    int N_its = tau / dtau;
    double E_T = 0.5;
    gsl_rng *U = init_gsl_rng(19);
    double gamma = 0.5;
    int its_eq = (tau /20.) / dtau;
    double ET_avg = 0;

    double* walkmen_pos = (double*)calloc(N_sprinters * 100, sizeof(double));
    init_walkmen(walkmen_pos,N_sprinters);

    FILE* fpET = fopen("ET_Nwalk_non_eq.csv","w");
    FILE* fpw = fopen("I_was_walkin_in_morse.csv","w");
    printf("eq %i\n",its_eq);


    for(int i = 0; i < N_its ; i++)
    {
        // printf("N_sprinters %i\n",N_sprinters);
        // calculating m for each walker
        // for(int j = 0; j < N_sprinters;j++)
        // {
        //     printf("walker_new %f \n",walkmen_pos[j]);
        // } 
            
        int* num_walkers = (int*)calloc(N_sprinters, sizeof(int));
        for(int j = 0; j< N_sprinters;j++)
        {
            //if(weight(dtau,E_T,walkmen_pos[j]) >= 1.)
            //{
                //printf("W(x) %i\n",(int)weight(dtau,E_T,walkmen_pos[j]));
            int m = (int)(weight(dtau,E_T,walkmen_pos[j]) + gsl_rng_uniform(U) * 1.);
            num_walkers[j] = m;
            // }
            // else
            // {
            //    num_walkers[j] = 0;
            // }
            //printf("m %i\t W %f\t x %f \n",num_walkers[j],weight(dtau,E_T,walkmen_pos[j]), walkmen_pos[j]);
        }
        // number of surviving walkers
        int N_sprinters_1 = int_sum(num_walkers, N_sprinters);
        //printf("N %i\n",(N_sprinters_1));

        // Giving birth to new walkers 
        double* walkmen_pos_new = (double*)malloc(N_sprinters_1 * sizeof(double));
        //int shift = 0;
        // for(int j = 0; j < N_sprinters;j++)
        // {
        //     if(num_walkers[j] == 0)
        //     {
        //         shift -= 1;
        //     }
        //     else if(num_walkers[j] == 1)
        //     {
        //         for(int m = 0; m < num_walkers[j];m++)
        //         {
        //             walkmen_pos_new[j+shift] = walkmen_pos[j];
        //         } 
        //     }     
        //     else if(num_walkers[j] > 1)
        //     {
        //         for(int m = 0; m < num_walkers[j];m++)
        //         {
        //             walkmen_pos_new[j+shift] = walkmen_pos[j];
        //             shift += 1;
        //         }      
        //     }
        //     //printf("Shift %i\t N_new - N %i\n",shift,N_sprinters_1 - N_sprinters);
        // }
        // Allocate memory for the new walkers
        //double* walkmen_pos_new = (double*)malloc(N_sprinters_1 * sizeof(double));

        // Index for the new array
        int M = 0;
        // Populate the new array with the walkers
        for (int j = 0; j < N_sprinters; j++) 
        {
            // If the walker survives (num_walkers[j] > 0)
            for (int m = 0; m < num_walkers[j]; m++) {
                walkmen_pos_new[M] = walkmen_pos[j]; // Copy walker's position
                M++; // Increment the new array index
            }
        }

        // generate new positions, x' = x + sqrt(dtau)*G

        for(int j = 0; j < N_sprinters_1;j++)
        {
            walkmen_pos[j + M] = 0;
            walkmen_pos[j] = walkmen_pos_new[j] + gsl_ran_gaussian(U, 1.) * sqrt(dtau); 


            if(i > its_eq)
            {
                fprintf(fpw,"%lf,\t",walkmen_pos[j]);
            }
        }  
        if(i > its_eq)
        {
            fprintf(fpw,"\n");
        }

        // Updating ET
        double sprinter_ratio = (double)N_sprinters_1 / (double)N0  ;

        if(i >= its_eq)
        {
            int j = i - its_eq;
            //printf("ET_avg = %f\n",ET_avg);
            E_T = ET_avg - gamma * log(sprinter_ratio);
            ET_avg = (ET_avg * j + E_T) / (j + 1);
            //E_T = ET_avg;
            //E_T -= gamma * log(sprinter_ratio);

            fprintf(fpET,"%lf,%i,%f\n",ET_avg,N_sprinters_1,(100. * (double)N_sprinters_1) / (double)N_sprinters);
        }
        else
        {
            E_T -= gamma * log(sprinter_ratio);
        }
        N_sprinters = N_sprinters_1;
        free(walkmen_pos_new);
        free(num_walkers);
    }
    fclose(fpET);
    fclose(fpw);
    
    return 0;
}

void
init_walkmen(double* walkers, int N_walkers)
{
    int j = 0;
    for(double i = -5; i < 5; i += 10./N_walkers)
    {
        walkers[j] = i;
        j+=1;
    }
}

double
Mooooooorse(double x)
{
    double V = 0.5 * pow((1  - exp(-x)),2);
    
    return V;
}

double
weight(double dtau,double ET, double x)
{
    double W = exp((ET - Mooooooorse(x)) * dtau);
    return W;
}

void
wavef_ground_state_anal(double *Psi,double *x, int n)
{
    for(int i = 0;i < n; i++)
    {
        Psi[i] = sqrt(2.) * exp(-exp(-x[i]) - x[i]/2.);
    }
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

