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
    int N_sprinters = 200;
    // double E0 = 3./8.;
    double dtau = 0.02; double tau = 350.;
    int N_its = tau / dtau;
    double E_T = 0.5;
    gsl_rng *U = init_gsl_rng(19);
    double gamma = 0.5;

    double* walkmen_pos = (double*)calloc(N_sprinters * 10, sizeof(double));
    init_walkmen(walkmen_pos,N_sprinters);

    FILE* fp = fopen("E_T_non_eq.csv","w");

    for(int i = 0; i < N_its ; i++)
    {
        //printf("N_sprinters %i\n",N_sprinters);
        // calculating m for each walker
        // for(int j = 0; j < N_sprinters;j++)
        // {
        //     printf("walker_new %f \n",walkmen_pos[j]);
        // }  
            
        int* num_walkers = (int*)calloc(N_sprinters, sizeof(int));
        for(int j = 0; j< N_sprinters;j++)
        {
            int m = (int) (weight(dtau,E_T,walkmen_pos[j]) + gsl_rng_uniform(U) * 1.);
            num_walkers[j] = m;
            //printf("m %i\t W %f\t x %f \n",num_walkers[j],weight(dtau,E_T,walkmen_pos[j]), walkmen_pos[j]);
        }
        // number of surviving walkers
        int N_sprinters_1 = int_sum(num_walkers, N_sprinters);
        // printf("Mortality rate %i\n",(100 * N_sprinters_1) / N_sprinters);

        // Giving birth to new walkers 
        double* walkmen_pos_new = (double*)malloc(N_sprinters_1 * sizeof(double));
        int M = 0;
        for(int j = 0; j < N_sprinters;j++)
        {
            if(num_walkers[j] == 0)
            {
                M -= 1;
            }
            else if(num_walkers[j] > 0)
            {
                for(int m = 0; m < num_walkers[j];m++)
                {
                    walkmen_pos_new[j+M] = walkmen_pos[j];
                    M += m;
                }      
            }
            //printf("M %i\t j %i\n",M,j);
        }

        // generate new positions, x' = x - sqrt(tau)*G
        for(int j = 0; j < N_sprinters_1;j++)
        {
            walkmen_pos[j] = walkmen_pos_new[j] + gsl_ran_gaussian(U, 1.) * sqrt(dtau); 
        }   

        free(walkmen_pos_new);

        // Updating ET
        double sprinter_ratio = (double)N_sprinters_1 / (double)N_sprinters;
        E_T -= gamma * log(sprinter_ratio);
        
        // N_sprinters = N_sprinters_1;
        fprintf(fp,"%lf,%i\n",E_T,N_sprinters_1);
    }
    fclose(fp);

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
    double W = exp(-(Mooooooorse(x) - ET) * dtau);
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

