#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "lattice.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng * init_gsl_rng(int seed);
double boundary_distance_between_vectors(double *v1, double *v2, int dim, double box_length);
void nearest_neighbors_bcc(double **pos_A, double **pos_B, int N_A, int N_B, 
                            double box_length, double cutoff, double closest_distance, 
                            FILE *fp_A, FILE *fp_B);

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
    double closest_distance_bcc = sqrt(3) * a / 2;

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

    int L = 10;
    int N_atoms = 2 * L * L * L;
    double **sub_A = create_2D_array(N_atoms / 2, 3);
    double **sub_B = create_2D_array(N_atoms / 2, 3);

    char filename[100];
    sprintf(filename, "data/task_2/neighborsA.csv");
    FILE *fp_NN_A = fopen(filename, "w");
    sprintf(filename, "data/task_2/neighborsB.csv");
    FILE *fp_NN_B = fopen(filename, "w");

    init_sc(sub_A, L, a, (double[3]){0, 0, 0});
    init_sc(sub_B, L, a, (double[3]){0.5, 0.5, 0.5});
    // format filename with variable L
    sprintf(filename, "data/task_2/sub_A_%i.csv", L);
    FILE *fp_A = fopen(filename, "w");
    sprintf(filename, "data/task_2/sub_B_%i.csv", L);
    FILE *fp_B = fopen(filename, "w");
    for (int i = 0; i < N_atoms / 2; i++)
    {
        fprintf(fp_A, "%f, %f, %f\n", sub_A[i][0], sub_A[i][1], sub_A[i][2]);
        fprintf(fp_B, "%f, %f, %f\n", sub_B[i][0], sub_B[i][1], sub_B[i][2]);
    }

    nearest_neighbors_bcc(sub_A, sub_B, 1000, 1000, 10 * a, 0.1, closest_distance_bcc, fp_NN_A, fp_NN_B);
    fclose(fp_NN_A);
    fclose(fp_NN_B);
    destroy_2D_array(sub_A);
    destroy_2D_array(sub_B);

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

double
boundary_distance_between_vectors(double *v1, double *v2, int dim, double box_length)
{
    double r = 0.0;
    double delta;
    for (int d = 0; d < dim; d++) 
    {
        delta = v1[d] - v2[d];
        // Apply boundary conditions
        delta -= round(delta / box_length) * box_length;
        r += delta * delta;
    }

    return sqrt(r);
}

void
nearest_neighbors_bcc(double **pos_A, double **pos_B, int N_A, int N_B, 
                        double box_length, double cutoff, double closest_distance, 
                        FILE *fp_A, FILE *fp_B)
{
    double r;
    for (int i = 0; i < N_A; i++)
    {   
        int count = 0;
        for (int j = 0; j < N_B; j++)
        {
            r = boundary_distance_between_vectors(pos_A[i], pos_B[j], 3, box_length);
            if (fabs(closest_distance - r) < cutoff)
            {   
                count++;
                if (count < 8)
                {
                    fprintf(fp_A, "%i, ", j);
                }
                else if (count == 8)
                {
                    fprintf(fp_A, "%i\n", j);
                }
            }
        }
    }
    for (int i = 0; i < N_B; i++)
    {   
        int count = 0;
        for (int j = 0; j < N_A; j++)
        {
            r = boundary_distance_between_vectors(pos_B[i], pos_A[j], 3, box_length);
            if (fabs(closest_distance - r) < cutoff)
            {
                count++;
                if (count < 8)
                {
                    fprintf(fp_B, "%i, ", j);
                }
                else if (count == 8)
                {
                    fprintf(fp_B, "%i\n", j);
                }
            }
        }
        fprintf(fp_B, "\n");
    }
}
