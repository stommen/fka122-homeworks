#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "lattice.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct
{
    int alpha;
    int beta;
    int valueA;
    int valueB;
} idx;

typedef struct
{
    int accepted;
    double Etot;
} metro;

gsl_rng * init_gsl_rng(int seed);
double boundary_distance_between_vectors(double *v1, double *v2, int dim, double box_length);
void nearest_neighbors_bcc(double **pos_A, double **pos_B, int N_atoms,
                           int **neighbors, double box_length,
                           double cutoff, double closest_distance);
metro metropolis(int its, int *N, int **neighbors, double k_B, double T,
                double E_cucu, double E_znzn, double E_cuzn, double E_tot, gsl_rng *r, FILE *fp);
idx swappie(int *N, gsl_rng *r);
double energy_bond(idx index, int *N, int **neighbors,
                   double E_cucu, double E_znzn, double E_cuzn);
void lattice_to_files(FILE *fp_atoms, FILE *fp_neighbors, int *N, int **neighbors, int N_atoms);

int
run(int argc, char *argv[])
{
    // ------------------------------- Constants --------------------------------- //
    double E_cucu = -0.436; // eV
    double E_znzn = -0.113; // eV
    double E_cuzn = -0.294; // eV
    double k_B = 8.617333262145e-5; // eV/K
    // double N = 2500000;
    double a = 2.949; // [Å]
    double closest_distance_bcc = sqrt(3) * a / 2;
    char filename[100];

    // ------------------------------- Task 1 --------------------------------- //

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

    // ------------------------------- Task 2 --------------------------------- //

    int L = 10;
    int N_atoms = 2 * L * L * L;
    int its_eq = 100000;
    double T = 1000;
    double E_tot = 2000*E_cuzn;
    gsl_rng *r = init_gsl_rng(1234);
    metro metro_result;

    double **sub_A = create_2D_array(N_atoms / 2, 3);
    double **sub_B = create_2D_array(N_atoms / 2, 3);
    int **neighbors = (int **)create_2D_array(N_atoms, 8);
    int *atoms = (int *)calloc(N_atoms, sizeof(int));
    for (int i = 0; i < N_atoms / 2; i++)
    {
        atoms[i] = 1;
    }

    init_sc(sub_A, L, a, (double[3]){0, 0, 0});
    init_sc(sub_B, L, a, (double[3]){0.5, 0.5, 0.5});
    nearest_neighbors_bcc(sub_A, sub_B, N_atoms, neighbors, 10 * a, 0.001, closest_distance_bcc);

    sprintf(filename, "data/task_2/equilibrium_%i_%.0f.csv", its_eq, T);
    FILE *fp_eq = fopen(filename, "w");
    fprintf(fp_eq, "accepted, E_tot\n");

    metro_result = metropolis(its_eq, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_tot, r, fp_eq);
    E_tot = metro_result.Etot;
    int accepted = metro_result.accepted;
    printf("Acceptance rate from equilibrium: %f\n", (double)accepted / its_eq);

    int its = 1000000;
    sprintf(filename, "data/task_2/energy_%i_%i_%.0f.csv", its_eq, its, T);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "Equilibrium iterations, Accepted ratio equilibrium, E_tot equilibrium\n");
    fprintf(fp, "%i, %f, %f\n", its_eq, (double)accepted / its_eq, E_tot);
    fprintf(fp, "accepted, E_tot\n");

    sprintf(filename, "data/task_2/lattice/atoms_%i_%i_%.0f.csv", its_eq, its, T);
    FILE *fp_atoms = fopen(filename, "w");
    sprintf(filename, "data/task_2/lattice/neighbors.csv");
    FILE *fp_neighbors = fopen(filename, "w");

    metro_result = metropolis(its, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_tot, r, fp);

    lattice_to_files(fp_atoms, fp_neighbors, atoms, neighbors, N_atoms);

    destroy_2D_array(sub_A);
    destroy_2D_array(sub_B);
    destroy_2D_array((double **)neighbors);
    free(atoms);
    fclose(fp);
    fclose(fp_atoms);
    gsl_rng_free(r);

    // ------------------------------- Task 3 --------------------------------- //

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
nearest_neighbors_bcc(double **pos_A, double **pos_B, int N_atoms, 
                      int **neighbors, double box_length, 
                      double cutoff, double closest_distance)
{
    double r;
    for (int i = 0; i < N_atoms / 2; i++)
    {   
        int count = 0;
        for (int j = 0; j < N_atoms / 2; j++)
        {
            r = boundary_distance_between_vectors(pos_A[i], pos_B[j], 3, box_length);
            if (fabs(closest_distance - r) < cutoff)
            {      
                neighbors[i][count] = j + 1000;
                count++;
            }
        }
    }
    for (int i = 0; i < N_atoms / 2; i++)
    {   
        int count = 0;
        for (int j = 0; j < N_atoms / 2; j++)
        {
            r = boundary_distance_between_vectors(pos_B[i], pos_A[j], 3, box_length);
            if (fabs(closest_distance - r) < cutoff)
            {   
                neighbors[i+1000][count] = j;
                count++;
            }
        }
    }
}

metro
metropolis(int its, int *N, int **neighbors, double k_B, double T,
           double E_cucu, double E_znzn, double E_cuzn, double E_tot, gsl_rng *r, FILE *fp)
{   
    metro metro_result;
    idx index;
    double E;
    double E_prime;
    double delta_E;
    double alpha;
    int accepted = 0;
    for (int i = 0; i < its; i++)
    {   
        index = swappie(N, r);
        int A = index.alpha; // Index of the atom in sub_A
        int B = index.beta; // Index of the atom in sub_B
        int value_A = index.valueA;
        int value_B = index.valueB;
        E = energy_bond(index, N, neighbors, E_cucu, E_znzn, E_cuzn);
        N[A] = value_B;
        N[B] = value_A;
        E_prime = energy_bond(index, N, neighbors, E_cucu, E_znzn, E_cuzn);
        delta_E = E_prime - E;
        alpha = exp(-delta_E / k_B / T);
        if (gsl_rng_uniform(r) < alpha)
        {
            accepted++;
            E_tot += delta_E;
            if (fp != NULL)
            {
                fprintf(fp, "%i, %f\n", accepted, E_tot);
            }
        }
        else
        {    
            if (fp != NULL)
            {
                fprintf(fp, "%i, %f\n", accepted, E_tot);
            }
            N[A] = value_A;
            N[B] = value_B;
        }
    }
    metro_result.accepted = accepted;
    metro_result.Etot = E_tot;

    return metro_result;
}

idx
swappie(int *N, gsl_rng *r)
{   
    idx index;
    int idx_1 = (int)(1999 * gsl_rng_uniform(r));
    int A = N[idx_1];
    int idx_2 = (int)(1999 * gsl_rng_uniform(r));
    int B = N[idx_2];
    while (A == B)
    {
        idx_2 = (int)(1999 * gsl_rng_uniform(r));
        B = N[idx_2];
    }
    index.alpha = idx_1;
    index.beta = idx_2;
    index.valueA = A;
    index.valueB = B;

    return index;
}

double
energy_bond(idx index, int *N, int **neighbors, 
            double E_cucu, double E_znzn, double E_cuzn)
{   
    double E = 0.;
    int atom_1_idx = index.alpha;
    int atom_2_idx = index.beta;
    int atom_1 = N[atom_1_idx];
    int atom_2 = N[atom_2_idx];
    // int *neighbors_A = neigh_A[A];
    // int *neighbors_B = neigh_B[B];
    int *neighbors_1 = neighbors[atom_1_idx];
    int *neighbors_2 = neighbors[atom_2_idx];
    for (int i = 0; i < 8; i++)
    {
        if (N[neighbors_1[i]] == 1 && atom_1 == 1)
        {
            E += E_cucu;
        }
        else if (N[neighbors_1[i]] == 0 && atom_1 == 0)
        {
            E += E_znzn;
        }
        else
        {
            E += E_cuzn;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (N[neighbors_2[i]] == 1 && atom_2 == 1)
        {
            E += E_cucu;
        }
        else if (N[neighbors_2[i]] == 0 && atom_2 == 0)
        {
            E += E_znzn;
        }
        else
        {
            E += E_cuzn;
        }
    }

    return E;
}

void
lattice_to_files(FILE *fp_atoms, FILE *fp_neighbors, int *N, int **neighbors, int N_atoms)
{
    for (int i = 0; i < N_atoms; i++)
    {
        fprintf(fp_atoms, "%i\n", N[i]);
        for (int j = 0; j < 8; j++)
        {   
            if (j == 7)
            {
                fprintf(fp_neighbors, "%i\n", neighbors[i][j]);
            }
            else
            {
                fprintf(fp_neighbors, "%i, ", neighbors[i][j]);
            }
        }
    }
}
