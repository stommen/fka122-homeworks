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

typedef struct
{
    int N_Cu_A;
    int N_CuZn;
} atom_count;

gsl_rng * init_gsl_rng(int seed);
double boundary_distance_between_vectors(double *v1, double *v2, int dim, double box_length);
void nearest_neighbors_bcc(double **pos_A, double **pos_B, int N_atoms,
                           int **neighbors, double box_length,
                           double cutoff, double closest_distance);
metro metropolis(int its, int *atoms, int **neighbors, double k_B, double T, double E_cucu, 
                 double E_znzn, double E_cuzn, double E_tot, gsl_rng *r, 
                 double *C_V, double *U, double *P, double *R, int N_atoms, FILE *fp_eq);
idx swappy(int *atoms, gsl_rng *r);
double energy_bond(idx index, int *atoms, int **neighbors,
                   double E_cucu, double E_znzn, double E_cuzn);
void lattice_to_files(FILE *fp_atoms, FILE *fp_neighbors, int *atoms, int **neighbors, int N_atoms);
atom_count lattice_props(int *atoms, int **neighbors, int N_atoms);

int
run(int argc, char *argv[])
{
    // ------------------------------- Constants --------------------------------- //
    double E_cucu = -0.436; // eV
    double E_znzn = -0.113; // eV
    double E_cuzn = -0.294; // eV
    double k_B = 8.617333262145e-5; // eV/K
    double a = 2.949; // [Å]
    double closest_distance_bcc = sqrt(3) * a / 2;
    int L = 10;
    int N_atoms = 2 * L * L * L;
    double E_initial = 2000*E_cuzn;
    char filename[100];
    gsl_rng *r = init_gsl_rng(19);

    // ------------------------------- Task 1 --------------------------------- //

    // double step_size = 0.0001;
    // int N = (int)(1 / step_size) + 100;
    // FILE *file = fopen("data/task_1/data.csv", "w");
    // double delta_E = E_cucu + E_znzn - 2 * E_cuzn;
    // double *U = (double *)malloc(N * sizeof(double));
    // double *T = (double *)malloc(N * sizeof(double));
    // double *P = (double *)malloc((N + 1) * sizeof(double));
    // P[0] = 1.;
    // for (int i = 1; i < N; i++)
    // {   
    //     if (i < (int)(1 / step_size))
    //     {
    //         P[i] = P[i-1] - step_size;
    //         T[i-1] = 4 * P[i] * delta_E / k_B / log((1 + P[i]) / (1 - P[i]));
    //         U[i-1] = 2 * N_atoms * (E_cucu + E_znzn + 2 * E_cuzn) - 2 * N_atoms * P[i] * P[i] * delta_E;
    //         double C_V = variance(U, i) / k_B / T[i-1] / T[i-1];
    //         fprintf(file, "%f, %f, %f, %f\n", P[i], T[i-1], U[i-1], C_V);
    //     }
    //     else if (i >= (int)(1 / step_size))
    //     {
    //         P[i] = 0.;
    //         T[i-1] = T[i-2] + 5;
    //         U[i-1] = 2 * N_atoms * (E_cucu + E_znzn + 2 * E_cuzn) - 2 * N_atoms * P[i] * P[i] * delta_E;
    //         double C_V = variance(U, i) / k_B / T[i-1] / T[i-1];
    //         fprintf(file, "%f, %f, %f, %f\n", P[i], T[i-1], U[i-1], C_V);
    //     }
    // }

    // fclose(file);
    // free(T);
    // free(P);

    // ------------------------------- Task 2 --------------------------------- //

    // int its_eq = 100000;
    // double T = 1000;
    // metro metro_result;

    // double **sub_A = create_2D_array(N_atoms / 2, 3);
    // double **sub_B = create_2D_array(N_atoms / 2, 3);
    // int **neighbors = (int **)create_2D_array(N_atoms, 8);
    // int *atoms = (int *)calloc(N_atoms, sizeof(int));
    // for (int i = 0; i < N_atoms / 2; i++)
    // {
    //     atoms[i] = 1;
    // }

    // init_sc(sub_A, L, a, (double[3]){0, 0, 0});
    // init_sc(sub_B, L, a, (double[3]){0.5, 0.5, 0.5});
    // nearest_neighbors_bcc(sub_A, sub_B, N_atoms, neighbors, 10 * a, 0.001, closest_distance_bcc);

    // sprintf(filename, "data/task_2/equilibrium_%i_%.0f.csv", its_eq, T);
    // FILE *fp_eq = fopen(filename, "w");
    // fprintf(fp_eq, "accepted, E_tot\n");

    // metro_result = metropolis(its_eq, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_initial, r, NULL, NULL, NULL, NULL, N_atoms, fp_eq);
    // E_tot = metro_result.Etot;
    // int accepted = metro_result.accepted;
    // printf("Acceptance rate from equilibrium: %f\n", (double)accepted / its_eq);

    // int its = 1000000;
    // sprintf(filename, "data/task_2/energy_%i_%i_%.0f.csv", its_eq, its, T);
    // FILE *fp = fopen(filename, "w");
    // fprintf(fp, "Equilibrium iterations, Accepted ratio equilibrium, E_tot equilibrium\n");
    // fprintf(fp, "%i, %f, %f\n", its_eq, (double)accepted / its_eq, E_tot);
    // fprintf(fp, "accepted, E_tot\n");

    // sprintf(filename, "data/task_2/lattice/atoms_%i_%i_%.0f.csv", its_eq, its, T);
    // FILE *fp_atoms = fopen(filename, "w");
    // sprintf(filename, "data/task_2/lattice/neighbors.csv");
    // FILE *fp_neighbors = fopen(filename, "w");

    // metro_result = metropolis(its, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_tot, r, fp);

    // lattice_to_files(fp_atoms, fp_neighbors, atoms, neighbors, N_atoms);

    // destroy_2D_array(sub_A);
    // destroy_2D_array(sub_B);
    // destroy_2D_array((double **)neighbors);
    // free(atoms);
    // fclose(fp);
    // fclose(fp_eq);
    // fclose(fp_atoms);
    // gsl_rng_free(r);

    // ------------------------------- Task 3 --------------------------------- //

    double T_start = 300; // [K]
    double T_end = 1000; // [K]
    metro metro_result;
    metro metro_result_eq;
    atom_count lat_props;
    
    // sprintf(filename, "data/task_3/U.csv");
    // FILE *fp_U = fopen(filename, "w");
    // fprintf(fp_U, "T, U, Accepted ratio, Accepted ratio equilibrium\n");

    // sprintf(filename, "data/task_3/P_r.csv");
    // FILE *fp_P_r = fopen(filename, "w");
    // fprintf(fp_P_r, "T, N_Cu_A, N_CuZn, Accepted ratio, Accepted ratio equilibrium\n");

    sprintf(filename, "data/task_3/data.csv");
    FILE *fp_data = fopen(filename, "w");
    fprintf(fp_data, "T, U, C_V, P, r, Accepted ratio, Accepted ratio equilibrium\n");

    int its_eq;
    int dt = 25;
    for (double T = T_start; T < T_end+1; T+=dt)
    {   
        if (T < 600)
        {
            its_eq = 300000;
        }
        else
        {
            its_eq = 100000;
        }
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

        metro_result_eq = metropolis(its_eq, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_initial, r, NULL, NULL, NULL, NULL, N_atoms, NULL);

        int its = 1000000;
        // if ((int)T % 100)
        // {   
        //     if (T < 601)
        //     {
        //         its = 2000000;
        //     }
        //     else
        //     {
        //         its = 3000000;
        //     }
        // }
        double *U = (double *)calloc(its, sizeof(double));
        double *C_V = (double *)calloc(its, sizeof(double));
        double *P = (double *)calloc(its, sizeof(double));
        double *R = (double *)calloc(its, sizeof(double));
        metro_result = metropolis(its, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, metro_result_eq.Etot, r, U, C_V, P, R, N_atoms, NULL);

        lat_props = lattice_props(atoms, neighbors, N_atoms);
        // fprintf(fp_P_r, "%f, %i, %i, %f, %f\n", T, lat_props.N_Cu_A, lat_props.N_CuZn, (double)metro_result.accepted / its, (double)metro_result_eq.accepted / its_eq);

        double U_avg= average(U, metro_result.accepted);
        double C_V_inst = variance(U, metro_result.accepted) / k_B / T / T;
        double p = (2. * lat_props.N_Cu_A / (N_atoms / 2) - 1);
        double Rr = (lat_props.N_CuZn - 4. * N_atoms / 2) / (4 * N_atoms / 2);
        // fprintf(fp_U, "%f, %f, %f, %f, %f\n", T, U_avg, C_V_inst, (double)metro_result.accepted / its, (double)metro_result_eq.accepted / its_eq);
        fprintf(fp_data, "%f, %f, %f, %f, %f, %f, %f\n", T, U_avg, C_V_inst, p, Rr, (double)metro_result.accepted / its, (double)metro_result_eq.accepted / its_eq);

    //     if ((int)T % 100 == 0 && T < 1001)
    //     {   
    //         sprintf(filename, "data/task_3/test.csv");
    //         FILE *fp_test = fopen(filename, "w");
    //         fprintf(fp_test, "%i\n", metro_result.accepted);
    //         for (int i = 0; i < its; i++)
    //         {
    //             fprintf(fp_test, "%f, %f, %f, %f\n", U[i], C_V[i], P[i], R[i]);
    //         }

    //         sprintf(filename, "data/task_3/auto_corr_%i.csv", (int)T);
    //         FILE *fp_auto_corr = fopen(filename, "w");
    //         fprintf(fp_auto_corr, "N, Var_U, Var_CV, Var_P, Var_R\n");
    //         fprintf(fp_auto_corr, "%i, %f, %f, %f, %f\n", metro_result.accepted, variance(U, metro_result.accepted), variance(C_V, metro_result.accepted), variance(P, metro_result.accepted), variance(R, metro_result.accepted));
    //         fprintf(fp_auto_corr, "Lag, U, C_V, P, R\n");

    //         sprintf(filename, "data/task_3/blocking_%i.csv", (int)T);
    //         FILE *fp_blocking = fopen(filename, "w");
    //         fprintf(fp_blocking, "Block_size, U, C_V, P, R\n");
    //         for (int i = 0; i < 100000; i+=100)
    //         {   
    //             addition_with_constant(U, U, -average(U, metro_result.accepted), metro_result.accepted);
    //             addition_with_constant(C_V, C_V, -average(C_V, metro_result.accepted), metro_result.accepted);
    //             addition_with_constant(P, P, -average(P, metro_result.accepted), metro_result.accepted);
    //             addition_with_constant(R, R, -average(R, metro_result.accepted), metro_result.accepted);
    //             fprintf(fp_auto_corr, "%i, %f, %f, %f, %f\n", i, autocorrelation(U, metro_result.accepted, i), autocorrelation(C_V, metro_result.accepted, i), autocorrelation(P, metro_result.accepted, i), autocorrelation(R, metro_result.accepted, i));
    //         }
    //         for (int b = 1; b < metro_result.accepted*0.5; b+=100)
    //         {
    //             fprintf(fp_blocking, "%i, %f, %f, %f, %f\n", b, block_average(U, metro_result.accepted, b), block_average(C_V, metro_result.accepted, b), block_average(P, metro_result.accepted, b), block_average(R, metro_result.accepted, b));
    //         }
    //     }
        
        free(U);
        free(C_V);
        free(P);
        free(R);
        destroy_2D_array(sub_A);
        destroy_2D_array(sub_B);
        destroy_2D_array((double **)neighbors);
        free(atoms);
    }

    // fclose(fp_U);
    // fclose(fp_P_r);
    fclose(fp_data);
    gsl_rng_free(r);

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
metropolis(int its, int *atoms, int **neighbors, double k_B, double T, double E_cucu, 
           double E_znzn, double E_cuzn, double E_tot, gsl_rng *r, double *U, 
           double *C_V, double *P, double *R, int N_atoms, FILE *fp_eq)
{   
    metro metro_result;
    // atom_count lat_props;
    idx index;
    double E;
    double E_prime;
    double delta_E;
    double alpha;
    int accepted = 0;
    for (int i = 0; i < its; i++)
    {   
        index = swappy(atoms, r);
        int A = index.alpha; // Index of the atom in sub_A
        int B = index.beta; // Index of the atom in sub_B
        int value_A = index.valueA;
        int value_B = index.valueB;
        E = energy_bond(index, atoms, neighbors, E_cucu, E_znzn, E_cuzn);
        atoms[A] = value_B;
        atoms[B] = value_A;
        E_prime = energy_bond(index, atoms, neighbors, E_cucu, E_znzn, E_cuzn);
        delta_E = E_prime - E;
        alpha = exp(-delta_E / k_B / T);
        if (gsl_rng_uniform(r) < alpha)
        {   
            accepted++;
            E_tot += delta_E;
            if (fp_eq != NULL)
            {
                fprintf(fp_eq, "%i, %f\n", accepted, E_tot);
            }
            if (U != NULL)
            {   
                U[accepted-1] = E_tot;
                // if ((int)T % 100 == 0 && T < 1001)
                // {
                //     double U_fluct = variance(U, accepted);
                //     C_V[accepted-1] = U_fluct / k_B / T / T;
                //     lat_props = lattice_props(atoms, neighbors, N_atoms);
                //     P[accepted-1] = (2. * lat_props.N_Cu_A / (N_atoms / 2) - 1);
                //     R[accepted-1] = (lat_props.N_CuZn - 4. * N_atoms / 2) / (4 * N_atoms / 2);
                // }
            }
        }
        else
        {    
            if (fp_eq != NULL)
            {
                fprintf(fp_eq, "%i, %f\n", accepted, E_tot);
            }
            atoms[A] = value_A;
            atoms[B] = value_B;
        }
    }
    metro_result.accepted = accepted;
    metro_result.Etot = E_tot;

    return metro_result;
}

idx
swappy(int *atoms, gsl_rng *r)
{   
    idx index;
    int idx_1 = (int)(1999 * gsl_rng_uniform(r));
    int A = atoms[idx_1];
    int idx_2 = (int)(1999 * gsl_rng_uniform(r));
    int B = atoms[idx_2];
    while (A == B)
    {
        idx_2 = (int)(1999 * gsl_rng_uniform(r));
        B = atoms[idx_2];
    }
    index.alpha = idx_1;
    index.beta = idx_2;
    index.valueA = A;
    index.valueB = B;

    return index;
}

double
energy_bond(idx index, int *atoms, int **neighbors, 
            double E_cucu, double E_znzn, double E_cuzn)
{   
    double E = 0.;
    int atom_1_idx = index.alpha;
    int atom_2_idx = index.beta;
    int atom_1 = atoms[atom_1_idx];
    int atom_2 = atoms[atom_2_idx];
    int *neighbors_1 = neighbors[atom_1_idx];
    int *neighbors_2 = neighbors[atom_2_idx];
    for (int i = 0; i < 8; i++)
    {
        if (atoms[neighbors_1[i]] == 1 && atom_1 == 1)
        {
            E += E_cucu;
        }
        else if (atoms[neighbors_1[i]] == 0 && atom_1 == 0)
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
        if (atoms[neighbors_2[i]] == 1 && atom_2 == 1)
        {
            E += E_cucu;
        }
        else if (atoms[neighbors_2[i]] == 0 && atom_2 == 0)
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
lattice_to_files(FILE *fp_atoms, FILE *fp_neighbors, int *atoms, int **neighbors, int N_atoms)
{
    for (int i = 0; i < N_atoms; i++)
    {
        fprintf(fp_atoms, "%i\n", atoms[i]);
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

atom_count
lattice_props(int *atoms, int **neighbors, int N_atoms)
{
    atom_count count;
    int N_Cu_A = 0;
    int N_CuZn = 0;
    for (int i = 0; i < N_atoms / 2; i++)
    {
        if (atoms[i] == 1)
        {
            N_Cu_A++;
        }
        int atom = atoms[i];
        for (int j = 0; j < 8; j++)
        {
            if (atom == 1 && atoms[neighbors[i][j]] == 0)
            {
                N_CuZn++;
            }
            else if (atom == 0 && atoms[neighbors[i][j]] == 1)
            {
                N_CuZn++;
            }
        }
    }
    count.N_Cu_A = N_Cu_A;
    count.N_CuZn = N_CuZn;

    return count;
}
