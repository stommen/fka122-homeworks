#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"

#include "tools.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void calculate(double *potential, double *virial, double **force,
           double **position, double cell_length, int nbr_atoms);
gsl_rng *init_gsl_rng(int seed);
void verlet(double **positions, double **velocities, double **forces, int its, 
            double a, double delta_t, double m, double k_B, int N, FILE *fp);
void equilibrium(int its, int N, double T_eq, double P_eq, double **positions, 
            double **velocities, double *T, double *P, double *a, double a_0, double m,
            double beta, double k_B, double delta_t, FILE *fp);

int
run(
    int argc,
    char *argv[]
   )
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <delta_t>\n", argv[0]);
        return 1;
    }
    // ******************************** Constants ********************************* //
    int N = 256; // Number of atoms
    double k_B = 8.617333262145 * 1e-5; // Boltzmann constant [eV/K]
    double m = 27.0 / 9649; // Aluminium mass [eV ps²/Å²]

    // ********************************** Task 1 ********************************** //

    // // Lattice parameters
    // double a0[] = {3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98,
    //                3.99, 4.0, 4.01, 4.02, 4.03, 4.04, 4.05, 4.06, 4.07};

    // const char *filename = "data/task_1/energies.csv";
    // FILE *fp = fopen(filename, "w");

    // for (int i = 0; i < 18; i++)
    // {   
    //     // Positions vector
    //     double **pos = create_2D_array(N, 3);
    //     // Initialize the fcc lattice
    //     init_fcc(pos, 4, a0[i]);
    //     // Calculate the energy
    //     double energy = get_energy_AL(pos, a0[i] * N, N);
    //     // Write the energy to the file
    //     fprintf(fp, "%f, %f\n", a0[i], energy);
    //     // Free the memory
    //     destroy_2D_array(pos);
    // }

    // fclose(fp);

    // ********************************* Task 2 ********************************** //

    int t_max = 5;
    const double delta_t = atof(argv[1]);
    int its = (int)(t_max / delta_t);

    double a_0 = 4.0;
    double **positions = create_2D_array(N, 3);
    double **forces = create_2D_array(N, 3);
    double **velocities = create_2D_array(N, 3);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            velocities[i][j] = 0.0;
        }
    }
    init_fcc(positions, 4, a_0);
    // Initialize the random number generator
    gsl_rng *r = init_gsl_rng(42);

    // char filename[50];
    // Format the filename with delta_t
    // sprintf(filename, "data/task_2/data_dt_%.3f_%i.csv", delta_t, its);
    // FILE *fp = fopen(filename, "w");
    // fprintf(fp, "%i, %i, %f, %f\n", its, t_max, delta_t, 0.);

    for (int i = 0; i < N; i++)
        {   
            double rng = 0.935 + 0.13 * gsl_rng_uniform(r);
            addition_with_constant(positions[i], positions[i], a_0 * (1 - rng), 3);
        }

    verlet(positions, velocities, forces, its, a_0, delta_t, m, k_B, N, NULL);

    // fclose(fp);
    gsl_rng_free(r);
    // destroy_2D_array(positions);
    destroy_2D_array(forces);
    // destroy_2D_array(velocities);

    // ********************************* Task 3 ********************************** //

    double beta = 1. / (76e9 * 6.2415103e-12); // In [Å³/eV]
    double T_eq = 500.;
    double P_eq = 0.1e6 * 6.2415103e-12; // In [eV/Å³]
    int its_new = 5000;
    // double t_max = 10;
    // double a_0 = 4.0;

    // double delta_t = 0.001;
    // int its = (int)(t_max / delta_t);
    FILE *fp_pos = fopen("data/task_3/positions.xyz", "w");
    write_xyz(fp_pos, "Al", positions, velocities, a_0 * N, N);

    // Format the filename with delta_t
    char filename[50];
    sprintf(filename, "data/task_3/data_dt_%.3f_%i.csv", delta_t, its);
    FILE *fp = fopen(filename, "w");

    // double **positions = create_2D_array(N, 3);
    // double **velocities = create_2D_array(N, 3);
    double *T = malloc(its_new * sizeof(double));
    double *P = malloc(its_new * sizeof(double));
    double *a = malloc(its_new * sizeof(double));
    // for (int i = 0; i < N; i++)
    // {
        // for (int j = 0; j < 3; j++)
        // {
        //     velocities[i][j] = 0.0;
        // }
    // }
    // init_fcc(positions, 4, a_0);
    equilibrium(its_new, N, T_eq, P_eq, positions, velocities, T, P, a, a_0, m, beta, k_B, delta_t, fp);

    fclose(fp);
    destroy_2D_array(positions);
    destroy_2D_array(velocities);
    free(T);
    free(P);

    // ********************************* Task 4 ********************************** //
    
    return 0;
}

gsl_rng *
init_gsl_rng(int seed){
    const gsl_rng_type * T;
    gsl_rng * r; // r is a pointer, so if returner
    gsl_rng_env_setup();
    T = gsl_rng_default; // default random number generator
    r = gsl_rng_alloc(T); // allocate memory for the random number generator

    if (!r) {
        fprintf(stderr, "Error: Could not allocate memory for RNG.\n");
        exit(EXIT_FAILURE); // Exit if allocation fails
    }

    gsl_rng_set(r, seed);

    return r;
}

void 
verlet(double **positions, double **velocities, double **forces, int its, 
            double a, double delta_t, double m, double k_B, int N, FILE *fp)
{   
    double E_pot = 0.0;
    calculate(&E_pot, NULL, forces, positions, a * N, N);
    double T_avg = 0.0;
    double *E_kin_tot = malloc(its * sizeof(double));
    for (unsigned int i = 0; i < its; i++)
    {   
        double E_kin = 0.0;
        // Perform the first half step
        for (unsigned int j = 0; j < N; j++)
        {
            double *vels = velocities[j];
            E_kin += 0.5 * m * vector_norm(vels, 3) * vector_norm(vels, 3);
            for (unsigned int k = 0; k < 3; k++)
            {
                vels[k] += 0.5 * delta_t * forces[j][k] / m;
                positions[j][k] += delta_t * vels[k];
            }
        }
        E_kin_tot[i] = E_kin;
        T_avg = average(E_kin_tot, i+1) * 2.0 / 3.0 / k_B / N;
        if (fp != NULL){
            fprintf(fp, "%f, %f, %f, %f\n", E_kin, E_pot, E_kin + E_pot, T_avg);
        }

        // Calculate new accelerations
        calculate(&E_pot, NULL, forces, positions, a * N, N);

        // Perform the second half step
        for (unsigned int i = 0; i < N; i++)
        {
            double *vels = velocities[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                vels[j] += 0.5 * delta_t * forces[i][j] / m;
            }
        }
    }
}

void
equilibrium(int its, int N, double T_eq, double P_eq, 
            double **positions, double **velocities, double *T, double *P, 
            double *a, double a_0, double m, double beta, double k_B, double delta_t, FILE *fp)
{
    double tau_T = delta_t * 100; double tau_P = delta_t * 500;
    double alpha_T;
    double alpha_P = 1.;
    a[0] = a_0;
    for (unsigned int i = 0; i < its; i++)
    {   
        double E_k = 0.;
        for (unsigned int j = 0; j < N; j++)
        {
            E_k += 0.5 * m * vector_norm(velocities[j], 3) * vector_norm(velocities[j], 3);
        }
        T[i] = 2.0 / 3.0 / k_B / N * E_k;
        P[i] = 1 / (3 * 64 * pow(a_0, 3) * alpha_P) * (E_k + get_virial_AL(positions, a[i] * N, N));

        alpha_T = 1 + 2 * delta_t * (T_eq - T[i]) / (T[i] * tau_T);
        alpha_P = 1 - beta * delta_t * (P_eq - P[i]) / tau_P;
        fprintf(fp, "%f, %f, %f, %f, %f\n", T[i], P[i], a[i], alpha_T, alpha_P);

        if (i < its - 1){
            a[i+1] = a[i] * pow(alpha_P, 1. / 3.);
        }

        for (unsigned int j = 0; j < N; j++)
        {
            multiplication_with_constant(velocities[j], velocities[j], sqrt(alpha_T), 3);
            multiplication_with_constant(positions[j], positions[j], pow(alpha_P, 1. / 3.), 3);
        }
    }
}
