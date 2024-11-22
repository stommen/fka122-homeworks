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
void verlet(double **positions, double **velocities, double **forces, int its, int its_eq,
            double delta_t, double m, double k_B, double a_0, double beta, int N, 
            int cell_length, double T_eq, double P_eq, FILE *fp, FILE *fp_traj);

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
    int cell_length = 4; // Number of atoms in each direction
    double k_B = 8.617333262145 * 1e-5; // Boltzmann constant [eV/K]
    double m = 27.0 / 9649; // Aluminium mass [eV ps²/Å²]
    double beta = 1. / (76e3); // Bulk Modulus inverse [1 / MPa]

    // ********************************** Task 1 ********************************** //

    // // Lattice parameters
    // double a0[] = {4.0, 4.005, 4.01, 4.015, 4.02, 4.025, 4.03, 4.035, 4.04, 4.045, 
    //                 4.05, 4.055, 4.06, 4.065, 4.07, 4.075, 4.08};

    // const char *filename = "data/task_1/energies.csv";
    // FILE *fp = fopen(filename, "w");

    // for (int i = 0; i < 17; i++)
    // {   
    //     // Positions vector
    //     double **pos = create_2D_array(N, 3);
    //     // Initialize the fcc lattice
    //     init_fcc(pos, 4, a0[i]);
    //     // Calculate the energy
    //     double energy = get_energy_AL(pos, a0[i] * cell_length, N);
    //     // Write the energy to the file
    //     fprintf(fp, "%f, %f\n", a0[i], energy);
    //     // Free the memory
    //     destroy_2D_array(pos);
    // }

    // fclose(fp);

    // ********************** Initialize Task 2 & 3 & 4 ************************** //
    
    int t_max = 25; // [ps]
    const double delta_t = atof(argv[1]);
    int its = (int)(t_max / delta_t);

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

    // Initialize the random number generator
    gsl_rng *r = init_gsl_rng(42);

    double a_0 = 4.05;
    init_fcc(positions, 4, a_0);
    for (int i = 0; i < N; i++)
        {   
            double rng = 0.935 + 0.13 * gsl_rng_uniform(r);
            addition_with_constant(positions[i], positions[i], a_0 * (1 - rng), 3);
        }

    // ********************************* Task 2 ********************************** //

    // int its_eq = 0;

    // char filename[50];
    // // Format the filename with delta_t
    // sprintf(filename, "data/task_2/data_%.3f_%i_%i.csv", delta_t, its, its_eq);
    // FILE *fp = fopen(filename, "w");
    // fprintf(fp, "its, t_max, delta_t, its_eq, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0);
    // fprintf(fp, "E_kin [eV], E_pot [eV], E_tot [eV], <T> [K], T [K], <P> [MPa], P [MPa], a [Å]\n");

    // verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, cell_length, 0., 0., fp, NULL);

    // ********************************* Task 3 ********************************** //

    double T_eq = 500. + 273.15; // [K]
    double P_eq = 0.1; // [MPa]
    int its_eq = 25000;

    // Format the filename with delta_t
    char filename[50];
    sprintf(filename, "data/task_3/data_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "its, t_max, delta_t, its_eq, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0);
    fprintf(fp, "E_kin [eV], E_pot [eV], E_tot [eV], <T> [K], T [K], <P> [MPa], P [MPa], a [Å]\n");

    // sprintf(filename, "data/task_3/trajs_%.3f_%i_%i.csv", delta_t, its, its_eq);
    // FILE *fp_traj = fopen(filename, "w");
    // fprintf(fp_traj, "its, t_max, delta_t, its_eq, -, -, -, -, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0, 0, 0, 0, 0);
    // fprintf(fp_traj, "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4\n");
    verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, cell_length, T_eq, P_eq, fp, NULL);
    // printf("%f, %f, %f\n", positions[0][0], positions[0][1], positions[0][2]);

    its_eq = 0;
    t_max = 15;
    its = (int)(t_max / delta_t);

    a_0 = 4.048232;
    verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, cell_length, T_eq, P_eq, fp, NULL);

    // ********************************* Task 4 ********************************** //

    // int its_eq = 20000;
    // double T_eq = 900. + 273.15; // [K]
    // double P_eq = 0.1; // [MPa]

    // char filename[50];
    // sprintf(filename, "data/task_4/trajs_%.3f_%i_%i.csv", delta_t, its, its_eq);
    // FILE *fp_1 = fopen(filename, "w");
    // fprintf(fp_1, "its, t_max, delta_t, its_eq, -, -, -, -, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0, 0, 0, 0, 0);
    // fprintf(fp_1, "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4\n");

    // verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, T_eq, P_eq, NULL, fp_1);

    // t_max = 40; // [ps]
    // its = (int)(t_max / delta_t);
    // its_eq = 25000;
    // T_eq = 700. + 273.15; // [K]

    // sprintf(filename, "data/task_4/trajs_%.3f_%i_%i.csv", delta_t, its, its_eq);
    // FILE *fp_2 = fopen(filename, "w");
    // fprintf(fp_2, "its, t_max [ps], delta_t [ps], its_eq, T_eq [K], P_eq [MPa], -, -, -, -, -, -\n%i, %i, %f, %i, %f, %f, %i, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, T_eq, P_eq, 0, 0, 0, 0, 0, 0);
    // fprintf(fp_2, "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4\n");

    // verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, T_eq, P_eq, NULL, fp_2);

    // fclose(fp_1);
    fclose(fp);
    gsl_rng_free(r);
    destroy_2D_array(positions);
    destroy_2D_array(forces);
    destroy_2D_array(velocities);
    
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
verlet(double **positions, double **velocities, double **forces, int its, int its_eq, 
            double delta_t, double m, double k_B, double beta, double a_0, int N, 
            int cell_length, double T_eq, double P_eq, FILE *fp, FILE *fp_traj)
{   
    double tau_T = delta_t * 200;
    double tau_P = delta_t * 1000;
    double alpha_T;
    double alpha_P = 1.;
    double E_pot = 0.0;
    double virial = 0.0;
    double *P = malloc(its * sizeof(double));
    double *T = malloc(its * sizeof(double));
    double T_avg = 0.0;
    double P_avg = 0.0;
    calculate(&E_pot, &virial, forces, positions, a_0 * cell_length, N);

    for (unsigned int i = 0; i < its; i++)
    {   
        double E_kin = 0.0;
        // Perform the first half step
        for (unsigned int j = 0; j < N; j++)
        {
            E_kin += 0.5 * m * vector_norm(velocities[j], 3) * vector_norm(velocities[j], 3);
            for (unsigned int k = 0; k < 3; k++)
            {
                velocities[j][k] += 0.5 * delta_t * forces[j][k] / m;
                positions[j][k] += delta_t * velocities[j][k];
            }
        }

        T[i] = 2.0 / 3.0 / k_B / N * E_kin; // [K]
        T_avg = average(T, i+1);
        P[i] = 1 / (3 * 64 * pow(a_0, 3) * alpha_P) * (E_kin + virial) / 6.2415 * 1e6; // [MPa]
        P_avg = average(P, i+1);
        if (i < its_eq){
            a_0 = a_0 * pow(alpha_P, 1. / 3.);
            T[0] = 0.1;

            alpha_T = 1 + 2 * delta_t * (T_eq - T[i]) / (T[i] * tau_T);
            alpha_P = 1 - beta * delta_t * (P_eq - P[i]) / tau_P;
        }

        if (fp != NULL){
            fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f\n", E_kin, E_pot, E_kin + E_pot, T_avg, T[i], P_avg, P[i], a_0);
        }
        if (fp_traj != NULL){
            fprintf(fp_traj, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
                            positions[0][0], positions[0][1], positions[0][2],
                            positions[15][0], positions[15][1], positions[15][2], 
                            positions[27][0], positions[27][1], positions[27][2], 
                            positions[35][0], positions[35][1], positions[35][2]);
        }

        // Calculate new accelerations
        calculate(&E_pot, &virial, forces, positions, a_0 * cell_length, N);

        // Perform the second half step
        for (unsigned int i = 0; i < N; i++)
        {
            for (unsigned int j = 0; j < 3; j++)
            {
                velocities[i][j] += 0.5 * delta_t * forces[i][j] / m;
            }
            if (i < its_eq){
                multiplication_with_constant(positions[i], positions[i], pow(alpha_P, 1. / 3.), 3);
                multiplication_with_constant(velocities[i], velocities[i], sqrt(alpha_T), 3);
            }
        }
    }
    free(T);
    free(P);
}
