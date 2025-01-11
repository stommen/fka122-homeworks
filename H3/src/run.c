#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"
#include "run.h"

int
run(
    int argc,
    char *argv[]
   )
{
    // -------------------------------- Constants -------------------------------- //
    gsl_rng *U = init_gsl_rng(19); // Random number generator
    double gamma = 0.5; // Energy damping factor
    DMC_results results; // Results struct
    char filename[100]; // Filename

    // --------------------------------- Task 1a --------------------------------- //
    int N_0 = 200; // Initial number of walkers
    double dtau = 0.02; // Time step
    double tau = 5000; // Total time
    int N_its = tau / dtau; // Number of iterations
    double E_T_start = 0.5; // Initial energy
    int its_eq = 20 / dtau; // Number of equilibration steps

    // Initialize the walkers
    double **walkmen_pos = create_2D_array(N_0 * 100, 1);
    init_walkmen_1D(walkmen_pos, N_0);

    // Equilibrate the walkers
    results = DMC(its_eq, N_0, N_0+1, dtau, E_T_start, U, gamma, NULL, NULL, 1, false, 1);
    E_T_start = results.E_T;
    int N_sprinters = results.N_sprinters;

    // Run the simulation
    FILE* fp_ET = fopen("data/task_1/1D/ET_Nwalk_non_eq.csv", "w");
    fprintf(fp_ET, "E_T_avg, N_sprinters, N_survived/N_sprinters (%%), E_T\n");
    FILE* fp_w = fopen("data/task_1/1D/I_was_walkin_in_morse.csv", "w");

    DMC(N_its, N_0, N_sprinters, dtau, E_T_start, U, gamma, fp_ET, fp_w, 1, false, 1);

    // Close the files
    fclose(fp_ET);
    fclose(fp_w);

    // --------------------------------- Task 1b --------------------------------- //
    int N_0 = 1000; // Initial number of walkers
    int N_sprinters_init = N_0; // Initial number of sprinters
    double dtau = 0.01; // Time step
    double tau = 5000; // Total time
    int N_its = tau / dtau; // Number of iterations
    double E_T_start_init = -3; // Initial energy
    int its_eq = 100 / dtau; // Number of equilibration steps

    // Initialize the walkers
    double **walkmen_pos_init = create_2D_array(N_0 * 10, 6); 
    init_walkmen_6D(walkmen_pos_init, N_0, U);

    // Equilibrate the walkers
    results = DMC(walkmen_pos_init, its_eq, N_0, N_sprinters_init, dtau, E_T_start_init, U, gamma, NULL, NULL, 6, false, 1);
    double E_T_start = results.E_T;
    int N_sprinters = results.N_sprinters;

    // Run the simulation
    sprintf(filename, "data/task_1/6D/ET_Nwalk_non_eq.csv");
    FILE* fp_ET = fopen(filename, "w");
    fprintf(fp_ET, "E_T_avg, N_sprinters, N_survived/N_sprinters (%%), E_T\n");
    FILE* fp_w = fopen("data/task_1/6D/I_was_walkin_in_morse.csv", "w");

    results = DMC(walkmen_pos_init, N_its, N_0, N_sprinters, dtau, E_T_start, U, gamma, fp_ET, NULL, 6, false, 1);

    // Close the files
    fclose(fp_ET);
    fclose(fp_w);

    // --------------------------------- Task 2a --------------------------------- //
    dtau = 0.1; // Time step
    its_eq = 1000 / dtau; // Number of equilibration steps
    int decomp = 2; // Decomposition
    N_its = 50000; // Number of iterations

    // Equilibrate the walkers
    results = DMC(its_eq, N_0, N_sprinters_init, dtau, E_T_start_init, U, gamma, NULL, NULL, 6, true, decomp);
    double E_T_start = results.E_T;
    int N_sprinters = results.N_sprinters;
    
    // Run the simulation
    sprintf(filename, "data/task_2/ET_Nwalk_non_eq_decomp_%i_dtau_%f.csv", decomp, dtau);
    FILE *fp_ET_IS = fopen(filename, "w");
    fprintf(fp_ET_IS, "E_T_avg, N_sprinters, N_survived/N_sprinters (%%), E_T\n");

    sprintf(filename, "data/task_2/I_was_walkin_in_morse_decomp_%i_dtau_%f.csv", decomp, dtau);
    FILE* fp_w_IS = fopen(filename, "w");
    results = DMC(N_its, N_0, N_sprinters, dtau, E_T_start, U, gamma, fp_ET_IS, fp_w_IS, 6, true, decomp);

    // Close the files
    fclose(fp_ET_IS);
    fclose(fp_w_IS);

    // --------------------------------- Task 2b --------------------------------- //
    double dtaus[10] = {0.01, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4}; // Time steps
    decomp = 1; // Decomposition

    // File to store the results
    sprintf(filename, "data/task_2/ET_Nwalks_%i.csv", decomp);
    FILE *fp_ET_IS = fopen(filename, "w");
    fprintf(fp_ET_IS, "dtau, E_T_avg\n");

    // Run the simulation
    for (int i = 0; i < 10; i++)
    {   
        dtau = dtaus[i]; // Time step
        printf("dtau: %f\n", dtau);

        // Equilibrate the walkers
        results = DMC(its_eq, N_0, N_sprinters_init, dtau, E_T_start_init, U, gamma, NULL, NULL, 6, true, decomp);
        double E_T_start = results.E_T;
        int N_sprinters = results.N_sprinters;

        // Run the simulation
        results = DMC(N_its, N_0, N_sprinters, dtau, E_T_start, U, gamma, NULL, NULL, 6, true, decomp);
        fprintf(fp_ET_IS, "%f, %f\n", dtau, results.E_T_avg);
        printf("E_T_avg: %f\n", results.E_T_avg);
    }

    // Close the file
    fclose(fp_ET_IS);

    return 0;
}

DMC_results
DMC(
    int N_its,
    int N_0,
    int N_sprinters,
    double dtau,
    double E_T_start,
    gsl_rng *U,
    double gamma,
    FILE *fp_ET,
    FILE *fp_w,
    int ndim,
    bool IS,
    int decomp
   )
{
    double E_T = E_T_start; // Initial energy
    double E_T_avg = E_T_start; // Initial average energy
    DMC_results results; // Results struct
    int N_survivers = N_sprinters; // Number of surviving walkers
    double **walkmen_pos = create_2D_array(N_sprinters, ndim); // Walkers positions
    // Initialize the walkers
    if (ndim == 1)
    {
        init_walkmen_1D(walkmen_pos, N_sprinters);
    }
    else
    {
        init_walkmen_6D(walkmen_pos, N_sprinters, U);
    }
    for (int i = 0; i < N_its ; i++)
    {   
        // Generate new positions
        update_positions(walkmen_pos, N_survivers, ndim, dtau, fp_w, U, IS, decomp);
        int *num_walkers = (int *)calloc(N_sprinters, sizeof(int)); // Array to store the offsprings of each walker
        // Calculate the offspring of each walker
        for(int j = 0; j < N_sprinters; j++)
        {   
            if (ndim == 1)
            {
                int m = (int)(weight(dtau, E_T, morse(walkmen_pos[j][0])) + gsl_rng_uniform(U) * 1.); // Number of offsprings
                num_walkers[j] = m; // Store the number of offsprings
            }
            else if (ndim == 6)
            {   
                if (IS)
                {
                    double *R = walkmen_pos[j]; // Walker's position
                    double E_L = local_energy(R, 0.15); // Local energy
                    int m = (int)(weight(dtau, E_T, E_L) + gsl_rng_uniform(U) * 1.); // Number of offsprings
                    num_walkers[j] = m; // Store the number of offsprings
                }
                else
                {   
                    // Convert to spherical coordinates
                    double r_1 = sqrt(pow(walkmen_pos[j][0], 2) + pow(walkmen_pos[j][1], 2) + pow(walkmen_pos[j][2], 2));
                    double r_2 = sqrt(pow(walkmen_pos[j][3], 2) + pow(walkmen_pos[j][4], 2) + pow(walkmen_pos[j][5], 2));
                    double r_12 = sqrt(pow(walkmen_pos[j][3] - walkmen_pos[j][0], 2) + pow(walkmen_pos[j][4] - walkmen_pos[j][1], 2) + pow(walkmen_pos[j][5] - walkmen_pos[j][2], 2));
                    double V_tot = - 2 / sqrt(pow(r_1, 2) + 1e-4) - 2 / sqrt(pow(r_2, 2) + 1e-4) + 1 / sqrt(pow(r_12, 2) + 1e-4); // Potential energy
                    int m = (int)(weight(dtau, E_T, V_tot) + gsl_rng_uniform(U) * 1.); // Number of offsprings
                    num_walkers[j] = m; // Store the number of offsprings
                }
            }
        }
        N_survivers = int_sum(num_walkers, N_sprinters); // Number of surviving walkers
        double **walkmen_pos_new = create_2D_array(N_survivers, ndim); // New array to store the surviving walkers

        int M = 0; // Index for the new array
        // Populate the new array with the walkers
        for (int k = 0; k < N_sprinters; k++) 
        {
            for (int m = 0; m < num_walkers[k]; m++) // If the walker survives (num_walkers[j] > 0)
            {   
                for(int n = 0; n < ndim; n++)
                {
                    walkmen_pos_new[M][n] = walkmen_pos[k][n]; // Copy walker's position
                }
                M++; // Increment the new array index
            }
        }
        destroy_2D_array(walkmen_pos); // Free the memory of the old array
        walkmen_pos = create_2D_array(N_survivers, ndim); // Create a new array to store the walkers
        // Update the walker's position
        for (int k = 0; k < N_survivers; k++)
        {
            for (int n = 0; n < ndim; n++)
            {
                walkmen_pos[k][n] = walkmen_pos_new[k][n];
            }
        }

        // Updating ET
        double sprinter_ratio = (double)N_survivers / (double)N_0; // Ratio of surviving walkers
        E_T = E_T_avg - gamma * log(sprinter_ratio); // E_T at i+1
        E_T_avg = E_T / (i+1) + E_T_avg * i / (i+1); // E_T_avg at i+1

        // Print the results to the file
        if (fp_ET != NULL)
        {
            fprintf(fp_ET, "%lf, %i, %f, %f\n", E_T_avg, N_survivers, (100. * (double)N_survivers) / (double)N_sprinters, E_T);
        }
        N_sprinters = N_survivers; // Update the number of sprinters

        // Free the memory
        destroy_2D_array(walkmen_pos_new);
        free(num_walkers);
    }
    // Store the results in the struct
    results.E_T = E_T;
    results.E_T_avg = E_T_avg;
    results.N_sprinters = N_sprinters;
    results.R = walkmen_pos;

    return results;
}

void
init_walkmen_1D(
             double **walkers, 
             int N_walkers
            )
{
    int j = 0; // Index
    // Initialize the walkers
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
    // Initialize the walkers
    for(int i = 0; i < N_walkers; i++)
    {
        double r_1 = 0.7 + gsl_rng_uniform(U);
        double theta_1 = acos(2 * gsl_rng_uniform(U) - 1);
        double phi_1 = 2 * M_PI * gsl_rng_uniform(U);
        double r_2 = 0.7 + gsl_rng_uniform(U);
        double theta_2 = acos(2 * gsl_rng_uniform(U) - 1);
        double phi_2 = 2 * M_PI * gsl_rng_uniform(U);

        walkers[i][0] = r_1 * sin(theta_1) * cos(phi_1); // x_1
        walkers[i][1] = r_1 * sin(theta_1) * sin(phi_1); // y_1
        walkers[i][2] = r_1 * cos(theta_1); // z_1
        walkers[i][3] = r_2 * sin(theta_2) * cos(phi_2); // x_2
        walkers[i][4] = r_2 * sin(theta_2) * sin(phi_2); // y_2
        walkers[i][5] = r_2 * cos(theta_2); // z_2
    }
}

double
morse(
      double x
     )
{
    double V = 0.5 * pow((1  - exp(-x)),2); // Morse potential

    return V;
}

double
weight(
       double dtau, 
       double E_T, 
       double V
      )
{   
    double W = exp((E_T - V) * dtau); // Weight of the walker

    return W;
}

double local_energy(
                    double *R,
                    double alpha
                   )
{
    double x_1 = R[0];
    double y_1 = R[1];
    double z_1 = R[2];
    double x_2 = R[3];
    double y_2 = R[4];
    double z_2 = R[5];

    double r_12 = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2)); // Distance between the electrons

    double r1_mag = sqrt(pow(x_1, 2) + pow(y_1, 2) + pow(z_1, 2)); // Magnitude of the first electron
    double r2_mag = sqrt(pow(x_2, 2) + pow(y_2, 2) + pow(z_2, 2)); // Magnitude of the second electron

    // Unit vector of the first electron
    double r1_hat_x = x_1 / r1_mag;
    double r1_hat_y = y_1 / r1_mag;
    double r1_hat_z = z_1 / r1_mag;

    // Unit vector of the second electron
    double r2_hat_x = x_2 / r2_mag;
    double r2_hat_y = y_2 / r2_mag;
    double r2_hat_z = z_2 / r2_mag;

    // Unit vector between the electrons
    double delta_hat_x = r2_hat_x - r1_hat_x;
    double delta_hat_y = r2_hat_y - r1_hat_y;
    double delta_hat_z = r2_hat_z - r1_hat_z;

    // Delta vector between the electrons
    double delta_x = x_2 - x_1;
    double delta_y = y_2 - y_1;
    double delta_z = z_2 - z_1;

    double energy_dot = delta_hat_x * delta_x + delta_hat_y * delta_y + delta_hat_z * delta_z;

    // Local energy
    double E_L = - 4 + energy_dot / (r_12 * pow(1 + alpha * r_12, 2)) - 1 / (r_12 * pow(1 + alpha * r_12, 3)) - 1 / (4 * pow(1 + alpha * r1_mag, 4)) + 1 / r_12;

    return E_L;
}

void
update_positions(
                 double **walkmen_pos,
                 int N_survivers,
                 int ndim,
                 double dtau,
                 FILE *fp_w,
                 gsl_rng *U,
                 bool IS,
                 int decomp
                )
{
    // Update the positions of the walkers
    for (int l = 0; l < N_survivers; l++)
    {
        for (int n = 0; n < ndim; n++)
        {   
            // Diffusive part
            walkmen_pos[l][n] = walkmen_pos[l][n] + gsl_ran_gaussian(U, 1.) * sqrt(dtau);
        }
        if (IS)
        {      
            // Drift velocity
            drift_velocity(walkmen_pos[l], 0.15, dtau, decomp, U);
        }
        // Write the walker's position to the file
        if (fp_w != NULL)
        {   
            for (int n = 0; n < ndim; n++)
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
}

void
drift_velocity(
                double *R,
                double alpha,
                double dtau,
                int decomp,
                gsl_rng *U
              )
{
    double x_1 = R[0];
    double y_1 = R[1];
    double z_1 = R[2];
    double x_2 = R[3];
    double y_2 = R[4];
    double z_2 = R[5];

    double r_1 = sqrt(pow(x_1, 2) + pow(y_1, 2) + pow(z_1, 2)); // Magnitude of the first electron
    double r_2 = sqrt(pow(x_2, 2) + pow(y_2, 2) + pow(z_2, 2)); // Magnitude of the second electron

    double r12_vec[3] = {x_2 - x_1, y_2 - y_1, z_2 - z_1}; // Vector between the electrons
    double r12 = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2)); // Distance between the electrons

    double r1_hat[3] = {x_1 / r_1, y_1 / r_1, z_1 / r_1}; // Unit vector of the first electron
    double r2_hat[3] = {x_2 / r_2, y_2 / r_2, z_2 / r_2}; // Unit vector of the second electron
    double r12_hat[3] = {r12_vec[0] / r12, r12_vec[1] / r12, r12_vec[2] / r12}; // Unit vector between the electrons

    // Compute drift velocity for the electrons
    double v_F1[3];
    double v_F2[3];
    for (int i = 0; i < 3; i++) 
    {
        v_F1[i] = - 2 * r1_hat[i] - (1.0 / (2. * pow(1 + alpha * r12, 2))) * r12_hat[i];
        v_F2[i] = - 2 * r2_hat[i] + (1.0 / (2. * pow(1 + alpha * r12, 2))) * r12_hat[i];
    }

    // Second decomposition
    if (decomp == 2)
    {   
        // First half step
        x_1 += v_F1[0] * dtau / 2;
        y_1 += v_F1[1] * dtau / 2;
        z_1 += v_F1[2] * dtau / 2;
        x_2 += v_F2[0] * dtau / 2;
        y_2 += v_F2[1] * dtau / 2;
        z_2 += v_F2[2] * dtau / 2;
        double r_1 = sqrt(pow(x_1, 2) + pow(y_1, 2) + pow(z_1, 2)); // Magnitude of the first electron
        double r_2 = sqrt(pow(x_2, 2) + pow(y_2, 2) + pow(z_2, 2)); // Magnitude of the second electron

        double r12_vec[3] = {x_2 - x_1, y_2 - y_1, z_2 - z_1}; // Vector between the electrons
        double r12 = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2)); // Distance between the electrons

        double r1_hat[3] = {x_1 / r_1, y_1 / r_1, z_1 / r_1}; // Unit vector of the first electron
        double r2_hat[3] = {x_2 / r_2, y_2 / r_2, z_2 / r_2}; // Unit vector of the second electron
        double r12_hat[3] = {r12_vec[0] / r12, r12_vec[1] / r12, r12_vec[2] / r12}; // Unit vector between the electrons

        // Compute drift velocity for the electrons
        for (int i = 0; i < 3; i++)
        {
            v_F1[i] = -2 * r1_hat[i] - (1.0 / (2. * pow(1 + alpha * r12, 2))) * r12_hat[i];
            v_F2[i] = -2 * r2_hat[i] + (1.0 / (2. * pow(1 + alpha * r12, 2))) * r12_hat[i];
        }
    }

    // Update the positions
    R[0] += v_F1[0] * dtau;
    R[1] += v_F1[1] * dtau;
    R[2] += v_F1[2] * dtau;
    R[3] += v_F2[0] * dtau;
    R[4] += v_F2[1] * dtau;
    R[5] += v_F2[2] * dtau;
}

gsl_rng *
init_gsl_rng(
             int seed
            )
{
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default; // Default random number generator
    r = gsl_rng_alloc(T); // Allocate memory for the random number generator

    if (!r) {
        fprintf(stderr, "Error: Could not allocate memory for RNG.\n");
        exit(EXIT_FAILURE); // Exit if allocation fails
    }

    // Set the seed
    gsl_rng_set(r, seed);

    return r;
}
