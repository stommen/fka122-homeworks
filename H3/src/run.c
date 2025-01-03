#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"

typedef struct
{
    double E_T;
    double E_T_avg;
    int N_sprinters;
    double **R;
} DMC_results;


DMC_results DMC(
                double **walkmen_pos,
                int N_its,
                int N_0,
                int N_sprinters,
                double dtau,
                double tau,
                double E_T_start,
                gsl_rng *U,
                double gamma,
                FILE *fp_ET,
                FILE *fp_w,
                int ndim,
                bool IS,
                int decomp
               );

double morse(
             double x
            );

void wavef_ground_state_analytic(
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

double local_energy(
                    double *R,
                    double alpha
                   );

void update_positions(
                      double **walkmen_pos,
                      double **walkmen_pos_new,
                      int N_survivers,
                      int ndim,
                      double dtau,
                      FILE *fp_w,
                      gsl_rng *U,
                      bool IS,
                      int decomp
                     );

void drift_velocity(
                    double *R,
                    double alpha,
                    double dtau,
                    int decomp
                   );

int
run(
    int argc,
    char *argv[]
   )
{
    // -------------------------------- Constants -------------------------------- //
    gsl_rng *U = init_gsl_rng(19);
    double gamma = 0.5;
    DMC_results results;
    char filename[100];

    // --------------------------------- Task 1a --------------------------------- //
    // int N_0 = 200;
    // double dtau = 0.02;
    // double tau = 5000;
    // int N_its = tau / dtau;
    // double E_T_start = 0.5;
    // int its_eq = 20 / dtau;

    // double **walkmen_pos = create_2D_array(N_0 * 100, 1);
    // init_walkmen_1D(walkmen_pos, N_0);

    // results = DMC(walkmen_pos, its_eq, N_0, dtau, tau, E_T_start, U, gamma, NULL, NULL, 1, false, 1);
    // E_T_start = results.E_T;

    // FILE* fp_ET = fopen("data/task_1/1D/ET_Nwalk_non_eq.csv", "w");
    // FILE* fp_w = fopen("data/task_1/1D/I_was_walkin_in_morse.csv", "w");

    // DMC(walkmen_pos, N_its, N_0, dtau, tau, E_T_start, U, gamma, fp_ET, fp_w, 1, false, 1);

    // free(walkmen_pos);
    // fclose(fp_ET);
    // fclose(fp_w);

    // --------------------------------- Task 1b --------------------------------- //
    int N_0 = 1000;
    int N_sprinters = N_0;
    double dtau = 0.01;
    double tau = 1000;
    int N_its = tau / dtau;
    double E_T_start = -3;
    int its_eq = 100 / dtau;

    double **walkmen_pos = create_2D_array(N_0 * 100, 6);
    init_walkmen_6D(walkmen_pos, N_0, U);

    results = DMC(walkmen_pos, its_eq, N_0, N_0, dtau, tau, E_T_start, U, gamma, NULL, NULL, 6, false, 1);
    E_T_start = results.E_T;
    N_sprinters = results.N_sprinters;

    // FILE* fp_ET = fopen("data/task_1/6D/ET_Nwalk_non_eq.csv", "w");
    // FILE* fp_w = fopen("data/task_1/6D/I_was_walkin_in_morse.csv", "w");

    results = DMC(walkmen_pos, N_its, N_0, N_sprinters, dtau, tau, E_T_start, U, gamma, NULL, NULL, 6, false, 1);

    // fclose(fp_ET);
    // fclose(fp_w);

    // --------------------------------- Task 2a --------------------------------- //
    E_T_start = results.E_T;
    N_sprinters = results.N_sprinters;
    double **walkmen_pos_init = results.R;
    dtau = 0.3;
    its_eq = 1000 / dtau;
    int decomp = 1;
    N_its = 50000;

    results = DMC(walkmen_pos, its_eq, N_0, N_sprinters, dtau, tau, E_T_start, U, gamma, NULL, NULL, 6, true, decomp);
    E_T_start = results.E_T;
    N_sprinters = results.N_sprinters;
    walkmen_pos = results.R;

    sprintf(filename, "data/task_2/ET_Nwalk_non_eq_decomp_%i_dtau_%f.csv", decomp, dtau);
    FILE *fp_ET_IS = fopen(filename, "w");
    // FILE* fp_w_IS = fopen("data/task_2/6D/I_was_walkin_in_morse.csv", "w");

    results = DMC(walkmen_pos, N_its, N_0, N_sprinters, dtau, tau, E_T_start, U, gamma, fp_ET_IS, NULL, 6, true, decomp);

    // --------------------------------- Task 2b --------------------------------- //
    // double dtaus[6] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4};
    // decomp = 1;

    // sprintf(filename, "data/task_2/ET_Nwalks_%i.csv", decomp);
    // FILE *fp_ET_IS = fopen(filename, "w");

    // for (int i = 0; i < 6; i++)
    // {   
    //     dtau = dtaus[i];
    //     printf("dtau: %f\n", dtau);
    //     results = DMC(walkmen_pos_init, its_eq, N_0, N_sprinters, dtau, tau, E_T_start, U, gamma, NULL, NULL, 6, true, decomp);
    //     E_T_start = results.E_T;
    //     N_sprinters = results.N_sprinters;
    //     walkmen_pos = results.R;

    //     results = DMC(walkmen_pos, N_its, N_0, N_sprinters, dtau, tau, E_T_start, U, gamma, NULL, NULL, 6, true, decomp);

    //     fprintf(fp_ET_IS, "%f, %f\n", dtau, results.E_T_avg);
    // }

    free(walkmen_pos_init);
    free(walkmen_pos);
    fclose(fp_ET_IS);
    // fclose(fp_w_IS);
    
    return 0;
}

DMC_results
DMC(
    double **walkmen_pos,
    int N_its,
    int N_0,
    int N_sprinters,
    double dtau,
    double tau,
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
    double E_T = E_T_start;
    double E_T_avg = E_T_start;
    DMC_results results;

    for(int i = 0; i < N_its ; i++)
    {   
        int *num_walkers = (int *)calloc(N_sprinters, sizeof(int));
        for(int j = 0; j < N_sprinters; j++)
        {   
            if (ndim == 1)
            {
                int m = (int)(weight(dtau, E_T, morse(walkmen_pos[j][0])) + gsl_rng_uniform(U) * 1.);
                num_walkers[j] = m;
            }
            else if (ndim == 6)
            {   
                if (IS)
                {
                    double *R = walkmen_pos[j];
                    double E_L = local_energy(R, 0.15);
                    int m = (int)(weight(dtau, E_T, E_L) + gsl_rng_uniform(U) * 1.);
                    num_walkers[j] = m;
                }
                else
                {
                    double r_1 = walkmen_pos[j][0];
                    double r_2 = walkmen_pos[j][3];
                    double r_12 = sqrt(pow(r_1, 2) + pow(r_2, 2) - 2 * r_1 * r_2 * (sin(walkmen_pos[j][1]) * sin(walkmen_pos[j][4]) * cos(walkmen_pos[j][2] - walkmen_pos[j][5]) + cos(walkmen_pos[j][1]) * cos(walkmen_pos[j][4])));
                    double V_tot = - 2 / sqrt(pow(r_1, 2) + 1e-1) - 2 / sqrt(pow(r_2, 2) + 1e-1) + 1 / sqrt(pow(r_12, 2) + 1e-1);
                    int m = (int)(weight(dtau, E_T, V_tot) + gsl_rng_uniform(U) * 1.);
                    num_walkers[j] = m;
                }
            }
        }
        // Number of surviving walkers
        int N_survivers = int_sum(num_walkers, N_sprinters);

        // Giving birth to new walkers 
        double **walkmen_pos_new = create_2D_array(N_survivers, ndim);

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

        // Generate new positions
        update_positions(walkmen_pos, walkmen_pos_new, N_survivers, ndim, dtau, fp_w, U, IS, decomp);

        // Updating ET
        double sprinter_ratio = (double)N_survivers / (double)N_0;
        E_T = E_T_avg - gamma * log(sprinter_ratio); // E_T at i+1
        E_T_avg = E_T / (i+1) + E_T_avg * i / (i+1); // E_T_avg at i+1

        if (fp_ET != NULL)
        {
            fprintf(fp_ET, "%lf, %i, %f, %f\n", E_T_avg, N_survivers, (100. * (double)N_survivers) / (double)N_sprinters, E_T);
        }
        
        N_sprinters = N_survivers;
        destroy_2D_array(walkmen_pos_new);
        free(num_walkers);
    }
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
        walkers[i][0] = 0.7 + gsl_rng_uniform(U); // r_1
        walkers[i][1] = acos(2 * gsl_rng_uniform(U) - 1); // theta_1
        walkers[i][2] = 2 * M_PI * gsl_rng_uniform(U); // phi_1
        walkers[i][3] = 0.7 + gsl_rng_uniform(U); // r_2
        walkers[i][4] = acos(2 * gsl_rng_uniform(U) - 1); // theta_2
        walkers[i][5] = 2 * M_PI * gsl_rng_uniform(U); // phi_2
    }
}

double
morse(
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
wavef_ground_state_analytic(
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

double local_energy(
                    double *R,
                    double alpha
                   )
{
    double r_1 = R[0];
    double theta_1 = R[1];
    double phi_1 = R[2];
    double r_2 = R[3];
    double theta_2 = R[4];
    double phi_2 = R[5];

    double x_1 = r_1 * sin(theta_1) * cos(phi_1);
    double y_1 = r_1 * sin(theta_1) * sin(phi_1);
    double z_1 = r_1 * cos(theta_1);

    double x_2 = r_2 * sin(theta_2) * cos(phi_2);
    double y_2 = r_2 * sin(theta_2) * sin(phi_2);
    double z_2 = r_2 * cos(theta_2);

    double r_12 = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));

    double r1_mag = sqrt(pow(x_1, 2) + pow(y_1, 2) + pow(z_1, 2));
    double r2_mag = sqrt(pow(x_2, 2) + pow(y_2, 2) + pow(z_2, 2));

    double r1_hat_x = x_1 / r1_mag;
    double r1_hat_y = y_1 / r1_mag;
    double r1_hat_z = z_1 / r1_mag;

    double r2_hat_x = x_2 / r2_mag;
    double r2_hat_y = y_2 / r2_mag;
    double r2_hat_z = z_2 / r2_mag;

    double delta_hat_x = r2_hat_x - r1_hat_x;
    double delta_hat_y = r2_hat_y - r1_hat_y;
    double delta_hat_z = r2_hat_z - r1_hat_z;

    double delta_x = x_2 - x_1;
    double delta_y = y_2 - y_1;
    double delta_z = z_2 - z_1;

    double energy_dot = delta_hat_x * delta_x + delta_hat_y * delta_y + delta_hat_z * delta_z;

    double E_L = - 4 + energy_dot / (r_12 * pow(1 + alpha * r_12, 2)) - 1 / (r_12 * pow(1 + alpha * r_12, 3)) - 1 / (4 * pow(1 + alpha * r1_mag, 4)) + 1 / r_12;

    return E_L;
}

void
update_positions(
                 double **walkmen_pos,
                 double **walkmen_pos_new,
                 int N_survivers,
                 int ndim,
                 double dtau,
                 FILE *fp_w,
                 gsl_rng *U,
                 bool IS,
                 int decomp
                )
{
    for (int l = 0; l < N_survivers; l++)
    {
        for (int n = 0; n < ndim; n++)
        {   
            if (IS)
            {   
                walkmen_pos[l][n] = walkmen_pos_new[l][n] + gsl_ran_gaussian(U, 1.) * sqrt(dtau);
            }
            else
            {
                walkmen_pos[l][n] = walkmen_pos_new[l][n] + gsl_ran_gaussian(U, 1.) * sqrt(dtau);
            }
            if (fp_w != NULL)
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
        drift_velocity(walkmen_pos[l], 0.15, dtau, decomp);
    }  
}

void
drift_velocity(
                double *R,
                double alpha,
                double dtau,
                int decomp
              )
{
    double r_1 = R[0];
    double theta_1 = R[1];
    double phi_1 = R[2];
    double r_2 = R[3];
    double theta_2 = R[4];
    double phi_2 = R[5];

    double x_1 = r_1 * sin(theta_1) * cos(phi_1);
    double y_1 = r_1 * sin(theta_1) * sin(phi_1);
    double z_1 = r_1 * cos(theta_1);

    double x_2 = r_2 * sin(theta_2) * cos(phi_2);
    double y_2 = r_2 * sin(theta_2) * sin(phi_2);
    double z_2 = r_2 * cos(theta_2);

    double r12_vec[3] = {x_2 - x_1, y_2 - y_1, z_2 - z_1};
    double r12 = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));

    double r1_hat[3] = {x_1 / r_1, y_1 / r_1, z_1 / r_1};
    double r2_hat[3] = {x_2 / r_2, y_2 / r_2, z_2 / r_2};
    double r12_hat[3] = {r12_vec[0] / r12, r12_vec[1] / r12, r12_vec[2] / r12};

    // Compute drift velocity for the first electron
    double v_F1[3];
    double v_F2[3];
    for (int i = 0; i < 3; i++) {
        v_F1[i] = -2 * r1_hat[i] - (1.0 / (2 * pow(1 + alpha * r12, 2))) * r12_hat[i];
        v_F2[i] = -2 * r2_hat[i] + (1.0 / (2 * pow(1 + alpha * r12, 2))) * r12_hat[i];
    }

    // Update positions with drift
    for (int i = 0; i < 3; i++)
    {
        if (decomp == 1)
        {   
            x_1 += v_F1[i] * dtau;
            y_1 += v_F1[i] * dtau;
            z_1 += v_F1[i] * dtau;
            x_2 += v_F2[i] * dtau;
            y_2 += v_F2[i] * dtau;
            z_2 += v_F2[i] * dtau;
        }
        else if (decomp == 2)
        {
            x_1 += v_F1[i] * dtau / 2;
            y_1 += v_F1[i] * dtau / 2;
            z_1 += v_F1[i] * dtau / 2;
            x_2 += v_F2[i] * dtau / 2;
            y_2 += v_F2[i] * dtau / 2;
            z_2 += v_F2[i] * dtau / 2; 
        }
    }

    if (decomp == 2)
    {   
        double r_1 = sqrt(pow(x_1, 2) + pow(y_1, 2) + pow(z_1, 2));
        double r_2 = sqrt(pow(x_2, 2) + pow(y_2, 2) + pow(z_2, 2));

        double r12_vec[3] = {x_2 - x_1, y_2 - y_1, z_2 - z_1};
        double r12 = sqrt(pow(x_2 - x_1, 2) + pow(y_2 - y_1, 2) + pow(z_2 - z_1, 2));

        double r1_hat[3] = {x_1 / r_1, y_1 / r_1, z_1 / r_1};
        double r2_hat[3] = {x_2 / r_2, y_2 / r_2, z_2 / r_2};
        double r12_hat[3] = {r12_vec[0] / r12, r12_vec[1] / r12, r12_vec[2] / r12};

        // Compute drift velocity for the first electron
        double v_F1[3];
        double v_F2[3];
        for (int i = 0; i < 3; i++) {
            v_F1[i] = -2 * r1_hat[i] - (1.0 / (2 * pow(1 + alpha * r12, 2))) * r12_hat[i];
            v_F2[i] = -2 * r2_hat[i] + (1.0 / (2 * pow(1 + alpha * r12, 2))) * r12_hat[i];
        }

        for (int i = 0; i < 3; i++)
        {
            x_1 += v_F1[i] * dtau;
            y_1 += v_F1[i] * dtau;
            z_1 += v_F1[i] * dtau;
            x_2 += v_F2[i] * dtau;
            y_2 += v_F2[i] * dtau;
            z_2 += v_F2[i] * dtau; 
        }
    }

    // Update the positions
    R[0] = sqrt(pow(x_1, 2) + pow(y_1, 2) + pow(z_1, 2));
    R[1] = acos(z_1 / R[0]);
    R[2] = atan2(y_1, x_1);
    R[3] = sqrt(pow(x_2, 2) + pow(y_2, 2) + pow(z_2, 2));
    R[4] = acos(z_2 / R[3]);
    R[5] = atan2(y_2, x_2);
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
