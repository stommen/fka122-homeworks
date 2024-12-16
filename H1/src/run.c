#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include "run.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int
run(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <delta_t>\n", argv[0]);
        return 1;
    }
    // ____________________________________________________________________________ //
    // ******************************** Constants ********************************* //
    // ____________________________________________________________________________ //
    int N = 256; // Number of atoms
    int cell_length = 4; // Number of unit cells in each direction
    double k_B = 8.617333262145 * 1e-5; // Boltzmann constant [eV/K]
    double m = 27.0 / 9649; // Aluminium mass [eV ps²/Å²]
    double beta = 1. / (76e3); // Bulk Modulus inverse [1 / MPa]

    // ____________________________________________________________________________ //
    // ********************************** Task 1 ********************************** //
    // ____________________________________________________________________________ //
    // Lattice parameters
    double a0[] = {4.0, 4.005, 4.01, 4.015, 4.02, 4.025, 4.03, 4.035, 4.04, 4.045, 
                    4.05, 4.055, 4.06, 4.065, 4.07, 4.075, 4.08};

    // Energy file
    char *filename = "data/task_1/energies.csv";
    FILE *fp = fopen(filename, "w");

    // Iterate over the lattice parameters and calculate the energy
    for (int i = 0; i < 17; i++)
    {   
        // Positions vector
        double **pos = create_2D_array(N, 3);
        // Initialize the fcc lattice
        init_fcc(pos, 4, a0[i]);
        // Calculate the energy
        double energy = get_energy_AL(pos, a0[i] * cell_length, N);
        // Write the energy to the file
        fprintf(fp, "%f, %f\n", a0[i], energy);
        // Free the memory
        destroy_2D_array(pos);
    }

    fclose(fp);

    // ____________________________________________________________________________ //
    // ********************** Initialize Task 2 & 3 & 4 ************************** //
    // ____________________________________________________________________________ //
    int t_max = 50; // The simulation time [ps]
    const double delta_t = atof(argv[1]); // The time step [ps]
    int its = (int)(t_max / delta_t); // Number of iterations

    double **positions = create_2D_array(N, 3); // Atom positions
    double **forces = create_2D_array(N, 3); // Forces
    double **velocities = create_2D_array(N, 3); // Velocities
    // Initialize the velocities to 0
    for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                velocities[i][j] = 0.0;
            }
        }

    // Initialize the random number generator
    gsl_rng *r = init_gsl_rng(42);

    double a_0 = 4.03; // Lattice constant at 0 K [Å]
    init_fcc(positions, 4, a_0); // Initialize lattice
    for (int i = 0; i < N; i++)
        {   
            double rng = 0.935 + 0.13 * gsl_rng_uniform(r); // Disturbance
            // Add the disturbance to all initial atom positions
            addition_with_constant(positions[i], positions[i], a_0 * (1 - rng), 3);
        }

    // ____________________________________________________________________________ //
    // ********************************* Task 2 *********************************** //
    // ____________________________________________________________________________ //
    int its_eq = 0; // Number of iterations for equilibration, 0 for task 2

    char filename[50];
    // Format the filename with delta_t, its and its_eq
    sprintf(filename, "data/task_2/data_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp = fopen(filename, "w");
    // Write the header
    fprintf(fp, "its, t_max, delta_t, its_eq, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0);
    fprintf(fp, "E_kin [eV], E_pot [eV], E_tot [eV], <T> [K], T [K], <P> [MPa], P [MPa], a [Å]\n");

    // Perform the Verlet algorithm
    a_0 = verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, cell_length, 0., 0., fp, NULL, NULL, NULL);
    
    // ____________________________________________________________________________ //
    // ********************************* Task 3 *********************************** //
    // ____________________________________________________________________________ //
    double T_eq = 500. + 273.15; // Temperature at equilibrium [K]
    double P_eq = 0.1; // Pressure at equilibrium [MPa]
    int its_eq = 25000; // Number of iterations for equilibration

    char filename[50];
    // Format the data filename with delta_t, its and its_eq
    sprintf(filename, "data/task_3/data_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp = fopen(filename, "w");
    // Write the data header
    fprintf(fp, "its, t_max, delta_t, its_eq, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0);
    fprintf(fp, "E_kin [eV], E_pot [eV], E_tot [eV], <T> [K], T [K], <P> [MPa], P [MPa], a [Å]\n");

    // Format the trajectory filename with delta_t, its and its_eq
    sprintf(filename, "data/task_3/trajs_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp_traj = fopen(filename, "w");
    // Write the trajectory header
    fprintf(fp_traj, "its, t_max, delta_t, its_eq, -, -, -, -, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0, 0, 0, 0, 0);
    fprintf(fp_traj, "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4\n");

    // Perform the Verlet algorithm
    a_0 = verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, cell_length, T_eq, P_eq, fp, fp_traj, NULL, NULL);

    // ____________________________________________________________________________ //
    // ********************************* Task 4 *********************************** //
    // ____________________________________________________________________________ //
    int its_eq = 50000; // Number of iterations for first equilibration phase, task 4
    double T_eq = 900. + 273.15; // Temperature at first equilibrium [K]
    double P_eq = 0.1; // Pressure at first equilibrium [MPa]

    char filename[50];
    // Format the data filename with delta_t, its and its_eq
    sprintf(filename, "data/task_4/data_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp_1 = fopen(filename, "w");
    // Write the data header
    fprintf(fp_1, "its, t_max, delta_t, its_eq, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0);
    fprintf(fp_1, "E_kin [eV], E_pot [eV], E_tot [eV], <T> [K], T [K], <P> [MPa], P [MPa], a [Å]\n");

    // Format the trajectory filename with delta_t, its and its_eq
    sprintf(filename, "data/task_4/trajs_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp_2 = fopen(filename, "w");
    // Write the trajectory header
    fprintf(fp_2, "its, t_max, delta_t, its_eq, -, -, -, -, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0, 0, 0, 0, 0);
    fprintf(fp_2, "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4\n");

    // Perform the Verlet algorithm for the first equilibration phase
    a_0 = verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, cell_length, T_eq, P_eq, fp_1, fp_2, NULL, NULL);

    t_max = 50; // The simulation time for the second phase [ps]
    its = (int)(t_max / delta_t); // Number of iterations for the second phase
    its_eq = 30000; // Number of iterations for the second equilibration phase
    T_eq = 700. + 273.15; // Temperature at the second equilibrium [K]

    // Format the data filename with delta_t, its and its_eq
    sprintf(filename, "data/task_4/data_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp_3 = fopen(filename, "w");
    // Write the data header
    fprintf(fp_3, "its, t_max, delta_t, its_eq, -, -, -, -\n%i, %i, %f, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, 0, 0, 0, 0);
    fprintf(fp_3, "E_kin [eV], E_pot [eV], E_tot [eV], <T> [K], T [K], <P> [MPa], P [MPa], a [Å]\n");

    // Format the trajectory filename with delta_t, its and its_eq
    sprintf(filename, "data/task_4/trajs_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp_4 = fopen(filename, "w");
    // Write the trajectory header
    fprintf(fp_4, "its, t_max [ps], delta_t [ps], its_eq, T_eq [K], P_eq [MPa], -, -, -, -, -, -\n%i, %i, %f, %i, %f, %f, %i, %i, %i, %i, %i, %i\n", its, t_max, delta_t, its_eq, T_eq, P_eq, 0, 0, 0, 0, 0, 0);
    fprintf(fp_4, "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4\n");

    // Format the radial distribution filename with delta_t, its and its_eq
    sprintf(filename, "data/task_4/rdist_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp_rdist = fopen(filename, "w");
    // Write the radial distribution header
    fprintf(fp_rdist, "its, its_eq, delta_t [ps], a_0 [Å]\n");

    // Format the structure factor filename with delta_t, its and its_eq
    sprintf(filename, "data/task_4/sfact_%.3f_%i_%i.csv", delta_t, its, its_eq);
    FILE *fp_sfact = fopen(filename, "a");
    // Write the structure factor header
    fprintf(fp_sfact, "its, its_eq, delta_t [ps], a_0 [Å]\n");
    
    // Perform the Verlet algorithm for the second phase
    a_0 = verlet(positions, velocities, forces, its, its_eq, delta_t, m, k_B, beta, a_0, N, cell_length, T_eq, P_eq, fp_3, fp_4, fp_rdist, fp_sfact);

    // Close files and free memory
    fclose(fp);
    fclose(fp_traj);
    fclose(fp_1);
    fclose(fp_2);
    fclose(fp_3);
    fclose(fp_4);
    fclose(fp_rdist);
    fclose(fp_sfact);
    gsl_rng_free(r);
    destroy_2D_array(positions);
    destroy_2D_array(forces);
    destroy_2D_array(velocities);
    
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
verlet(double **positions, double **velocities, double **forces, int its, int its_eq, 
            double delta_t, double m, double k_B, double beta, double a_0, int N, 
            int cell_length, double T_eq, double P_eq, FILE *fp, FILE *fp_traj, 
            FILE *fp_rdist, FILE *fp_sfact)
{   
    double tau_T = delta_t * 200; // Time constant for temperature [ps]
    double tau_P = delta_t * 1000; // Time constant for pressure [ps]
    double alpha_T; // Scaling factor for temperature
    double alpha_P; // Scaling factor for pressure
    double E_pot = 0.0; // Potential energy (instantaneous)
    double E_kin = 0.0; // Kinetic energy (instantaneous)
    double virial = 0.0; // Virial
    double *P = malloc(its * sizeof(double)); // Instantaneous pressure
    double *T = malloc(its * sizeof(double)); // Instantaneous temperature
    double T_avg = 0.0; // Average temperature
    double P_avg = 0.0; // Average pressure

    // Calculate the initial forces, potential energy and virial
    calculate(&E_pot, &virial, forces, positions, a_0 * cell_length, N);
    for (unsigned int i = 0; i < its; i++)
    {   
        E_kin = 0.0;
        // Perform the first half step
        for (unsigned int j = 0; j < N; j++)
        {
            for (unsigned int k = 0; k < 3; k++)
            {
                velocities[j][k] += 0.5 * delta_t * forces[j][k] / m;
                positions[j][k] += delta_t * velocities[j][k];
            }
        }
        // Calculate new accelerations
        calculate(&E_pot, &virial, forces, positions, a_0 * cell_length, N);
        // Perform the second half step
        for (unsigned int j = 0; j < N; j++)
    {
            for (unsigned int k = 0; k < 3; k++)
            {
                velocities[j][k] += 0.5 * delta_t * forces[j][k] / m;
            }
            E_kin += 0.5 * m * vector_norm(velocities[j], 3) * vector_norm(velocities[j], 3);
        }

        T[i] = 2.0 / 3.0 / k_B / N * E_kin; // Calculate the instantaneous temperature [K]
        T_avg = average(T, i+1); // Calculate the average temperature [K]
        P[i] = 1 / (3 * 64 * pow(a_0, 3)) * (E_kin + virial) / 6.2415 * 1e6; // Calculate the instantaneous pressure [MPa]
        P_avg = average(P, i+1); // Calculate the average pressure [MPa]

        // Scale the temperature and pressure if the equilibration phase is not over
        if (i < its_eq)
        {   
            // Calculate the scaling factors
            alpha_T = 1 + 2 * delta_t * (T_eq - T[i]) / (T[i] * tau_T);
            alpha_P = 1 - beta * delta_t * (P_eq - P[i]) / tau_P;

            // Scale the lattice constant
            a_0 = a_0 * pow(alpha_P, 1. / 3.);

            // Scale the temperature and velocities
            for (unsigned int j = 0; j < N; j++)
            {
                multiplication_with_constant(positions[j], positions[j], pow(alpha_P, 1. / 3.), 3);
                multiplication_with_constant(velocities[j], velocities[j], sqrt(alpha_T), 3);
            }
        }
        // Write the data (energy, temperature, pressure and lattice constant) to the data file
        if (fp != NULL)
        {
            fprintf(fp, "%f, %f, %f, %f, %f, %f, %f, %f\n", E_kin, E_pot, E_kin + E_pot, T_avg, T[i], P_avg, P[i], a_0);
        }
        // Write the trajectory data to the trajectory file
        if (fp_traj != NULL)
        {
            fprintf(fp_traj, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
                            positions[15][0], positions[15][1], positions[15][2], 
                            positions[27][0], positions[27][1], positions[27][2], 
                            positions[35][0], positions[35][1], positions[35][2],
                            positions[49][0], positions[49][1], positions[49][2]);
        }
        // Calculate the radial distribution and write it to the radial distribution file
        if (fp_rdist != NULL && i >= its_eq)
        {
            double bin_width = 0.05; // Width of the bins
            double L = a_0 * 4; // Length of the supercell
            int N_bins = (int)(L / 2 / bin_width); // Number of bins
            double *bins = calloc(N_bins, sizeof(double)); // Bins for the radial distribution
            // Write the header
            if (i == its_eq){
                fprintf(fp_rdist, "%i, %i, %f, %f\n", its, its_eq, delta_t, a_0);
            }

            // Calculate the radial distribution
            radial_dist(bins, positions, N_bins, bin_width, N, L);

            // Write the radial distribution to the file
            for (int j = 0; j < N_bins; j++)
            {   
                if (j == N_bins - 1)
                {
                    fprintf(fp_rdist, "%f", bins[j]);
                }
                else
                {
                    fprintf(fp_rdist, "%f, ", bins[j]);
                }
            }
            fprintf(fp_rdist, "\n");
        }
        // Calculate the structure factor and write it to the structure factor file
        if (fp_sfact != NULL && i >= its_eq)
        {
            int N_max = 10; // Maximum value in each direction for n_x, n_y and n_z
            double L = a_0 * 4; // Length of the supercell
            int n = 2*N_max + 1; // Number of grid points in each direction
            int n_points = n * n * n; // Total number of grid points
            double *S_q = (double *)malloc(n_points * sizeof(double)); // Structure factor
            int n_bins = 500; // Number of bins
            double q_max = 2 * M_PI / L * N_max * sqrt(3); // Maximum value of the magnitude of q
            double bin_width = q_max / n_bins; // Width of the bins
            
            // Initialize the grid and keep it throughout the iterations
            static double ****grid = NULL;
            if (i == its_eq)
            {
                grid = init_grid(N_max, L);
            }
            // Calculate the structure factor
            structure_factor(grid, positions, S_q, N, n);
            // Calculate the spherical average of the structure factor
            spherical_avg(grid, S_q, n, n_bins, bin_width, fp_sfact);
            // Free memory
            if (i == its - 1)
            {
                destroy_grid(grid, N_max);
                free(S_q);
            }
        }
    }
    // Free memory
    free(T);
    free(P);

    return a_0;
}

void
radial_dist(double *bins, double **positions, int N_bins, double bin_width, int N, double L)
{
    double r; // Distance between atoms
    double V = L * L * L; // Volume of the supercell
    double norm_factor; // Normalization factor
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {   
            if (i != j)
            {   
                // Calculate the distance between the atoms
                r = boundary_distance_between_vectors(positions[i], positions[j], 3, L);
                for (int l = 0; l < N_bins; l++)
                {   
                    // Check if the distance is within the bin, and increment the 
                    // bin if true
                    if (l * bin_width < r && r < (l + 1) * bin_width)
                    {
                        bins[l] += 1;
                    }
                }
            }
        }
    }
    for (int l = 0; l < N_bins; l++)
    {   
        bins[l] /= (N - 1); // Average over the number of atoms
        norm_factor = (N - 1) * 4 * M_PI * (3 * pow((l + 1), 2) - 3 * (l + 1) + 1) * pow(bin_width, 3) / 3 / V;
        bins[l] /= norm_factor; // Apply the scale factor
    }
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

double ****
init_grid(int N_max, double L)
{
    int N = 2*N_max + 1; // Number of grid points in each direction
    double ****grid = (double ****)malloc(N * sizeof(double ***)); // Initialize the grid
    for (int i = 0; i < N; i++)
    {   
        // Allocate memory for the grid
        grid[i] = (double ***)malloc(N * sizeof(double **));
        for (int j = 0; j < N; j++)
        {
            grid[i][j] = (double **)malloc(N * sizeof(double *));
            for (int k = 0; k < N; k++)
            {
                grid[i][j][k] = (double *)malloc(3 * sizeof(double));
                // Calculate the grid points
                grid[i][j][k][0] = 2 * M_PI / L * (i - N_max);
                grid[i][j][k][1] = 2 * M_PI / L * (j - N_max);
                grid[i][j][k][2] = 2 * M_PI / L * (k - N_max);
            }
        }
    }
    return grid;
}

void
destroy_grid(double ****grid, int N_max)
{
    int N = 2*N_max + 1;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                free(grid[i][j][k]);
            }
            free(grid[i][j]);
        }
        free(grid[i]);
    }
    free(grid);
}

void
structure_factor(double ****grid, double **positions, double *Sq, int N, int n)
{   
    double S_q_cos; // S(q) cosine term
    double S_q_sin; // S(q) sine term
    int it = 0; // Iterator
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {   
                // The q vector components
                double qx = grid[i][j][k][0];
                double qy = grid[i][j][k][1];
                double qz = grid[i][j][k][2];

                S_q_cos = 0.0;
                S_q_sin = 0.0;
                // Calculate dot product for each atom
                for (int l = 0; l < N; l++)
                {   
                    double dot_product = (positions[l][0] * qx + 
                                          positions[l][1] * qy + 
                                          positions[l][2] * qz);
                    S_q_cos += cos(dot_product);
                    S_q_sin += sin(dot_product);
                }
                // Calculate the structure factor at the given grid point
                Sq[it] = 1.0 / N * (S_q_cos * S_q_cos + S_q_sin * S_q_sin);
                it++; // Increment the iterator
            }
        }
        
    }
}

void
spherical_avg(double ****grid, double *Sq, int n, int n_bins, double bin_width, FILE *fp_sfact)
{   
    double *S_avg = calloc(n_bins, sizeof(double)); // Spherical average
    int *counts = calloc(n_bins, sizeof(double)); // Counts for each bin
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {   
                // Calculate the magnitude of the q vector
                double q = sqrt(grid[i][j][k][0] * grid[i][j][k][0] + 
                                grid[i][j][k][1] * grid[i][j][k][1] + 
                                grid[i][j][k][2] * grid[i][j][k][2]);
                int bin = (int)(q / bin_width); // Bin index
                // Bin index cannot exceed the number of bins
                if (bin == n_bins){
                    bin = n_bins - 1;
                }
                // Calculate the spherical average and increment the count
                S_avg[bin] += Sq[i * n * n + j * n + k];
                counts[bin]++;
            }
        }
    }
    // Write the spherical average to the file if the bin is not empty
    for (int b = 0; b < n_bins; b++) {
        if (counts[b] > 0) {
            if (b == n_bins - 1){
                fprintf(fp_sfact, "%f\n", S_avg[b] / counts[b]);
            }
            else{
                fprintf(fp_sfact, "%f, ", S_avg[b] / counts[b]);
            }
        }
    }
    free(S_avg);
    free(counts);
}
