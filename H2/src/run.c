#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tools.h"
#include "lattice.h"
#include "run.h"

int
run(
    int argc, 
    char *argv[]
   )
{
    // ------------------------------- Constants --------------------------------- //

    double E_cucu = -0.436; // Bonding energy Cu-Cu [eV]
    double E_znzn = -0.113; // Bonding energy Zn-Zn [eV]
    double E_cuzn = -0.294; // Bonding energy Cu-Zn [eV]
    double k_B = 8.617333262145e-5; // Boltzmann constant [eV/K]
    double a = 2.949; // Lattice parameter [Å]
    double closest_distance_bcc = sqrt(3) * a / 2; // Closest distance between atoms in BCC [Å]
    int L = 10; // Number of unit cells in each direction
    int N_atoms = 2 * L * L * L; // Number of atoms
    double E_initial = N_atoms / 2 * E_cuzn * 8; // Initial energy
    char filename[100];
    gsl_rng *r = init_gsl_rng(19);

    // ------------------------------- Task 1 --------------------------------- //

    double step_size = 0.0001; // Step size for P
    int N = (int)(1 / step_size) + 100; // Number of points (add 100 for the end, where P = 0)
    FILE *file = fopen("data/task_1/data.csv", "w");
    double delta_E = E_cucu + E_znzn - 2 * E_cuzn; 
    double *Us = (double *)malloc(N * sizeof(double)); // Energy vector
    double *Ts = (double *)malloc(N * sizeof(double)); // Temperature vector
    double *Ps = (double *)malloc((N + 1) * sizeof(double)); // Pressure vector
    double C_V;
    Ps[0] = 1.; // Start at P = 1
    // Loop over P and calculate T, U and C_V
    for (int i = 1; i < N; i++)
    {   
        if (i < (int)(1 / step_size))
        {
            Ps[i] = Ps[i-1] - step_size; // Update P
            Ts[i-1] = 4 * Ps[i] * delta_E / k_B / log((1 + Ps[i]) / (1 - Ps[i])); // Calculate T
            Us[i-1] = N_atoms * (E_cucu + E_znzn + 2 * E_cuzn) - N_atoms * Ps[i] * Ps[i] * delta_E; // Calculate U
            if (i >= 2)
            {
                double dU = Us[i-1] - Us[i-2]; // dU
                double dT = Ts[i-1] - Ts[i-2]; // dT
                C_V = dU / dT; // Calculate C_V
            }
            else
            {
                C_V = 0.; // Set C_V to 0 for the first iteration
            }
            fprintf(file, "%f, %f, %f, %f\n", Ps[i], Ts[i-1], Us[i-1], C_V);
        }
        // When P = 0, set P = 0 and continue to calculate T, U and C_V
        else if (i >= (int)(1 / step_size))
        {
            Ps[i] = 0.;
            Ts[i-1] = Ts[i-2] + 5; // Use a step size of 5 K
            Us[i-1] = N_atoms * (E_cucu + E_znzn + 2 * E_cuzn) - N_atoms * Ps[i] * Ps[i] * delta_E;
            double C_V = (Us[i-1] - Us[i-2]) / (Ts[i-1] - Ts[i-2]);
            fprintf(file, "%f, %f, %f, %f\n", Ps[i], Ts[i-1], Us[i-1], C_V);
        }
    }

    // Free memory and close files
    fclose(file);
    free(Us);
    free(Ts);
    free(Ps);

    // ------------------------------- Task 2 --------------------------------- //

    int its_eq = 250000; // Number of iterations for equilibrium
    double T = 400; // Temperature [K]
    metro metro_result;

    double **sub_A = create_2D_array(N_atoms / 2, 3); // Sublattice A
    double **sub_B = create_2D_array(N_atoms / 2, 3); // Sublattice B
    int **neighbors = (int **)create_2D_array(N_atoms, 8); // Neighbors matrix
    int *atoms = (int *)calloc(N_atoms, sizeof(int)); // Atoms vector
    // Initialize all atoms in sublattice A to be Cu (1) (= Cold start)
    for (int i = 0; i < N_atoms / 2; i++)
    {
        atoms[i] = 1;
    }

    init_sc(sub_A, L, a, (double[3]){0, 0, 0}); // Initialize sublattice A
    init_sc(sub_B, L, a, (double[3]){0.5, 0.5, 0.5 }); // Initialize sublattice B
    nearest_neighbors_bcc(sub_A, sub_B, N_atoms, neighbors, 10 * a, 0.001, closest_distance_bcc); // Find nearest neighbors

    // File for equilibrium
    sprintf(filename, "data/task_2/equilibrium_%i_%.0f.csv", its_eq, T);
    FILE *fp_eq = fopen(filename, "w");
    fprintf(fp_eq, "accepted, E_tot\n");

    // Run metropolis for equilibrium
    metro_result = metropolis(its_eq, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_initial, r, NULL, NULL, NULL, NULL, N_atoms, fp_eq);
    double E_tot = metro_result.Etot;
    int accepted = metro_result.accepted;
    printf("Acceptance rate from equilibrium: %f\n", (double)accepted / its_eq);

    int its = 1000000; // Number of MC iterations after equilibrium
    // File for MC iterations
    sprintf(filename, "data/task_2/energy_%i_%i_%.0f.csv", its_eq, its, T);
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "Equilibrium iterations, Accepted ratio equilibrium, E_tot equilibrium\n");
    fprintf(fp, "%i, %f, %f\n", its_eq, (double)accepted / its_eq, E_tot);
    fprintf(fp, "accepted, E_tot\n");

    // File for atoms and neighbors
    sprintf(filename, "data/task_2/lattice/atoms_%i_%i_%.0f.csv", its_eq, its, T);
    FILE *fp_atoms = fopen(filename, "w");
    sprintf(filename, "data/task_2/lattice/neighbors.csv");
    FILE *fp_neighbors = fopen(filename, "w");

    // Run metropolis for MC iterations
    metro_result = metropolis(its, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_tot, r, NULL, NULL, NULL, NULL, N_atoms, fp);
    // Calculate lattice properties
    lattice_to_files(fp_atoms, fp_neighbors, atoms, neighbors, N_atoms);

    // Free memory and close files
    destroy_2D_array(sub_A);
    destroy_2D_array(sub_B);
    destroy_2D_array((double **)neighbors);
    free(atoms);
    fclose(fp);
    fclose(fp_eq);
    fclose(fp_atoms);
    gsl_rng_free(r);

    // ------------------------------- Task 3 --------------------------------- //

    double T_start = 300; // T value at MC simulation start [K]
    double T_end = 1000; // T value at MC simulation end [K]
    int dt = 25; // T step size [K]
    metro metro_result; 
    metro metro_result_eq;
    atom_count lat_props;
    
    // File for data
    sprintf(filename, "data/task_3/data.csv");
    FILE *fp_data = fopen(filename, "w");
    fprintf(fp_data, "T, U, C_V, P, r, Accepted ratio, Accepted ratio equilibrium, Iterations\n");

    int its_eq; 
    // Loop over T values
    for (double T = T_start; T < T_end+1; T+=dt)
    {   
        printf("T: %f\n", T); // Print T value, for debugging
        // Let the number of iterations for equilibrium depend on T
        if (T < 600)
        {
            its_eq = 250000;
        }
        else
        {
            its_eq = 100000;
        }
        double **sub_A = create_2D_array(N_atoms / 2, 3); // Sublattice A
        double **sub_B = create_2D_array(N_atoms / 2, 3); // Sublattice B
        int **neighbors = (int **)create_2D_array(N_atoms, 8); // Neighbors matrix
        int *atoms = (int *)calloc(N_atoms, sizeof(int)); // Atoms vector
        // Initialize all atoms in sublattice A to be Cu (1) (= Cold start)
        for (int i = 0; i < N_atoms / 2; i++)
        {
            atoms[i] = 1;
        }

        init_sc(sub_A, L, a, (double[3]){0, 0, 0}); // Initialize sublattice A
        init_sc(sub_B, L, a, (double[3]){0.5, 0.5, 0.5}); // Initialize sublattice B
        nearest_neighbors_bcc(sub_A, sub_B, N_atoms, neighbors, 10 * a, 0.001, closest_distance_bcc); // Find nearest neighbors

        // Run metropolis for equilibrium
        metro_result_eq = metropolis(its_eq, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, E_initial, r, NULL, NULL, NULL, NULL, N_atoms, NULL);

        int its = 100000; // Number of accepted MC iterations
        double *U = (double *)calloc(its, sizeof(double)); // Energy vector
        double *C_V = (double *)calloc(its, sizeof(double)); // Heat capacity vector
        double *P = (double *)calloc(its, sizeof(double)); // Pressure vector
        double *R = (double *)calloc(its, sizeof(double)); // Order parameter vector

        // Run metropolis for MC iterations
        metro_result = metropolis(its, atoms, neighbors, k_B, T, E_cucu, E_znzn, E_cuzn, metro_result_eq.Etot, r, U, C_V, P, R, N_atoms, NULL);
        lat_props = lattice_props(atoms, neighbors, N_atoms); // Calculate lattice properties, e.g., N_Cu_A and N_CuZn
        its = metro_result.its; // Number of iterations

        double U_avg = average(U, metro_result.accepted); // Average energy
        double C_V_inst = variance(U, metro_result.accepted) / k_B / T / T; // Instantaneous heat capacity
        double p = (2. * lat_props.N_Cu_A / (N_atoms / 2) - 1); // Long-range order parameter
        double Rr = (lat_props.N_CuZn - 4. * N_atoms / 2) / (4 * N_atoms / 2); // Short-range order parameter
        fprintf(fp_data, "%f, %f, %f, %f, %f, %f, %f, %i\n", T, U_avg, C_V_inst, p, Rr, (double)metro_result.accepted / its, (double)metro_result_eq.accepted / its_eq, its); 

        // Files for autocorrelation and blocking
        sprintf(filename, "data/task_3/auto_corr_%i.csv", (int)T);
        FILE *fp_auto_corr = fopen(filename, "w");
        fprintf(fp_auto_corr, "N, Var_U, Var_CV, Var_P, Var_R\n");
        fprintf(fp_auto_corr, "%i, %f, %f, %f, %f\n", metro_result.accepted, variance(U, metro_result.accepted), variance(C_V, metro_result.accepted), variance(P, metro_result.accepted), variance(R, metro_result.accepted));
        fprintf(fp_auto_corr, "Lag, U, C_V, P, R\n");

        sprintf(filename, "data/task_3/blocking_%i.csv", (int)T);
        FILE *fp_blocking = fopen(filename, "w");
        fprintf(fp_blocking, "Block_size, U, C_V, P, R\n");

        // Loop over block sizes
        for (int b = 1; b < metro_result.accepted*0.2; b+=10)
        {
            fprintf(fp_blocking, "%i, %f, %f, %f, %f\n", b, block_average(U, metro_result.accepted, b), block_average(C_V, metro_result.accepted, b), block_average(P, metro_result.accepted, b), block_average(R, metro_result.accepted, b));
        }
        // Loop over lags
        for (int i = 0; i < metro_result.accepted*0.5; i+=(metro_result.accepted*0.5/1000+1))
        {   
            // Subtract the average value from the data
            addition_with_constant(U, U, -average(U, metro_result.accepted), metro_result.accepted);
            addition_with_constant(C_V, C_V, -average(C_V, metro_result.accepted), metro_result.accepted);
            addition_with_constant(P, P, -average(P, metro_result.accepted), metro_result.accepted);
            addition_with_constant(R, R, -average(R, metro_result.accepted), metro_result.accepted);
            fprintf(fp_auto_corr, "%i, %f, %f, %f, %f\n", i, autocorrelation(U, metro_result.accepted, i), autocorrelation(C_V, metro_result.accepted, i), autocorrelation(P, metro_result.accepted, i), autocorrelation(R, metro_result.accepted, i));
        }
        
        // Free memory and close files
        free(U);
        free(C_V);
        free(P);
        free(R);
        destroy_2D_array(sub_A);
        destroy_2D_array(sub_B);
        destroy_2D_array((double **)neighbors);
        free(atoms);
    }

    // Free memory and close files
    fclose(fp_data);
    gsl_rng_free(r);

    return 0;
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

double
boundary_distance_between_vectors(
                                  double *v1, 
                                  double *v2, 
                                  int dim, 
                                  double box_length
                                 )
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
nearest_neighbors_bcc(
                      double **pos_A, 
                      double **pos_B, 
                      int N_atoms, 
                      int **neighbors, 
                      double box_length, 
                      double cutoff, 
                      double closest_distance
                     )
{
    double r; // Distance between atoms
    // Loop over all atoms in sublattice A
    for (int i = 0; i < N_atoms / 2; i++)
    {   
        int count = 0; // Counter for neighbors
        for (int j = 0; j < N_atoms / 2; j++)
        {
            r = boundary_distance_between_vectors(pos_A[i], pos_B[j], 3, box_length);
            // If the distance between atoms is within the cutoff distance, store the 
            // index of the atom in sublattice B (index > 1000 in atoms vector)
            // 
            if (fabs(closest_distance - r) < cutoff)
            {      
                neighbors[i][count] = j + 1000;
                count++;
            }
        }
    }
    // Loop over all atoms in sublattice B
    for (int i = 0; i < N_atoms / 2; i++)
    {   
        int count = 0; // Counter for neighbors
        for (int j = 0; j < N_atoms / 2; j++)
        {
            r = boundary_distance_between_vectors(pos_B[i], pos_A[j], 3, box_length);
            // If the distance between atoms is within the cutoff distance, store the
            // index of the atom in sublattice A (index < 1000 in atoms vector)
            if (fabs(closest_distance - r) < cutoff)
            {   
                neighbors[i+1000][count] = j;
                count++;
            }
        }
    }
}

metro
metropolis(
           int its, 
           int *atoms, 
           int **neighbors, 
           double k_B, 
           double T, 
           double E_cucu, 
           double E_znzn, 
           double E_cuzn, 
           double E_tot, 
           gsl_rng *r, 
           double *U, 
           double *C_V, 
           double *P, 
           double *R, 
           int N_atoms, 
           FILE *fp
          )
{   
    metro metro_result;
    atom_count lat_props;
    idx index;
    double E; // Energy before swap
    double E_prime; // Energy after swap
    double delta_E; // Energy difference
    double alpha; // Acceptance probability
    int accepted = 0; // Number of accepted swaps
    int i = 0; // Number of iterations
    // Loop over the number of iterations (task 2) or until its swaps are accepted (task 3)
    // for (int i = 0; i < its; i++)
    while (accepted < its)
    {   
        index = swappy(atoms, r); // Get the indices and values of the atoms to swap
        int A = index.alpha; // Index of the atom in sub_A
        int B = index.beta; // Index of the atom in sub_B
        int value_A = index.valueA; // Value of the atom in sub_A
        int value_B = index.valueB; // Value of the atom in sub_B
        E = energy_bond(index, atoms, neighbors, E_cucu, E_znzn, E_cuzn); // Energy before swap
        atoms[A] = value_B; // Swap the atoms
        atoms[B] = value_A; // Swap the atoms
        E_prime = energy_bond(index, atoms, neighbors, E_cucu, E_znzn, E_cuzn); // Energy after swap
        delta_E = E_prime - E; // Energy difference
        alpha = exp(-delta_E / k_B / T); // Acceptance probability
        // If the swap is accepted, increment the number of accepted swaps and update the energy
        if (gsl_rng_uniform(r) < alpha)
        {   
            accepted++;
            E_tot += delta_E;
            // Files for task 2
            if (fp != NULL)
            {
                fprintf(fp, "%i, %f\n", accepted, E_tot);
            }
            // For task 3, store the quantity values in the vectors (for autocorrelation and blocking)
            if (U != NULL)
            {   
                U[accepted-1] = E_tot;
                double U_fluct = variance(U, accepted);
                C_V[accepted-1] = U_fluct / k_B / T / T;
                lat_props = lattice_props(atoms, neighbors, N_atoms);
                P[accepted-1] = (2. * lat_props.N_Cu_A / (N_atoms / 2) - 1);
                R[accepted-1] = (lat_props.N_CuZn - 4. * N_atoms / 2) / (4 * N_atoms / 2);
            }
        }
        else
        {    
            // File for task 2
            if (fp != NULL)
            {
                fprintf(fp, "%i, %f\n", accepted, E_tot);
            }
            // If the swap is not accepted, swap the atoms back
            atoms[A] = value_A;
            atoms[B] = value_B;
        }
        i++; // Increment the number of iterations
    }
    // Store the values in the struct and return it
    metro_result.accepted = accepted;
    metro_result.Etot = E_tot;
    metro_result.its = i;

    return metro_result;
}

idx
swappy(
       int *atoms, 
       gsl_rng *r
      )
{   
    idx index;
    int idx_1 = (int)(1999 * gsl_rng_uniform(r)); // Random index for atom A
    int A = atoms[idx_1]; // Value of atom A
    int idx_2 = (int)(1999 * gsl_rng_uniform(r)); // Random index for atom B
    int B = atoms[idx_2]; // Value of atom B
    // Make sure that A and B are not the same atom type
    while (A == B)
    {
        idx_2 = (int)(1999 * gsl_rng_uniform(r));
        B = atoms[idx_2];
    }
    // Store the indices and values of the atoms and return the resulting struct
    index.alpha = idx_1;
    index.beta = idx_2;
    index.valueA = A;
    index.valueB = B;

    return index;
}

double
energy_bond(
            idx index, 
            int *atoms, 
            int **neighbors, 
            double E_cucu, 
            double E_znzn, 
            double E_cuzn
           )
{   
    double E = 0.; // Energy
    int atom_1_idx = index.alpha; // Index of atom A
    int atom_2_idx = index.beta; // Index of atom B
    int atom_1 = atoms[atom_1_idx]; // Value of atom A
    int atom_2 = atoms[atom_2_idx]; // Value of atom B
    int *neighbors_1 = neighbors[atom_1_idx]; // Neighbors of atom A
    int *neighbors_2 = neighbors[atom_2_idx]; // Neighbors of atom B
    // Loop over the neighbors of atom A
    for (int i = 0; i < 8; i++)
    {   
        // Add the energy of the bond between atom A and its neighbors
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
    // Loop over the neighbors of atom B
    for (int i = 0; i < 8; i++)
    {
        // Add the energy of the bond between atom B and its neighbors
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

atom_count
lattice_props(
              int *atoms, 
              int **neighbors, 
              int N_atoms
             )
{
    atom_count count;
    int N_Cu_A = 0; // Number of Cu atoms in sublattice A
    int N_CuZn = 0; // Number of Cu-Zn bonds
    // Loop over all atoms in sublattice A
    for (int i = 0; i < N_atoms / 2; i++)
    {   
        // If the atom is Cu, increment N_Cu_A
        if (atoms[i] == 1)
        {
            N_Cu_A++;
        }
        int atom = atoms[i]; // Value of the atom
        // Loop over the neighbors of the atom
        for (int j = 0; j < 8; j++)
        {
            // If the atom is Cu and the neighbor is Zn or vice versa, increment N_CuZn
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
    // Store the values in the struct and return it
    count.N_Cu_A = N_Cu_A;
    count.N_CuZn = N_CuZn;

    return count;
}

void
lattice_to_files(
                 FILE *fp_atoms, 
                 FILE *fp_neighbors, 
                 int *atoms, 
                 int **neighbors, 
                 int N_atoms
                )
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
