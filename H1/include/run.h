#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void calculate(double *potential, double *virial, double **force,
           double **position, double cell_length, int nbr_atoms);

/* **********************************************************************************   
*
* Init_gsl_rng - Function to initialize the random number generator
*
* Parameters
* ----------
* seed - the seed for the random number generator
*
* Returns
* -------
* A pointer to the random number generator
*
* ******************************************************************************** */   
gsl_rng *init_gsl_rng(int seed);

/* **********************************************************************************
*
* Verlet - Function to perform the Verlet algorithm and calculate the properties
*          of the system
*
* Parameters
* ----------
* positions - the positions of the atoms, N x 3 matrix
* velocities - the velocities of the atoms, N x 3 matrix
* forces - the forces on the atoms, N x 3 matrix
* its - the number of iterations
* its_eq - the number of iterations for equilibration
* delta_t - the time step
* m - the mass of the atoms
* k_B - the Boltzmann constant
* beta - the inverse of the bulk modulus
* a_0 - the lattice constant
* N - the number of atoms
* cell_length - the number of unit cells in each direction
* T_eq - the equilibrium temperature
* P_eq - the equilibrium pressure
* fp - the file pointer for the output file
* fp_traj - the file pointer for the trajectory file
* fp_rdist - the file pointer for the radial distribution file
* fp_sfact - the file pointer for the structure factor file
*
* Returns
* -------
* The lattice constant
*
* ******************************************************************************** */
double verlet(double **positions, double **velocities, double **forces, int its, int its_eq,
            double delta_t, double m, double k_B, double a_0, double beta, int N, 
            int cell_length, double T_eq, double P_eq, FILE *fp, FILE *fp_traj, 
            FILE *fp_rdist, FILE *fp_sfact);

/* **********************************************************************************
*
* Radial_dist - Function to calculate the radial distribution of the atoms
*
* Parameters
* ----------
* bins - the bins for the radial distribution
* positions - the positions of the atoms
* N_bins - the number of bins
* bin_width - the width of the bins
* N - the number of atoms
* L - the length of the supercell
*
* ******************************************************************************** */
void radial_dist(double *bins, double **positions, int N_bins, double bin_width, int N, double L);

/* **********************************************************************************
*
* Boundary_distance_between_vectors - Function to calculate the distance between two
*                                     vectors with periodic boundary conditions
*
* Parameters
* ----------
* v1 - the first vector
* v2 - the second vector
* dim - the dimension of the vectors
* box_length - the length of the supercell
*
* Returns
* -------
* The distance between the vectors
*
* ******************************************************************************** */
double boundary_distance_between_vectors(double *v1, double *v2, int dim, double box_length);

/* **********************************************************************************
*
* Init_grid - Function to initialize the grid for the structure factor
*
* Parameters
* ----------
* N_max - the maximum number of grid points
* L - the length of the supercell
*
* Returns
* -------
* The grid, a 3D grid, where every point is a 3D vector
*
* ******************************************************************************** */
double ****init_grid(int N_max, double L);

/* **********************************************************************************
*
* Destroy_grid - Function to destroy the grid for the structure factor
*
* Parameters
* ----------
* grid - the grid to be destroyed
* N_max - the maximum number of grid points
*
* ******************************************************************************** */
void destroy_grid(double ****grid, int N_max);

/* **********************************************************************************
*
* Structure_factor - Function to calculate the structure factor of the atoms
*
* Parameters
* ----------
* grid - the grid for the structure factor
* positions - the positions of the atoms
* Sq - the structure factor
* N - the number of atoms
* n - the number of grid points
*
* ******************************************************************************** */
void structure_factor(double ****grid, double **positions, double *Sq, int N, int n);

/* **********************************************************************************
*
* Spherical_avg - Function to calculate the spherical average of the structure factor
*
* Parameters
* ----------
* grid - the grid for the structure factor
* Sq - the structure factor
* n - the number of grid points
* n_bins - the number of bins
* bin_width - the width of the bins
* fp_sfact - the file pointer for the structure factor file
*
* ******************************************************************************** */
void spherical_avg(double ****grid, double *Sq, int n, int n_bins, double bin_width, FILE *fp_sfact);
