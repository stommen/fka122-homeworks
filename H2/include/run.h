#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Struct for storing the indices and values of the atoms to swap
typedef struct
{
    int alpha; // Index of atom A
    int beta; // Index of atom B
    int valueA; // Value of atom A
    int valueB; // Value of atom B
} idx;

// Struct for storing the results of the metropolis algorithm
typedef struct
{
    int accepted; // Number of accepted swaps
    double Etot; // Total energy at the end of the MC iterations
    int its; // Number of iterations (task 3)
} metro;

// Struct for storing the number of Cu atoms in sublattice A and the number of Cu-Zn bonds
typedef struct
{
    int N_Cu_A; // Number of Cu atoms in sublattice A
    int N_CuZn; // Number of Cu-Zn bonds in sublattice A
} atom_count;

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
gsl_rng * init_gsl_rng(int seed);

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
* Nearest_neighbors_bcc - Function to find the nearest neighbors in a BCC lattice
*
* Parameters
* ----------
* pos_A - the positions of the atoms in sublattice A
* pos_B - the positions of the atoms in sublattice B
* N_atoms - the number of atoms in the supercell
* neighbors - the neighbors of the atoms (matrix)
* box_length - the length of the supercell
* cutoff - the cutoff from the closest distance to consider as a neighbor
* closest_distance - the closest distance between atoms
*
* ******************************************************************************** */
void nearest_neighbors_bcc(double **pos_A, double **pos_B, int N_atoms,
                           int **neighbors, double box_length,
                           double cutoff, double closest_distance);

/* **********************************************************************************
*
* Metropolis - Function to perform the Metropolis algorithm
*
* Parameters
* ----------
* its - the number of iterations (task 2) or the number of accepted swaps (task 3)
* atoms - the atoms in the supercell
* neighbors - the neighbors of the atoms (matrix)
* k_B - the Boltzmann constant [eV/K]
* T - the temperature [K]
* E_cucu - the energy of a Cu-Cu bond [eV]
* E_znzn - the energy of a Zn-Zn bond [eV]
* E_cuzn - the energy of a Cu-Zn bond [eV]
* E_tot - the total energy of the system [eV]
* r - random number generator
* C_V - the heat capacity [eV/K]
* U - the energy of the system [eV]
* P - the long-range order parameter
* R - the short-range order parameter
* N_atoms - the number of atoms in the supercell
* fp - the data file pointer
*
* Returns
* -------
* The results of the Metropolis algorithm (number of accepted swaps, total energy, and number of iterations)
*
* ******************************************************************************** */
metro metropolis(int its, int *atoms, int **neighbors, double k_B, double T, double E_cucu, 
                 double E_znzn, double E_cuzn, double E_tot, gsl_rng *r, 
                 double *C_V, double *U, double *P, double *R, int N_atoms, FILE *fp);

/* **********************************************************************************
*
* Swappy - Function to swap two atoms in the supercell
*
* Parameters
* ----------
* atoms - the atoms in the supercell
* r - random number generator
*
* Returns
* -------
* The indices and values of the atoms to swap
*
* ******************************************************************************** */
idx swappy(int *atoms, gsl_rng *r);

/* **********************************************************************************
*
* Energy_bond - Function to calculate the energy after swapping two atoms
*
* Parameters
* ----------
* index - the indices and values of the atoms to swap
* atoms - the atoms in the supercell
* neighbors - the neighbors of the atoms (matrix)
* E_cucu - the energy of a Cu-Cu bond [eV]
* E_znzn - the energy of a Zn-Zn bond [eV]
* E_cuzn - the energy of a Cu-Zn bond [eV]
*
* Returns
* -------
* The energy of the system after swapping two atoms
*
* ******************************************************************************** */
double energy_bond(idx index, int *atoms, int **neighbors,
                   double E_cucu, double E_znzn, double E_cuzn);

/* **********************************************************************************
*
* Lattice_props - Function to calculate the number of Cu atoms in sublattice A and the number of Cu-Zn bonds
*
* Parameters
* ----------
* atoms - the atoms in the supercell
* neighbors - the neighbors of the atoms (matrix)
* N_atoms - the number of atoms in the supercell
*
* Returns
* -------
* The number of Cu atoms and the number of Cu-Zn bonds in sublattice A
*
* ******************************************************************************** */
atom_count lattice_props(int *atoms, int **neighbors, int N_atoms);

/* **********************************************************************************
*
* Lattice_to_files - Function to write the atoms and neighbors to files
*
* Parameters
* ----------
* fp_atoms - the file pointer for the atoms
* fp_neighbors - the file pointer for the neighbors
* atoms - the atoms in the supercell
* neighbors - the neighbors of the atoms (matrix)
* N_atoms - the number of atoms in the supercell
*
* ******************************************************************************** */
void lattice_to_files(FILE *fp_atoms, FILE *fp_neighbors, int *atoms, int **neighbors, int N_atoms);
