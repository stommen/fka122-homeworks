#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Struct to store the results of the DMC simulation
typedef struct
{
    double E_T; // The total energy of the system
    double E_T_avg; // The average total energy of the system
    int N_sprinters; // The number of sprinters
    double **R; // The positions of the walkers
} DMC_results;


/* ********************************************************************************** 
*
* DMC - Function to perform a Diffusion Monte Carlo simulation
*
* Parameters
* ----------
* N_its - the number of iterations
* N_0 - the initial number of walkers
* N_sprinters - the instantaneous number of walkers
* dtau - the time step
* E_T_start - the initial total energy of the system
* U - the random number generator
* gamma - the branching parameter
* fp_ET - the file pointer to the file to store the total energy
* fp_w - the file pointer to the file to store the walkers
* ndim - the number of dimensions
* IS - whether to use importance sampling
* decomp - the decomposition parameter
*
* Returns
* -------
* The results of the DMC simulation
*
* ******************************************************************************** */
DMC_results DMC(
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
               );

/* **********************************************************************************
*
* morse - Function to calculate the Morse potential
*
* Parameters
* ----------
* x - the distance between the two particles
*
* Returns
* -------
* The Morse potential
*
* ******************************************************************************** */
double morse(
             double x
            );

/* **********************************************************************************
*
* weight - Function to calculate the weight of a walker
*
* Parameters
* ----------
* dtau - the time step
* E_T - the total energy of the system
* V - the potential energy of the system
*
* Returns
* -------
* The weight of the walker
*
* ******************************************************************************** */
double weight(
              double dtau, 
              double E_T, 
              double V
             ); 

/* **********************************************************************************
*
* init_gsl_rng - Function to initialize the random number generator
*
* Parameters
* ----------
* seed - the seed for the random number generator
*
* Returns
* -------
* The random number generator
*
* ******************************************************************************** */
gsl_rng * init_gsl_rng(
                       int seed
                      ); 

/* **********************************************************************************
*
* init_walkmen_1D - Function to initialize the walkers in 1D
*
* Parameters
* ----------
* walkers - the array to store the walkers
* N_walkers - the number of walkers
*
* ******************************************************************************** */
void init_walkmen_1D(
                     double **walkers, 
                     int N_walkers
                    );

/* **********************************************************************************
*
* init_walkmen_6D - Function to initialize the walkers in 6D
*
* Parameters
* ----------
* walkers - the array to store the walkers
* N_walkers - the number of walkers
* U - the random number generator
*
* ******************************************************************************** */
void init_walkmen_6D(
                     double **walkers,
                     int N_walkers,
                     gsl_rng *U
                    );

/* **********************************************************************************
*
* local_energy - Function to calculate the local energy of the system
*
* Parameters
* ----------
* R - the positions of the walkers
* alpha - the variational parameter
*
* Returns
* -------
* The local energy of the system
*
* ******************************************************************************** */
double local_energy(
                    double *R,
                    double alpha
                   );

/* **********************************************************************************
*
* update_positions - Function to update the positions of the walkers
*
* Parameters
* ----------
* walkmen_pos - the positions of the walkers
* N_survivers - the number of walkers
* ndim - the number of dimensions
* dtau - the time step
* fp_w - the file pointer to the file to store the walkers
* U - the random number generator
* IS - whether to use importance sampling
* decomp - the decomposition parameter
*
* ******************************************************************************** */
void update_positions(
                      double **walkmen_pos,
                      int N_survivers,
                      int ndim,
                      double dtau,
                      FILE *fp_w,
                      gsl_rng *U,
                      bool IS,
                      int decomp
                     );

/* **********************************************************************************
*
* drift_velocity - Function to calculate the drift velocity of the walkers
*
* Parameters
* ----------
* R - the positions of the walkers
* alpha - the variational parameter
* dtau - the time step
* decomp - the decomposition parameter
* U - the random number generator
*
* ******************************************************************************** */
void drift_velocity(
                    double *R,
                    double alpha,
                    double dtau,
                    int decomp,
                    gsl_rng *U
                   );
