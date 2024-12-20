#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double Mooooooorse(double x);
void wavef_ground_state_anal(double *Psi,double *x, int n);
double weight(double dtau,double ET, double x);


int
run(
    int argc,
    char *argv[]
   )
{
    int N_sprinters = 300;
    double E0 = 3./8.;
    double dtau = 0.02;
    double E_T = 0.5;



    return 0;
}


double
Mooooooorse(double x)
{
    double V = 0.5 * pow((1  - exp(-x)),2);
    
    return V;
}

double
weight(double dtau,double ET, double x)
{
    double W = exp(-(Mooooooorse(x) - ET) * dtau);
    return W;
}

void
wavef_ground_state_anal(double *Psi,double *x, int n)
{
    for(int i = 0;i < n; i++)
    {
        Psi[i] = sqrt(2.) * exp(-exp(-x[i]) - x[i]/2.);
    }
}

