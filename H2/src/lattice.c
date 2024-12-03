#include <stdio.h>

void init_sc(double **positions, int N, double lattice_param, double base[1][3], double **neighbors, int N_neighbors, FILE *f_neighbors)
{   
    int ind = 0;
    double help[1][3] = {0, 10, 100};
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            for (int z = 0; z < N; z++)
            {
                positions[ind][0] = lattice_param * (base[0][0] + x);
                positions[ind][1] = lattice_param * (base[0][1] + y);
                positions[ind][2] = lattice_param * (base[0][2] + z);
                ind += 1;
		    }
        }
    }
    double help_plz[1][8] = {0, 1, 10, 11, 100, 101, 110, 111};
    for (int j = 0; j < N_neighbors; j++)
    {
        for (int i = 0; i < N*N*N; i++)
        {
            neighbors[i][j] = help_plz[1][j] + i*10;
            fprintf(f_neighbors, "%f\n", neighbors[i][j]);
        } 
    }
}

