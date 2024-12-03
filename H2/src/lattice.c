#include <stdio.h>

void init_sc(double **positions, int N, double lattice_param, double base[3], 
             int **neighbors, int N_neighbors, FILE *f_neighbors)
{   
    int ind = 0;
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            for (int z = 0; z < N; z++)
            {
                positions[ind][0] = lattice_param * (base[0] + x);
                positions[ind][1] = lattice_param * (base[1] + y);
                positions[ind][2] = lattice_param * (base[2] + z);
                ind += 1;
		    }
        }
    }
    double help_plz[8] = {0, 1, 10, 11, 100, 101, 110, 111};
    for (int j = 0; j < N_neighbors; j++)
    {
        for (int i = 0; i < N*N*N; i++)
        {
            neighbors[i][j] = (int)(help_plz[j] + i);
        } 
    }
    for (int i = 0; i < N*N*N; i++)
    {
        for (int j = 0; j < N_neighbors; j++)
        {
            fprintf(f_neighbors, "%i, ", neighbors[i][j]);
        }
        fprintf(f_neighbors, "\n");
    }
}

