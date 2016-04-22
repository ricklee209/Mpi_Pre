


#include "Resolution.h"

double (*cube_size_pre) = new double[ncube];
int (*csl_pre) = new int[ncube];

double (*Xcube_pre) = new double[ncube];
double (*Ycube_pre) = new double[ncube];
double (*Zcube_pre) = new double[ncube];

int (*adj_number_pre)[5][7] = new int[ncube][5][7];
int (*wallcube_pre) = new int[n_wallcube];

int (*rank_map)[ncube] = new int[2][ncube];    // ---- [Cubes(Global)] [Ranks] [Cubes(local)] ---- // 

int (*cube_map)[np] = new int[2][np];

