



#include "Resolution.h"

extern int NZbc_l;
extern int NZbc_u;


void BCM_Grid_Pre
(
// =================================================== //
double (*cube_size) = new double[ncube],
int (*csl) = new int[ncube],

double (*Xcube) = new double[ncube],
double (*Ycube) = new double[ncube],
double (*Zcube) = new double[ncube],

int (*adj_number)[5][7] = new int[ncube][5][7],

int (*wallcube) = new int[ncube]
// ================================================== //
);

void Boundary_detection 
(
// =================================================== //
int *NXbc_l,
int *NXbc_u,

int *NYbc_l,
int *NYbc_u,

int *NZbc_l,
int *NZbc_u,

double *Xmax,
double *Xmin,

double *Ymax,
double *Ymin,

double *Zmax,
double *Zmin,

double (*cube_size) = new double[ncube],
int (*csl) = new int[ncube],

double (*Xcube) = new double[ncube],
double (*Ycube) = new double[ncube],
double (*Zcube) = new double[ncube],

int (*adj_number)[5][7] = new int[ncube][5][7]

// =================================================== //
);

void Z_periodic_condition
(
// =================================================== //
int nZbc_l,
int nZbc_u,

double Zmin,
double Zmax,

int (*Zbc_u) = new int[NZbc_u+1],
int (*Zbc_l) = new int[NZbc_l+1],


double (*cube_size) = new double[ncube],

double (*Xcube) = new double[ncube],
double (*Ycube) = new double[ncube],
double (*Zcube) = new double[ncube],

int (*adj_number)[5][7] = new int[ncube][5][7]
// =================================================== //
);

void Mpi_Pre
(
// ============================================================================ //

double (*cube_size) = new double[ncube],
int (*csl) = new int[ncube],

double (*Xcube) = new double[ncube],
double (*Ycube) = new double[ncube],
double (*Zcube) = new double[ncube],

int (*adj_number)[5][7] = new int[ncube][5][7],

int (*wallcube) = new int[ncube],

int (*rank_map)[ncube] = new int[2][ncube],

int (*cube_map)[np] = new int[2][np]

// ============================================================================ //
);



void Communication_table
(
// ============================================================================ //

double (*cube_size) = new double[ncube],
int (*csl) = new int[ncube],

double (*Xcube) = new double[ncube],
double (*Ycube) = new double[ncube],
double (*Zcube) = new double[ncube],

int (*adj_number)[5][7] = new int[ncube][5][7],

int (*wallcube) = new int[ncube],

int (*rank_map)[ncube] = new int[2][ncube]
// ============================================================================ //
);

