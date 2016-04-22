



#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <direct.h>

#include "main.h"
#include "prm.h"

int main(int argc, char **argv)
{   


	
	int NBC;

	int NXbc_l, NXbc_u, NYbc_l = 0, NYbc_u = 0, NZbc_l, NZbc_u;
	double Xmax, Xmin, Ymax, Ymin, Zmax, Zmin;
 
	double (*er) = new double[5];
	
	#include "BCM.h"
	#include "Resolution.h"
	#include "Array_Pre.h"    /**** allocate the memory ****/

	
	_mkdir(".\\MPI_files");

	
	BCM_Grid_Pre(cube_size_pre, csl_pre, Xcube_pre, Ycube_pre, Zcube_pre, adj_number_pre, wallcube_pre);


	Boundary_detection(&NXbc_l, &NXbc_u, &NYbc_l, &NYbc_u, &NZbc_l, &NZbc_u,
					   &Xmax, &Xmin, &Ymax, &Ymin, &Zmax, &Zmin,
		cube_size_pre, csl_pre, Xcube_pre, Ycube_pre, Zcube_pre, adj_number_pre);

	//printf("%f\t%f\t%f\t%f\t%f\t%f\n",Xmax,Xmin,Ymax,Ymin,Zmax,Zmin);

	int (*Xbc_u) = new int[NXbc_u+1]; 
	int (*Xbc_l) = new int[NXbc_l+1]; 

	int (*Ybc_u) = new int[NYbc_u+1]; 
	int (*Ybc_l) = new int[NYbc_l+1]; 

	int (*Zbc_u) = new int[NZbc_u+1]; 
	int (*Zbc_l) = new int[NZbc_l+1]; 


	//Z_periodic_condition(NZbc_l, NZbc_u, Zmin, Zmax, Zbc_u, Zbc_l, cube_size_pre, Xcube_pre, Ycube_pre, Zcube_pre, adj_number_pre);


	Mpi_Pre(cube_size_pre, csl_pre, Xcube_pre, Ycube_pre, Zcube_pre, adj_number_pre, wallcube_pre, rank_map, cube_map);

	Communication_table(cube_size_pre, csl_pre, Xcube_pre, Ycube_pre, Zcube_pre, adj_number_pre, wallcube_pre, rank_map);

	system("pause");


}