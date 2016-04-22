#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int NZbc_l;
extern int NZbc_u;

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
)
{

	#include "BCM.h"
	#include "prm.h"

	int adj;

	int lcube, ucube;

	//printf("%d\t%d\n",nZbc_l,nZbc_u);
	
	i = 0;
	for (icube = 1; icube < ncube; icube++) {  
	
		if (fabs(Zcube[icube]-Zmin) < minimum) {  i = i+1; Zbc_l[i] = icube; }
		
	};

	i = 0;

	for (icube = 1; icube < ncube; icube++) {  

		if (fabs(Zcube[icube]+cube_size[icube]-Zmax) < minimum) { i = i+1; Zbc_u[i] = icube; }

	};
	

	
	for (icube = 1; icube <= nZbc_l; icube++) {  

		lcube = Zbc_l[icube];
		
		for (iicube = 1; iicube <= nZbc_u; iicube++) {  

			ucube = Zbc_u[iicube];

			if (fabs(Xcube[lcube]-Xcube[ucube]) < minimum && fabs(Ycube[lcube]-Ycube[ucube]) < minimum && lcube != ucube) {

					//printf("%f\t%f\t",Zcube[lcube],Zcube[ucube]+cube_size[ucube]);

					//printf("%d\t%d\t",adj_number[lcube][1][5],adj_number[ucube][1][6]);

					adj_number[lcube][1][5] = ucube;
					adj_number[lcube][2][5] = ucube;
					adj_number[lcube][3][5] = ucube;
					adj_number[lcube][4][5] = ucube;

					adj_number[ucube][1][6] = lcube;
					adj_number[ucube][2][6] = lcube;
					adj_number[ucube][3][6] = lcube;
					adj_number[ucube][4][6] = lcube;

			}
				
				
				
		}
		

	}			
	



}