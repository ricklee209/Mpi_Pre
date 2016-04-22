#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <omp.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"


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
)
{

	#include "BCM.h"
	#include "prm.h"

	
	*Xmax = MIN;    /**** initialize ****/
	*Xmin = MAX;    /**** initialize ****/

	*Ymax = MIN;    /**** initialize ****/
	*Ymin = MAX;    /**** initialize ****/

	*Zmax = MIN;    /**** initialize ****/
	*Zmin = MAX;    /**** initialize ****/

	for (icube = 1; icube < ncube; icube++) {  

		*Xmin = min(*Xmin,Xcube[icube]);
		*Xmax = max(*Xmax,Xcube[icube]+cube_size[icube]);

		*Ymin = min(*Ymin,Ycube[icube]);
		*Ymax = max(*Ymax,Ycube[icube]+cube_size[icube]);

		*Zmin = min(*Zmin,Zcube[icube]);
		*Zmax = max(*Zmax,Zcube[icube]+cube_size[icube]);
		
	}

	*NXbc_l = *NXbc_u = *NYbc_l = *NYbc_u = *NZbc_l = *NZbc_u = 0;

	for (icube = 1; icube < ncube; icube++) {  

		if (fabs(Xcube[icube]+cube_size[icube]-*Xmax) < minimum) *NXbc_u = *NXbc_u+1;
		if (fabs(Ycube[icube]+cube_size[icube]-*Ymax) < minimum) *NYbc_u = *NYbc_u+1;
		if (fabs(Zcube[icube]+cube_size[icube]-*Zmax) < minimum) *NZbc_u = *NZbc_u+1;

		if (fabs(Xcube[icube]-*Xmin) < minimum) *NXbc_l = *NXbc_l+1;
		if (fabs(Ycube[icube]-*Ymin) < minimum) *NYbc_l = *NYbc_l+1;
		if (fabs(Zcube[icube]-*Zmin) < minimum) *NZbc_l = *NZbc_l+1;
		
	};

	//printf("%d\t%d\t%d\t%d\t%d\t%d\n",*NXbc_u,*NXbc_l,*NYbc_u,*NYbc_l,*NZbc_u,*NZbc_l);




}