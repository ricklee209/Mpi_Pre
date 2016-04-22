



#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

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
)

{

#include "BCM.h"
#include "prm.h"

	printf("ncube\tnx\tn_wallcube\n");
	printf("%d\t%d\t%d\n",ncube,nx,n_wallcube);

	int Ncell = NcubeX*NcubeY*NcubeZ;

	int *jkl_switch;
	jkl_switch = new int [Ncell+1];
	
	int ii, i1, i2, pp, qq, n_change, flag, ch_pos, count_index;


	char str[1024];
	int BCM_ncube, BCM_n_wallcube, BCM_n;

	FILE *fptr;
	fptr = fopen("BCM_Grid.dat","r");

	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);

// ======================================================= //
	fscanf(fptr,"%d\n",&BCM_ncube);    //  how many cubes  //
// ======================================================= //

	fscanf(fptr,"%[^\n]\n",str);

// ======================================================================================================= //
	for (icube = 1; icube < ncube; icube++) {                                                              //
		                                                                                                   //
		fscanf(fptr,"%lf\t%lf\t%lf\t%lf\n",&cube_size[icube],&Xcube[icube],&Ycube[icube],&Zcube[icube]);   //
		                                                                                                   //
	}                                                                                                      //
// ======================================================================================================= //

	fscanf(fptr,"%[^\n]\n",str);

// =================================================================================== //
	for (icube = 1; icube < ncube; icube++) {                                          //
		for (int direction_index = 1; direction_index <= 6;  direction_index++) {      //
			                                                                           //
			fscanf(fptr,"%d\t%d\t%d\t%d\n",&adj_number[icube][1][direction_index],     //
				                           &adj_number[icube][2][direction_index],     //
										   &adj_number[icube][3][direction_index],     //
										   &adj_number[icube][4][direction_index]);    //
		}                                                                              //
	}                                                                                  //
// =================================================================================== //

	fscanf(fptr,"%[^\n]\n",str);

// ============================================================= //
	fscanf(fptr,"%d\n",&BCM_n);    //  how many cells in a cube  //
// ============================================================= //

	fscanf(fptr,"%[^\n]\n",str);

// ================================================================ //
	fscanf(fptr,"%d\n",&BCM_n_wallcube);    //  how many wallcubes  //
// ================================================================ //

	fscanf(fptr,"%[^\n]\n",str);

// =============================================================== //
	for (iwallcube = 1; iwallcube < n_wallcube; iwallcube++) {     //
		                                                           //
		fscanf(fptr,"%d\n",&wallcube[iwallcube]);                  //
		                                                           //
	}                                                              //
// =============================================================== //

	fscanf(fptr,"%[^\n]\n",str);
	



double size_min = MAX;

// ===================================================================================== //
	for (icube = 1; icube < ncube; icube++) {                                            //
		                                                                                 //
		size_min = min(size_min,cube_size[icube]);                                       //
		                                                                                 //
	}																					 //
	                                                                                     //
	for (icube = 1; icube < ncube; icube++) {                                            //
																					     //
	  csl[icube] = static_cast<int>(0.5+log( cube_size[icube]/size_min ) / log(2.0));    //
																						 //
	  //-- Plus 0.5 is for preventing truncation error --//                              //
																						 //
	}																					 //
		                                                                                 //
// ===================================================================================== //


	fclose(fptr);

	printf("%d\t%d\t%d\n\n",BCM_ncube,BCM_n,BCM_n_wallcube);

}