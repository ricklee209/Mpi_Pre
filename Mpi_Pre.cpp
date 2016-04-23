#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>


#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

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
)
{
	#include "BCM.h"
	
	
	int n_edge = 0;

	char file_name[100];


	for (icube = 1; icube < ncube; icube++) {  

		for (i = 1; i <= 6; i++ ) {

			if (adj_number[icube][1][i] > 0 |
				adj_number[icube][2][i] > 0 |
				adj_number[icube][3][i] > 0 |
				adj_number[icube][4][i] > 0) {

					if (adj_number[icube][1][i] != adj_number[icube][2][i] &
						adj_number[icube][1][i] != adj_number[icube][3][i] &
						adj_number[icube][1][i] != adj_number[icube][4][i]) {
							
							
							if(adj_number[icube][1][i] > icube) n_edge=n_edge+1;
							if(adj_number[icube][2][i] > icube) n_edge=n_edge+1;
							if(adj_number[icube][3][i] > icube) n_edge=n_edge+1;
							if(adj_number[icube][4][i] > icube) n_edge=n_edge+1;
							

							//n_edge=n_edge+4;

						}

					else if(adj_number[icube][1][i] > icube) n_edge=n_edge+1;
						//n_edge = n_edge+1;

			}

		}    // ---- for (i = 1; i <= 6; i++ ) ---- //

	}// ---- for (icube = 1; icube < ncube; icube++) ---- //

	printf("n_edge = %d\n",n_edge);


	FILE *fptr;

	fptr = fopen("metisfn.dat","wb"); 

	fprintf(fptr,"%d\t%d\n",ncube-1,n_edge);    // ---- Metis format on line 1 ---- //

	int iedge = 0;

	
	for (icube = 1; icube < ncube; icube++) {  

		for (i = 1; i <= 6; i++ ) {

			if (adj_number[icube][1][i] > 0 |
				adj_number[icube][2][i] > 0 |
				adj_number[icube][3][i] > 0 |
				adj_number[icube][4][i] > 0) {

					if (adj_number[icube][1][i] != adj_number[icube][2][i] &
						adj_number[icube][1][i] != adj_number[icube][3][i] &
						adj_number[icube][1][i] != adj_number[icube][4][i]) {

						  fprintf(fptr,"%d\t%d\t%d\t%d\t",abs(adj_number[icube][1][i]),abs(adj_number[icube][2][i]),abs(adj_number[icube][3][i]),abs(adj_number[icube][4][i]));

						}

					else {

						  fprintf(fptr,"%d\t",abs(adj_number[icube][1][i]));

						}

			}

		}    // ---- for (i = 1; i <= 6; i++ ) ---- //

		fprintf(fptr,"\n");

	}// ---- for (icube = 1; icube < ncube; icube++) ---- //

	fclose(fptr);    // ---- close metis file ---- //


	//sprintf(file_name,"gpmetis.exe metisfn.dat ""%d",np);

	//system(file_name);    // ---- creating the metis file ---- //




// ---- reading the metis file to build rank_map ---- //

// ---- rank_map [Cubes(Global)] [Ranks] [Cubes(local)] ---- //

// =========================================================== //
															   
	//sprintf(file_name,"metisfn.dat.part.""%d",np);             
	//
	////sprintf(file_name,"metisfn_self.dat");        
	//			
	//fptr = fopen(file_name,"r");                               
	//                                                          
	//for (icube = 1; icube < ncube; icube++) {
	//
	//	fscanf(fptr,"%d\n",&rank_map[0][icube]);

	//	//printf("%d\t",rank_map[icube]);
	//
	//	}

	//fclose(fptr);    // ---- close metis file ---- //


	int Xdv = int( (ncube-1)/np);
	int Xr =  (ncube-1)-Xdv*np;

	for (int irank = 0; irank < np; irank++) {

		if(irank < Xr) {

			for (int i_mts = irank*(Xdv+1)+1; i_mts < irank*(Xdv+1)+Xdv+2; i_mts++ )
				rank_map[0][i_mts] = (np-1)-irank;

		}
		else {

			for (int i_mts =irank*Xdv+Xr+1; i_mts < irank*Xdv+Xr+Xdv+1; i_mts++ )
			rank_map[0][i_mts] = (np-1)-irank;

		}
	}





	int rank_cubes = 0;
	int rank_wall_cubes = 0;

	for (int irank = 0; irank < np; irank++) {

		rank_cubes = 0;

		for (icube = 1; icube < ncube; icube++) {


			if (irank == rank_map[0][icube]) {
				
				rank_cubes = rank_cubes+1;

				rank_map[1][icube] = rank_cubes;

				if(rank_cubes == 1) cube_map[0][irank] = icube;

				}



		}    // ---- for (icube = 1; icube < ncube; icube++) ---- //

		cube_map[1][irank] = cube_map[0][irank]+rank_cubes;

	}    // ---- for (int irank = 0; irank < np; irank++) ---- //

	
	sprintf(file_name,".\\MPI_files\\rank_map.dat");             
															   
	fptr = fopen(file_name,"wb");   

	fprintf(fptr,"Cube_N(global)\t CPU_rank\t Cube_N(local)\n");
	                                                       
	for (icube = 1; icube < ncube; icube++) {

		fprintf(fptr,"%d\t%d\t%d\n",icube,rank_map[0][icube],rank_map[1][icube]);

		}
		
	fclose(fptr);    // ---- close rank_map.dat ---- //
	

// =========================================================== //	




	FILE *fptr_runlength;
	fptr_runlength = fopen("BCM_Grid.dat","r");

	char runlength[1024] = "JKL_runlength for cell information >>";
	char str[1024];

	
	rewind(fptr_runlength);
	
	fscanf(fptr_runlength,"%[^\n]\n",str);

	while( strcmp(str,runlength) != 0)  { 
		fscanf(fptr_runlength,"%[^\n]\n",str); 
	}



	int n_change, flag, ch_pos;

	int count_index = 0;

	int Ncell = NcubeX*NcubeY*NcubeZ;

	int adj;

	int (*tmp_runlength) = new int[Ncell+1];


	#pragma omp parallel
	for (int i = 1; i <= Ncell; i++) {

		tmp_runlength[i] = 0;

	}

	/*
	
	#pragma omp parallel for private(\
	fptr,icube,rank_cubes,adj,\
	iwallcube,str,n_change,flag,ch_pos,rank_wall_cubes,fptr_runlength\
	)
	
	*/

	for (int irank = np-1; irank >= 0; irank--) {


		sprintf(file_name,".\\MPI_files\\BCMgrid""%0.5d"".dat",irank);    
		fptr = fopen(file_name,"wb"); 

		fprintf(fptr,"BCM GridGene Output Data Version 1.4 - Single Run-Length\n");
		fprintf(fptr,"#_of_cubes in the computational region >>\n");


// =================================== Cube size and corner information =================================== //

		rank_cubes = 0;

		
		fprintf(fptr,"%d\n",cube_map[1][irank]-cube_map[0][irank]);    // ---- cubes number for each CPU ---- //



		fprintf(fptr,"size_&_left-corner_location (x,y,z) of cube(i) >>\n");

		for (icube = cube_map[0][irank]; icube < cube_map[1][irank]; icube++) {

			fprintf(fptr,"%.9f\t%.9f\t%.9f\t%.9f\n",cube_size[icube],Xcube[icube],Ycube[icube],Zcube[icube]);

		}    // ---- for (icube = 1; icube < ncube; icube++) ---- //

// =================================== Cube size and corner information =================================== //



// =================================== Cube neighbor information =================================== //

		fprintf(fptr,"Cube_neighboring_information >>\n");

		i = 0;
		for (icube = cube_map[0][irank]; icube < cube_map[1][irank]; icube++) {

			for (int direction_index = 1; direction_index <= 6;  direction_index++) {    



				for (int iadj = 1; iadj <= 4; iadj++) {

					if (irank == rank_map[0][abs(adj_number[icube][iadj][direction_index])]) {

						if ( adj_number[icube][1][direction_index] >= 0 )

							adj = rank_map[1][adj_number[icube][iadj][direction_index]];

						else 
							adj = -rank_map[1][-adj_number[icube][iadj][direction_index]];

						}

					else adj = 0;

					fprintf(fptr,"%d\t",adj);
						

					}    // ---- for (int iadj = 1; iadj <= 4; iadj++) ---- //


			}    // ---- for (int direction_index = 1; direction_index <= 6;  direction_index++) ---- //

			fprintf(fptr,"\n");

		}    // ---- for (icube = 1; icube < ncube; icube++) ---- //
		fprintf(fptr,"\n");

// =================================== Cube neighbor information =================================== //




		fprintf(fptr,"#_of_cells on a cube_edge >>\n");
		fprintf(fptr,"%d\n",NcubeX);




// =================================== Wall cube number =================================== //

		fprintf(fptr,"#_of_wall_cubes >>\n");

		rank_wall_cubes = 0;
		#pragma omp parallel for reduction(+:rank_wall_cubes) 
		for (icube = cube_map[0][irank]; icube < cube_map[1][irank]; icube++) {

			for (iwallcube = 1; iwallcube < n_wallcube; iwallcube++) {

				if (icube == wallcube[iwallcube]) 

					rank_wall_cubes = rank_wall_cubes+1;

			}

		}    // ---- for (icube = 1; icube < ncube; icube++) ---- //
		#pragma omp barrier

		fprintf(fptr,"%d\n",rank_wall_cubes);

// =================================== Wall cube number =================================== //




		fprintf(fptr,"Wall_cube_# >>\n");

		for (icube = cube_map[0][irank]; icube < cube_map[1][irank]; icube++) {

			for (iwallcube = 1; iwallcube < n_wallcube; iwallcube++) {

				if (icube == wallcube[iwallcube]) 

					fprintf(fptr,"%d\t",rank_map[1][wallcube[iwallcube]]);

			}

		}    // ---- for (icube = 1; icube < ncube; icube++) ---- //
		fprintf(fptr,"\n");


		

// =================================== Run length information =================================== //

		fprintf(fptr,"JKL_runlength for cell information >>\n");

		

		for (iwallcube = 1; iwallcube < n_wallcube; iwallcube++) {


			if (irank == rank_map[0][wallcube[iwallcube]]) {

				fscanf(fptr_runlength,"%d\n",&n_change);    // ---- how many cells needed to be changed the flag ---- //

				fscanf(fptr_runlength,"%d\n",&flag);    // ---- the first cell's flag ---- //


				for (int iflag = 1; iflag <= n_change; iflag++) {

					fscanf(fptr_runlength,"%d\n",&tmp_runlength[iflag]);

				}

				fprintf(fptr,"%d\t",n_change);
				fprintf(fptr,"%d\t",flag);


				for (int iflag = 1; iflag <= n_change; iflag++) {

					fprintf(fptr,"%d\t",tmp_runlength[iflag]);

					tmp_runlength[iflag] = 0;

				}


				fprintf(fptr,"\n");

			}    // ---- if (icube == wallcube[iwallcube]) { ---- //

		}

		//fprintf(fptr_runlength,"\n");

// =================================== Run length information =================================== //

		
		fclose(fptr);    // ---- close BCMgrid file ---- //    

	}    // ---- for (int irank = 0; irank < np; irank++) ---- //

	
	#pragma omp barrier


	fclose(fptr_runlength);    // ---- close original BCMgrid file ---- 


}

