#include <omp.h>
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <omp.h>



struct Prop
{

	int CPU;
	int Cube;
	int CPU_neig;
	int Cube_neig;
	int Dir;
	int Adj;
	int Inter;
	int i_sort;
};

struct Array
{
	Prop *content;
	unsigned int size;
};


int compare_Inter( const void *a, const void *b )
{

	return ( (Prop*)a)->Inter > ( (Prop*)b)->Inter ?1:-1;

};

int compare_CPU( const void *a, const void *b )
{

	return ( (Prop*)a)->i_sort > ( (Prop*)b)->i_sort ?1:-1;

};


#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

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
)
{
	#include "BCM.h"
	
	
	int n_edge = 0;

	char file_name[100];

	FILE *fptr;
	sprintf(file_name,"MPI_communication_table.dat");    
	fptr = fopen(file_name,"wb"); 

	//fprintf(fptr,"Cube\t CPU\t Cube_neighbor\t CPU_neighbor\t direction_N\t Interface_N\n");

	fprintf(fptr,"Cube\t CPU\t Cube_neig\t CPU_neig\t Dir\t Adj\t Inter\n");




	int rank_cubes = 0;
	int rank_wall_cubes = 0;
	int irank, index;
	
	int acube;    // ---- adjacent cube ---- //

	int acube1, acube2, acube3, acube4;    // ---- adjacent cube ---- //

	int iiacube;

	int count_index = 0;

	for (icube = 1; icube < ncube; icube++) { 

		irank = rank_map[0][icube];

// ======================================== [X] direction ======================================== //


		if (adj_number[icube][1][1] > 0 |
			adj_number[icube][2][1] > 0 |
			adj_number[icube][3][1] > 0 |
			adj_number[icube][4][1] > 0) {

// ----------------------------------------- X equal ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][1];

					if (irank != rank_map[0][acube] & csl[icube] == csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],1,0);

						count_index = count_index+1;

						fprintf(fptr,"0\n");

						break;

					}

				}

// ----------------------------------------- X equal ----------------------------------------- //


// ----------------------------------------- X B 2 S ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][1];

					if (irank != rank_map[0][acube] & csl[icube] > csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],1,iacube);
						
						count_index = count_index+1;

						fprintf(fptr,"1\n");
					}

				}

				
// ----------------------------------------- X B 2 S ----------------------------------------- //


// ----------------------------------------- X S 2 B ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][1];

					if (irank != rank_map[0][acube] & csl[icube] < csl[acube]) {

						//if (iacube < 4 & rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][1]]) continue;

						if (iacube < 4)	{ if( rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][1]] ) continue; }    // -- the same neighbor cube -- //
						
						for (iiacube = 1; iiacube <= 4; iiacube++ ) if (icube == adj_number[acube][iiacube][2]) break;    // -- ???? -- //

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],1,iiacube);

						count_index = count_index+1;

						fprintf(fptr,"-1\n");
					}

				}

// ----------------------------------------- X S 2 B ----------------------------------------- //
		}



		if (adj_number[icube][1][2] > 0 |
			adj_number[icube][2][2] > 0 |
			adj_number[icube][3][2] > 0 |
			adj_number[icube][4][2] > 0) {

// ----------------------------------------- X equal ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][2];

					if (irank != rank_map[0][acube] & csl[icube] == csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],2,0);

						count_index = count_index+1;

						fprintf(fptr,"0\n");

						break;

					}

				}

// ----------------------------------------- X equal ----------------------------------------- //


// ----------------------------------------- X B 2 S ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][2];

					if (irank != rank_map[0][acube] & csl[icube] > csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],2,iacube);
						
						count_index = count_index+1;

						fprintf(fptr,"1\n");
					}

				}
				
// ----------------------------------------- X B 2 S ----------------------------------------- //


// ----------------------------------------- X S 2 B ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][2];

					if (irank != rank_map[0][acube] & csl[icube] < csl[acube]) {
						
						//if (iacube < 4 & rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][2]]) continue;

						if (iacube < 4)	{ if( rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][2]] ) continue; }
						

						for (iiacube = 1; iiacube <= 4; iiacube++ ) if (icube == adj_number[acube][iiacube][1]) break;
						
						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],2,iiacube);

						count_index = count_index+1;

						fprintf(fptr,"-1\n");
					}

				}

// ----------------------------------------- X S 2 B ----------------------------------------- //

		}

// ======================================== [X] direction ======================================== //




// ======================================== [Y] direction ======================================== //

		if (adj_number[icube][1][3] > 0 |
			adj_number[icube][2][3] > 0 |
			adj_number[icube][3][3] > 0 |
			adj_number[icube][4][3] > 0) {

// ----------------------------------------- Y equal ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][3];

					if (irank != rank_map[0][acube] & csl[icube] == csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],3,0);
						
						count_index = count_index+1;

						fprintf(fptr,"0\n");

						break;

					}

				}

// ----------------------------------------- Y equal ----------------------------------------- //


// ----------------------------------------- Y B 2 S ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][3];

					if (irank != rank_map[0][acube] & csl[icube] > csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],3,iacube);
						
						count_index = count_index+1;

						fprintf(fptr,"1\n");
					}

				}


// ----------------------------------------- Y B 2 S ----------------------------------------- //


// ----------------------------------------- Y S 2 B ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][3];

					if (irank != rank_map[0][acube] & csl[icube] < csl[acube]) {

						//if (iacube < 4 & rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][3]]) continue;

						if (iacube < 4)	{ if( rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][3]] ) continue; }
						
						for (iiacube = 1; iiacube <= 4; iiacube++ ) if (icube == adj_number[acube][iiacube][4]) break;
						
						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],3,iiacube);
						
						count_index = count_index+1;

						fprintf(fptr,"-1\n");
					}

				}

// ----------------------------------------- Y S 2 B ----------------------------------------- //
		}


		if (adj_number[icube][1][4] > 0 |
			adj_number[icube][2][4] > 0 |
			adj_number[icube][3][4] > 0 |
			adj_number[icube][4][4] > 0) {

// ----------------------------------------- Y equal ----------------------------------------- //


				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][4];

					if (irank != rank_map[0][acube] & csl[icube] == csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],4,0);
						
						count_index = count_index+1;

						fprintf(fptr,"0\n");

						break;

					}

				}

// ----------------------------------------- Y equal ----------------------------------------- //


// ----------------------------------------- Y B 2 S ----------------------------------------- //

				
				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][4];

					if (irank != rank_map[0][acube] & csl[icube] > csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],4,iacube);
						
						count_index = count_index+1;

						fprintf(fptr,"1\n");
					}

				}

// ----------------------------------------- Y B 2 S ----------------------------------------- //


// ----------------------------------------- Y S 2 B ----------------------------------------- //

				/*
				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][4];

					if (irank != rank_map[0][acube] & csl[icube] < csl[acube]) {

						printf("%d\t",rank_map[0][adj_number[icube][iacube][4]]);


					}

				}
				*/

				
				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][4];

					if (irank != rank_map[0][acube] & csl[icube] < csl[acube]) {

						//if (iacube < 4 & rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][4]]) continue;

						if (iacube < 4)	{ if( rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][4]] ) continue; }
						
							
						for (iiacube = 1; iiacube <= 4; iiacube++ ) if (icube == adj_number[acube][iiacube][3]) break;

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],4,iiacube);

						count_index = count_index+1;

						fprintf(fptr,"-1\n");
					}

				}


// ----------------------------------------- Y S 2 B ----------------------------------------- //

		}

// ======================================== [Y] direction ======================================== //




// ======================================== [Z] direction ======================================== //



		
		if (adj_number[icube][1][5] > 0 |
			adj_number[icube][2][5] > 0 |
			adj_number[icube][3][5] > 0 |
			adj_number[icube][4][5] > 0) {

// ----------------------------------------- Z equal ----------------------------------------- //


				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][5];

					if (irank != rank_map[0][acube] & csl[icube] == csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],5,0);
						
						count_index = count_index+1;

						fprintf(fptr,"0\n");

						break;

					}

				}

// ----------------------------------------- Z equal ----------------------------------------- //

// ----------------------------------------- Z B 2 S ----------------------------------------- //

				

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][5];

					if (irank != rank_map[0][acube] & csl[icube] > csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],5,iacube);
						
						count_index = count_index+1;

						fprintf(fptr,"1\n");
					}

				}

// ----------------------------------------- Z B 2 S ----------------------------------------- //


// ----------------------------------------- Z S 2 B ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][5];

					if (irank != rank_map[0][acube] & csl[icube] < csl[acube]) {
						
						//if (iacube < 4 & rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][5]]) continue;

						if (iacube < 4)	{ if( rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][5]] ) continue; }
						

						for (iiacube = 1; iiacube <= 4; iiacube++ ) if (icube == adj_number[acube][iiacube][6]) break;
						
						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],5,iiacube);

						count_index = count_index+1;

						fprintf(fptr,"-1\n");
					}

				}

// ----------------------------------------- Z S 2 B ----------------------------------------- //

		}


		if (adj_number[icube][1][6] > 0 |
			adj_number[icube][2][6] > 0 |
			adj_number[icube][3][6] > 0 |
			adj_number[icube][4][6] > 0) {

// ----------------------------------------- Z equal ----------------------------------------- //


				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][6];

					if (irank != rank_map[0][acube] & csl[icube] == csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],6,0);
						
						count_index = count_index+1;

						fprintf(fptr,"0\n");

						break;

					}

				}

// ----------------------------------------- Z equal ----------------------------------------- //

// ----------------------------------------- Z B 2 S ----------------------------------------- //

				

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][6];

					if (irank != rank_map[0][acube] & csl[icube] > csl[acube]) {

						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],6,iacube);
						
						count_index = count_index+1;

						fprintf(fptr,"1\n");
					}

				}

// ----------------------------------------- Z B 2 S ----------------------------------------- //


// ----------------------------------------- Z S 2 B ----------------------------------------- //

				for (int iacube = 1; iacube <= 4; iacube++ ) {

					acube = adj_number[icube][iacube][6];

					if (irank != rank_map[0][acube] & csl[icube] < csl[acube]) {
						
						//if (iacube < 4 & rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][6]]) continue;

						if (iacube < 4)	{ if( rank_map[0][acube] == rank_map[0][adj_number[icube][iacube+1][6]] ) continue; }
						
						for (iiacube = 1; iiacube <= 4; iiacube++ ) if (icube == adj_number[acube][iiacube][5]) break;
						
						fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t",rank_map[1][icube],irank,rank_map[1][acube],rank_map[0][acube],6,iiacube);

						count_index = count_index+1;

						fprintf(fptr,"-1\n");
					}

				}

// ----------------------------------------- Z S 2 B ----------------------------------------- //

		}



// ======================================== [Z] direction ======================================== //



	}    // ---- for (icube = 1; icube < ncube; icube++) ---- //

	fprintf(fptr,"Number of Cubes for connection  >>\n");
	fprintf(fptr,"%d\n",count_index);

	fclose(fptr);    // ---- close communication_table.dat ---- //




	char str[1024];
	int nadj;
	int MPI_Nadj = count_index+1;


	Array MPI_table;

	MPI_table.size = count_index;

	MPI_table.content = new Prop[MPI_table.size];

	int (*MPI_cube_) = new int[MPI_Nadj];
	int (*MPI_cpu_) = new int[MPI_Nadj];

	int (*MPI_cube_adj_) = new int[MPI_Nadj];
	int (*MPI_cpu_adj_) = new int[MPI_Nadj];

	int (*MPI_direction_) = new int[MPI_Nadj];
	int (*MPI_interface_) = new int[MPI_Nadj];

	int (*MPI_adjn_) = new int[MPI_Nadj];


	int (*MPI_cube) = new int[MPI_Nadj];
	int (*MPI_cpu) = new int[MPI_Nadj];

	int (*MPI_cube_adj) = new int[MPI_Nadj];
	int (*MPI_cpu_adj) = new int[MPI_Nadj];

	int (*MPI_direction) = new int[MPI_Nadj];
	int (*MPI_interface) = new int[MPI_Nadj];

	int (*MPI_adjn) = new int[MPI_Nadj];

	
	int index0 = 0;
	int index_0 = 0;
	int index0_ = 0;

	int iindex0 = 0;
	int iindex_0 = 0;
	int iindex0_ = 0;

	int count_CPU = 0;

	int t1,t2,t3,t4,t5,t6,t7;



	fptr = fopen("MPI_communication_table.dat","r"); 

	fscanf(fptr,"%[^\n]\n",str);

	for (nadj = 0; nadj < MPI_Nadj-1; nadj++) {
		
		fscanf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t%d\n", &t1, &t2, &t3, &t4, &t5, &t6, &t7);

		MPI_table.content[nadj].Cube = t1;
		MPI_table.content[nadj].CPU = t2;
		MPI_table.content[nadj].Cube_neig = t3;
		MPI_table.content[nadj].CPU_neig = t4;
		MPI_table.content[nadj].Dir = t5;
		MPI_table.content[nadj].Adj = t6;
		MPI_table.content[nadj].Inter = t7;

	}
	
	fclose(fptr);


	
	
	for (nadj = 0; nadj < MPI_Nadj-1; nadj++) {

		
		t1 = MPI_table.content[nadj].Cube;
		t2 = MPI_table.content[nadj].CPU;
		t4 = MPI_table.content[nadj].CPU_neig;
		t7 = MPI_table.content[nadj].Inter;

		
		if(t7 == 1)  MPI_table.content[nadj].i_sort = ( t1 + (t4+1)*np + (t2+1)*np*np ) + 4*np*np*np;
		if(t7 == 0)  MPI_table.content[nadj].i_sort = ( t1 + (t4+1)*np + (t2+1)*np*np );
		if(t7 == -1) MPI_table.content[nadj].i_sort = ( t1 + (t4+1)*np + (t2+1)*np*np ) + 2*np*np*np;
		
	}

	qsort( MPI_table.content, MPI_table.size, sizeof( Prop ), compare_CPU );
	


	#pragma omp parallel for reduction(+:index0, index0_) 

	for (nadj = 1; nadj < MPI_Nadj; nadj++) {

		if(MPI_table.content[nadj-1].Inter == 0) index0 = index0 + 1;
		if(MPI_table.content[nadj-1].Inter == 1) index0_ = index0_ + 1;

	}


	index_0 = MPI_Nadj-index0-index0_;

	for (nadj = 1; nadj < MPI_Nadj; nadj++) {


		MPI_cpu_[nadj] = MPI_table.content[nadj-1].CPU;
		     MPI_cube_[nadj] = MPI_table.content[nadj-1].Cube;
		     MPI_cpu_adj_[nadj] = MPI_table.content[nadj-1].CPU_neig;
		     MPI_cube_adj_[nadj] = MPI_table.content[nadj-1].Cube_neig;
		     MPI_direction_[nadj] = MPI_table.content[nadj-1].Dir;
		     MPI_adjn_[nadj] = MPI_table.content[nadj-1].Adj;
		     MPI_interface_[nadj] = MPI_table.content[nadj-1].Inter;
			 

		/*
		if(MPI_table.content[nadj].Inter == 100) {

			 iindex0 = iindex0 + 1;

			 MPI_cpu_[iindex0] = MPI_table.content[nadj-1].CPU;
		     MPI_cube_[iindex0] = MPI_table.content[nadj-1].Cube;
		     MPI_cpu_adj_[iindex0] = MPI_table.content[nadj-1].CPU_neig;
		     MPI_cube_adj_[iindex0] = MPI_table.content[nadj-1].Cube_neig;
		     MPI_direction_[iindex0] = MPI_table.content[nadj-1].Dir;
		     MPI_adjn_[iindex0] = MPI_table.content[nadj-1].Adj;
		     MPI_interface_[iindex0] = MPI_table.content[nadj-1].Inter;
			 
		}
		
		
		if(MPI_table.content[nadj].Inter == 100) {
			
			iindex0 = iindex0 + 1;

			 MPI_cpu_[iindex0] = MPI_table.content[nadj-1].CPU;
		     MPI_cube_[iindex0] = MPI_table.content[nadj-1].Cube;
		     MPI_cpu_adj_[iindex0] = MPI_table.content[nadj-1].CPU_neig;
		     MPI_cube_adj_[iindex0] = MPI_table.content[nadj-1].Cube_neig;
		     MPI_direction_[iindex0] = MPI_table.content[nadj-1].Dir;
		     MPI_adjn_[iindex0] = MPI_table.content[nadj-1].Adj;
		     MPI_interface_[iindex0] = MPI_table.content[nadj-1].Inter;

		    }

		
		if(MPI_table.content[nadj].Inter == 10) {
			
			iindex_0 = iindex_0 + 1;

			 MPI_cpu_[index0 + iindex_0] = MPI_table.content[nadj-1].CPU;
		     MPI_cube_[index0 + iindex_0] = MPI_table.content[nadj-1].Cube;
		     MPI_cpu_adj_[index0 + iindex_0] = MPI_table.content[nadj-1].CPU_neig;
		     MPI_cube_adj_[index0 + iindex_0] = MPI_table.content[nadj-1].Cube_neig;
		     MPI_direction_[index0 + iindex_0] = MPI_table.content[nadj-1].Dir;
		     MPI_adjn_[index0 + iindex_0] = MPI_table.content[nadj-1].Adj;
		     MPI_interface_[index0 + iindex_0] = MPI_table.content[nadj-1].Inter;

		    }

		
		
		if(MPI_table.content[nadj].Inter == 1) {
			
			iindex0_ = iindex0_ + 1;

			 MPI_cpu_[index0+index_0+iindex0_] = MPI_table.content[nadj-1].CPU;
		     MPI_cube_[index0+index_0+iindex0_] = MPI_table.content[nadj-1].Cube;
		     MPI_cpu_adj_[index0+index_0+iindex0_] = MPI_table.content[nadj-1].CPU_neig;
		     MPI_cube_adj_[index0+index_0+iindex0_] = MPI_table.content[nadj-1].Cube_neig;
		     MPI_direction_[index0+index_0+iindex0_] = MPI_table.content[nadj-1].Dir;
		     MPI_adjn_[index0+index_0+iindex0_] = MPI_table.content[nadj-1].Adj;
		     MPI_interface_[index0+index_0+iindex0_] = MPI_table.content[nadj-1].Inter;

		    }
			*/
	}




	sprintf(file_name,".\\MPI_files\\MPI_communication_table_in_turn.dat");    
	fptr = fopen(file_name,"wb"); 

	fprintf(fptr,"Number of Cubes for connection  >>\n");
	fprintf(fptr,"%d\n",count_index);

	fprintf(fptr,"CPU\t Cube\t CPU_neig\t Cube_neig\t Dir\t Adj\t Inter\n");


	for (nadj = 1; nadj < MPI_Nadj; nadj++) {

		t1 = MPI_cpu_[nadj];
		t2 = MPI_cube_[nadj];
		t3 = MPI_cpu_adj_[nadj];
		t4 = MPI_cube_adj_[nadj];
		t5 = MPI_direction_[nadj];
		t6 = MPI_adjn_[nadj];
		t7 = MPI_interface_[nadj];
					

		if(t7 == 100) t7 = 0;
		if(t7 == 10) t7 = -1;


		fprintf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t%d\n",t1, t2, t3, t4, t5, t6, t7);

	}

	fclose(fptr);



 // ============================== Reading the MPI_communication_table_in_turn ============================== //

	sprintf(file_name,".\\MPI_files\\MPI_communication_table_in_turn.dat");    
	fptr = fopen(file_name,"r"); 

	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);
	fscanf(fptr,"%[^\n]\n",str);

	for (nadj = 1; nadj < MPI_Nadj; nadj++) {

		fscanf(fptr,"%d\t%d\t%d\t%d\t%d\t%d\t%d\n", &MPI_cpu[nadj], &MPI_cube[nadj], &MPI_cpu_adj[nadj], &MPI_cube_adj[nadj], &MPI_direction[nadj], &MPI_adjn[nadj], &MPI_interface[nadj]);

	}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

	fclose(fptr);

 // ============================== Reading the MPI_communication_table_in_turn ============================== //




// ============================== communication table for each CPU ============================== //


	
	int indicator = 0;

	for (int irank = 0; irank < np; irank ++) {

		sprintf(file_name,".\\MPI_files\\Communication_irank.""%05d"".dat",irank);             

		fptr = fopen(file_name,"wb"); 




// -------------------------- The same size adjancent -------------------------- //

		fprintf(fptr,"# The same size adjancent\n");
		fprintf(fptr,"0\n");

		count_CPU = 1;

		indicator = 0;
		for (nadj = 1; nadj < index0+1; nadj++) {

			//if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 

			if (irank == MPI_cpu[nadj]) {

				indicator = 1;

				if (MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1] & MPI_cpu[nadj] == MPI_cpu[nadj+1]) 

				count_CPU = count_CPU + 1;

			}

			//if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1] & MPI_cpu[nadj] == MPI_cpu[nadj+1])  
			//    printf("%d\t%d\n",MPI_cpu_adj[nadj],MPI_cpu_adj[nadj+1]);

		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"# Numbers of neighbor CPUs\n");

		if (indicator == 1) fprintf(fptr,"%d\n", count_CPU);
		else fprintf(fptr,"%d\n", 0);

		//printf("%d\n",count_CPU);



		fprintf(fptr,"# Neighbor CPU rank \n");

		indicator = 0;
		for (nadj = 1; nadj < index0+1; nadj++) {

			//if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 

			if (irank == MPI_cpu[nadj]) {

				indicator = 1;

				if (MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1] & MPI_cpu[nadj] == MPI_cpu[nadj+1]) {

					fprintf(fptr,"%d\t", MPI_cpu_adj[nadj]);

				}

				i = nadj;

			}    // ---- if (irank == MPI_cpu[nadj]) ---- //

		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		if (indicator == 1) fprintf(fptr,"%d\t", MPI_cpu_adj[i]);    // ---- write MPI_cpu_adj[nadj+1] ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Numbers of cubes for each CPU \n");
		
		count_index = 1;
		index = 0;
		for (nadj = 1; nadj < index0+1; nadj++) {
			if (irank == MPI_cpu[nadj]) {

				if (irank == MPI_cpu[nadj-1] && irank == MPI_cpu[nadj+1]) {

					if (MPI_cpu_adj[nadj] == MPI_cpu_adj[nadj+1]) 

						count_index = count_index+1;

					else  { 
						fprintf(fptr,"%d\t", count_index);	

						count_index = 1;

					}

				}
				

				// ---- The first term for each Current CPU ---- //
				//if (irank != MPI_cpu[nadj-1]) {

				if (irank != MPI_cpu[nadj-1] && irank == MPI_cpu[nadj+1]) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;

					}
					else
						count_index = 2;

				}


				// ---- The last term for each Current CPU ---- //
				if (irank != MPI_cpu[nadj+1] && irank == MPI_cpu[nadj-1] ) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj-1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;
					
					}
				}


				if (irank != MPI_cpu[nadj+1] && irank != MPI_cpu[nadj-1] ) 

					fprintf(fptr,"%d\t", 1);	


			}    // ---- if (irank == MPI_cpu[nadj]) ---- //

		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		if (count_index > 1) fprintf(fptr,"%d\t", count_index); 

		//if (index == 0) fprintf(fptr,"%d\t", count_index);	

		fprintf(fptr,"\n");





		fprintf(fptr,"# Current CPU Cube number for sending \n");

		for (nadj = 1; nadj < index0+1; nadj++) {
			
			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_cube[nadj]);

			
			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Current CPU Cube number for recving \n");

		for (nadj = 1; nadj < index0+1; nadj++) {

			/*
			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_cube_adj[nadj]);

			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");
				*/

			
			if (irank == MPI_cpu_adj[nadj]) 
				fprintf(fptr,"%d\t", MPI_cube_adj[nadj]);

			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Direction index in current CPU for sending \n");

		for (nadj = 1; nadj < index0+1; nadj++) {

			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_direction[nadj]);

			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");



		fprintf(fptr,"# Direction index in current CPU for recving \n");

		for (nadj = 1; nadj < index0+1; nadj++) {

			if (irank == MPI_cpu_adj[nadj]) 
				fprintf(fptr,"%d\t", MPI_direction[nadj]);
				
			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");



		

// -------------------------- The same size adjancent -------------------------- //		


fprintf(fptr,"\n\n");

		
// --------------------------  Small to big -------------------------- //

		fprintf(fptr,"# Small to big\n");
		fprintf(fptr,"-1\n");

		count_CPU = 1;

		indicator = 0;
		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			if (irank == MPI_cpu[nadj]) {

				indicator = 1;

				if (MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1] & MPI_cpu[nadj] == MPI_cpu[nadj+1]) 

					count_CPU = count_CPU + 1;

			}
			
			
		}    // ---- for (nadj = index0+1; nadj < index0+index_0+1; nadj++) ---- //

		fprintf(fptr,"# Numbers of neighbor CPUs\n");

		
		if (indicator == 1) fprintf(fptr,"%d\n", count_CPU);
		else fprintf(fptr,"%d\n", 0);





		fprintf(fptr,"# Neighbor CPU rank \n");

		indicator = 0;
		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			
			//printf("%d\t%d\n",irank,MPI_cpu_adj[nadj]);


			if (irank == MPI_cpu[nadj]) {

				indicator = 1;

				if (MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1] & MPI_cpu[nadj] == MPI_cpu[nadj+1]) {

					fprintf(fptr,"%d\t", MPI_cpu_adj[nadj]);

				}

				i = nadj;

			}    // ---- if (irank == MPI_cpu[nadj]) ---- //

		}    // ---- for (nadj = index0+1; nadj < index0+index_0+1; nadj++) ---- //
		
		if (indicator == 1) fprintf(fptr,"%d\t", MPI_cpu_adj[i]);    // ---- write MPI_cpu_adj[nadj+1] ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Numbers of cubes for each CPU \n");
		
		count_index = 1;
		index = 0;
		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			if (irank == MPI_cpu[nadj]) {

				if (irank == MPI_cpu[nadj-1] && irank == MPI_cpu[nadj+1]) {

					if (MPI_cpu_adj[nadj] == MPI_cpu_adj[nadj+1]) 

						count_index = count_index+1;

					else  { 
						fprintf(fptr,"%d\t", count_index);	

						count_index = 1;

						index = 1;
					}

				}

				/*if (irank != MPI_cpu[nadj-1]) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;

					}

					else 
						count_index = 2;

				}

				if (irank != MPI_cpu[nadj+1] && irank == MPI_cpu[nadj-1] ) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj-1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;
					
					}
				}*/

				// ---- The first term for each Current CPU ---- //
				//if (irank != MPI_cpu[nadj-1]) {

				if (irank != MPI_cpu[nadj-1] && irank == MPI_cpu[nadj+1]) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;

					}
					else
						count_index = 2;

				}


				// ---- The last term for each Current CPU ---- //
				if (irank != MPI_cpu[nadj+1] && irank == MPI_cpu[nadj-1] ) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj-1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;
					
					}
				}


				if (irank != MPI_cpu[nadj+1] && irank != MPI_cpu[nadj-1] ) 

					fprintf(fptr,"%d\t", 1);	




			}    // ---- if (irank == MPI_cpu[nadj]) ---- //

		}    // ---- for (nadj = index0+1; nadj < index0+index_0+1; nadj++) ---- //

		
		if (count_index > 1) fprintf(fptr,"%d\t", count_index); 


		//if (index == 0) fprintf(fptr,"%d\t", count_index);	

		fprintf(fptr,"\n");




		fprintf(fptr,"# Current CPU Cube number for sending \n");

		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_cube[nadj]);

			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");



		fprintf(fptr,"#  Current CPU Cube number for recving \n");

		//for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			/*
			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_cube_adj[nadj]);

			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");
			*/

			if (irank == MPI_cpu_adj[nadj]) 
				fprintf(fptr,"%d\t", MPI_cube_adj[nadj]);

			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Direction index in current CPU for sending \n");

		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_direction[nadj]);

			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Direction index in current CPU for recving \n");

		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			if (irank == MPI_cpu_adj[nadj]) 
				fprintf(fptr,"%d\t", MPI_direction[nadj]);

			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Adjacent order of current cube(small) for adjacent cube(big) for recving \n");

		//for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			if (irank == MPI_cpu_adj[nadj]) 
				fprintf(fptr,"%d\t", MPI_adjn[nadj]);

			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");





// --------------------------  Small to big -------------------------- //		


		fprintf(fptr,"\n\n");

		
// --------------------------  Big to small -------------------------- //

		fprintf(fptr,"# Big to small \n");
		fprintf(fptr,"1\n");

		count_CPU = 1;

		indicator = 0;

		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			if (irank == MPI_cpu[nadj]) {

				indicator = 1;

				if (MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1] & MPI_cpu[nadj] == MPI_cpu[nadj+1]) 

				count_CPU = count_CPU + 1;

			}

		}    // ---- for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) ---- //

		fprintf(fptr,"# Numbers of neighbor CPUs\n");

		if (indicator == 1) fprintf(fptr,"%d\n", count_CPU);
		else fprintf(fptr,"%d\n", 0);


		

		fprintf(fptr,"# Neighbor CPU rank \n");

		indicator = 0;
		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			if (irank == MPI_cpu[nadj]) {

				indicator = 1;

				if (MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1] & MPI_cpu[nadj] == MPI_cpu[nadj+1]) {

					fprintf(fptr,"%d\t", MPI_cpu_adj[nadj]);

				}

				i = nadj;

			}    // ---- if (irank == MPI_cpu[nadj]) ---- //

		}    // ---- for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) ---- //
			
		if (indicator == 1) fprintf(fptr,"%d\t", MPI_cpu_adj[i]);    // ---- write MPI_cpu_adj[nadj+1] ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Numbers of cubes for each CPU \n");
		
		count_index = 1;
		index = 0;

		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {


			if (irank == MPI_cpu[nadj]) {

				if (irank == MPI_cpu[nadj-1] && irank == MPI_cpu[nadj+1]) {

					if (MPI_cpu_adj[nadj] == MPI_cpu_adj[nadj+1]) 

						count_index = count_index+1;

					else  { 
						fprintf(fptr,"%d\t", count_index);	

						count_index = 1;

						index = 1;
					}

				}

				/*if (irank != MPI_cpu[nadj-1]) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;

					}

					else 
						count_index = 2;

				}

				if (irank != MPI_cpu[nadj+1] && irank == MPI_cpu[nadj-1] ) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj-1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;
					
					}
				}*/

				// ---- The first term for each Current CPU ---- //
				//if (irank != MPI_cpu[nadj-1]) {

				if (irank != MPI_cpu[nadj-1] && irank == MPI_cpu[nadj+1]) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;

					}
					else
						count_index = 2;

				}


				// ---- The last term for each Current CPU ---- //
				if (irank != MPI_cpu[nadj+1] && irank == MPI_cpu[nadj-1] ) {

					if ((MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj-1])) {
						
						fprintf(fptr,"%d\t", 1);	
						index = 1;
					
					}
				}


				if (irank != MPI_cpu[nadj+1] && irank != MPI_cpu[nadj-1] ) 

					fprintf(fptr,"%d\t", 1);	


			}    // ---- if (irank == MPI_cpu[nadj]) ---- //

		}    // ---- for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) ---- //
		
		if (count_index > 1) fprintf(fptr,"%d\t", count_index); 

		//if (index == 0) fprintf(fptr,"%d\t", count_index);	


		fprintf(fptr,"\n");




		fprintf(fptr,"# Current CPU Cube number for sending \n");

		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_cube[nadj]);

			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Current CPU Cube number for recving \n");

		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			if (irank == MPI_cpu_adj[nadj])  {

				fprintf(fptr,"%d\t", MPI_cube_adj[nadj]);
				
			}


			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Direction index in current CPU for sending \n");

		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			if (irank == MPI_cpu[nadj]) 
				fprintf(fptr,"%d\t", MPI_direction[nadj]);

			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");


		fprintf(fptr,"# Direction index in current CPU for recving \n");

		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			if (irank == MPI_cpu_adj[nadj]) 
				fprintf(fptr,"%d\t", MPI_direction[nadj]);


			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Adjacent order of adjacent cube(small) for current cube(big) for recving \n");

		for (nadj = index0+1; nadj < index0+index_0+1; nadj++) {

			if (irank == MPI_cpu_adj[nadj]) 
				fprintf(fptr,"%d\t", MPI_adjn[nadj]);


			if (irank == MPI_cpu_adj[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");




		fprintf(fptr,"# Adjacent order of adjacent cube(small) for current cube(big) for Sending \n");

		
		for (nadj = index0+index_0+1; nadj < MPI_Nadj; nadj++) {

			if (irank == MPI_cpu[nadj])
				fprintf(fptr,"%d\t", MPI_adjn[nadj]);


			if (irank == MPI_cpu[nadj] & MPI_cpu_adj[nadj] != MPI_cpu_adj[nadj+1]) 
				fprintf(fptr,"\n");


		}    // ---- for (nadj = 1; nadj < index0+1; nadj++) ---- //

		fprintf(fptr,"\n");



// --------------------------  Big to small -------------------------- //



		fclose(fptr);

	}    // ---- for (int irank = 0; irank < np; irank ++) ---- //



// ============================== communication table for each CPU ============================== //

	
	







	//printf("%d\n",index0+index_+_index);

	//printf("%d\n",_index);

	delete []MPI_table.content;

	
	delete []MPI_cube_;
	delete []MPI_cpu_;
	delete []MPI_cube_adj_;
	delete []MPI_cpu_adj_;
	delete []MPI_direction_;
	delete []MPI_interface_;
	

	delete []MPI_cube;
	delete []MPI_cpu;
	delete []MPI_cube_adj;
	delete []MPI_cpu_adj;
	delete []MPI_direction;
	delete []MPI_interface;



}