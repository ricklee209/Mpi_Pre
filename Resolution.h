



#define np 1568	           // ---- how many CPUs ---- //

#define ncube 150496          // ---- cube number + 1 ---- //

#define n_wallcube 57886      // ---- wall cube number + 1 ---- // 

#define n_buffer 2         // ---- buffer region for BCM data communication (EVEN number) ---- //

#define NcubeX 16    
#define NcubeY 16   
#define NcubeZ 16    

#define nx NcubeX + n_buffer - 1    // ---- nx =17 ---- //
#define ny NcubeY + n_buffer - 1    
#define nz NcubeZ + n_buffer - 1    

#define X_size 20    // ---- NcubeX + 2 * n_buffer ---- //
#define Y_size 20    // ---- NcubeY + 2 * n_buffer ---- //
#define Z_size 20    // ---- NcubeZ + 2 * n_buffer ---- //

#define nxx nx+1
#define nyy ny+1
#define nzz nz+1


#define nxxx nxx+1
#define nyyy nyy+1
#define nzzz nzz+1 
